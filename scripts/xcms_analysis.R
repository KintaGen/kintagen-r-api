# Suppress package startup messages
suppressPackageStartupMessages({
  library(jsonlite)
  library(Spectra)
  library(mzR)
  library(dplyr)
  library(tidyr)
  library(MsCoreUtils)
  library(xcms)
  library(MsExperiment) 
  library(ProtGenerics)
})

# --- HELPERS ---
source("./scripts/helpers.R")


# ==============================================================================
# LOGIC A: SINGLE FILE (Original Logic)
# ==============================================================================

empty_peaks_df <- function() {
  data.frame(rt_apex = numeric(), rt_start = numeric(), rt_end = numeric(), peak_height = numeric(), peak_area = numeric(), stringsAsFactors = FALSE)
}

detect_and_integrate_peaks_original <- function(rt, intensity, noise_thresh) {
    if (length(rt) != length(intensity)) {
        min_len <- min(length(rt), length(intensity))
        rt <- rt[1:min_len]; intensity <- intensity[1:min_len]
    }
    is_peak <- which(diff(sign(diff(intensity))) == -2) + 1
    apexes <- is_peak[intensity[is_peak] > noise_thresh]
    if (length(apexes) == 0) return(empty_peaks_df())
    
    get_valley <- function(idx, dir){
      j <- idx
      while (j > 1 && j < length(intensity) && ((dir < 0 && intensity[j-1] <= intensity[j]) || (dir > 0 && intensity[j+1] <= intensity[j]))) { j <- j + dir }
      j
    }
    
    peak_list <- lapply(apexes, function(apex_idx) {
      left_idx  <- get_valley(apex_idx, -1)
      right_idx <- get_valley(apex_idx,  1)
      if (left_idx >= right_idx) return(NULL)
      
      peak_indices <- left_idx:right_idx
      peak_rt <- rt[peak_indices]; peak_int <- intensity[peak_indices]
      if (length(peak_rt) < 2) return(NULL)
      peak_area <- sum(diff(peak_rt) * (peak_int[-1] + peak_int[-length(peak_int)]) / 2)
      
      data.frame(rt_apex = rt[apex_idx], rt_start = rt[left_idx], rt_end = rt[right_idx], peak_height = intensity[apex_idx], peak_area = peak_area, stringsAsFactors = FALSE)
    })
    valid_peaks <- do.call(rbind, peak_list)
    if (is.null(valid_peaks) || nrow(valid_peaks) == 0) return(empty_peaks_df())
    distinct(valid_peaks, rt_start, rt_end, .keep_all = TRUE)
}

process_single_file <- function(mzMl_file_path, noise_threshold_percent) {
    log_message("Single file detected. Using standard detection...")
    output <- list(results = list())
    
    exp_spectra <- Spectra(mzMl_file_path, backend = MsBackendMzR())
    raw_rt <- rtime(exp_spectra)
    raw_int <- tic(exp_spectra)
    
    sg_coeffs <- MsCoreUtils::coefSG(hws = 5L, k = 3L)
    smooth_int <- MsCoreUtils::smooth(raw_int, cf = sg_coeffs)
    
    max_val <- max(smooth_int, na.rm = TRUE)
    abs_thresh <- max_val * (noise_threshold_percent / 100)
    peaks <- detect_and_integrate_peaks_original(raw_rt, smooth_int, abs_thresh)
    
    if (nrow(peaks) > 0) {
        peaks <- peaks %>% arrange(rt_apex) %>% mutate(peak_number = row_number())
        total_area <- sum(peaks$peak_area, na.rm=TRUE)
        output$results$quantitative_report <- peaks %>%
            mutate(area_percent = if(total_area > 0) (peak_area/total_area)*100 else 0, rt_minutes = round(rt_apex/60, 3), peak_area = round(peak_area, 0), area_percent = round(area_percent, 2)) %>%
            select(peak_number, rt_minutes, peak_area, area_percent)
        output$results$integrated_peaks_details <- peaks %>% mutate(rt_apex_min = rt_apex/60, rt_start_min = rt_start/60, rt_end_min = rt_end/60)
            
        top_peaks <- peaks %>% arrange(desc(peak_height)) %>% head(50)
        output$results$top_spectra_data <- lapply(1:nrow(top_peaks), function(i) {
             p <- top_peaks[i,]
             idx <- which.min(abs(raw_rt - p$rt_apex))
             vals <- peaksData(exp_spectra[idx])[[1]]
             df <- data.frame(mz=numeric(), relative_intensity=numeric())
             if(!is.null(vals) && nrow(vals) > 0) {
                 max_i <- max(vals[,2], na.rm=TRUE)
                 if(max_i > 0) df <- data.frame(mz = vals[,1], intensity = vals[,2]) %>% mutate(relative_intensity = (intensity/max_i)*100) %>% filter(relative_intensity >= 1.0) %>% select(mz, relative_intensity)
             }
             list(peak_number = p$peak_number, spectrum_data = df)
        })
    } else {
        output$results$quantitative_report <- list()
        output$results$integrated_peaks_details <- list()
        output$results$top_spectra_data <- list()
    }
    
    min_len <- min(length(raw_rt), length(smooth_int))
    output$results$raw_chromatogram_data <- data.frame(rt_min = raw_rt[1:min_len]/60, intensity = raw_int[1:min_len])
    output$results$smoothed_chromatogram_data <- data.frame(rt_min = raw_rt[1:min_len]/60, intensity = smooth_int[1:min_len])
    return(output)
}

# ==============================================================================
# LOGIC B: MULTI-FILE (XCMS ALIGNMENT + MULTI-EXPORT)
# ==============================================================================

# ==============================================================================
# LOGIC B: MULTI-FILE (XCMS ALIGNMENT + MULTI-EXPORT)
# ==============================================================================

process_multi_file <- function(file_paths, noise_threshold_percent) {
    log_message("Multiple files detected. Using XCMS for alignment...")
    output <- list(results = list())
    
    # 1. SETUP EXPERIMENT
    pd <- data.frame(sample_name = paste0("File_", 1:length(file_paths)), sample_group = "Study")
    mse <- readMsExperiment(spectraFiles = file_paths, sampleData = pd)
    
    # Auto-detect MS Level (usually Level 1 for chromatograms)
    lvls <- unique(msLevel(spectra(mse)))
    target_lvl <- if(1L %in% lvls) 1L else 2L
    
    # Sort spectra by file and time
    sps <- spectra(mse)
    spectra(mse) <- sps[order(sps$dataStorage, rtime(sps))]

    # 2. EXTRACT "BEFORE" DATA (Raw BPC)
    # The documentation uses BPC (Base Peak Chromatogram, agg="max") for alignment plots, not TIC ("sum")
    log_message("Extracting Pre-Alignment Chromatograms...")
    
    # Create the BPC object from raw data
    bpc_raw_obj <- chromatogram(mse, aggregationFun = "max", msLevel = target_lvl)
    
    # Helper to extract dataframes from chromatogram objects
    extract_chrom_data <- function(chrom_obj) {
        lapply(1:length(file_paths), function(i) {
            chr <- chrom_obj[1, i]
            # Convert to minutes and filter NAs
            df <- data.frame(
                rt_min = rtime(chr)/60, 
                intensity = intensity(chr)
            ) %>% filter(!is.na(intensity))
            
            list(filename = basename(file_paths[i]), data = df)
        })
    }
    
    output$results$chromatograms_before_alignment <- extract_chrom_data(bpc_raw_obj)

    # 3. PEAK DETECTION
    # Calculate noise based on BPC max value
    max_val <- max(sapply(bpc_raw_obj, function(x) max(intensity(x), na.rm=TRUE)), na.rm=TRUE)
    abs_noise <- if(is.finite(max_val)) (max_val/5)*(noise_threshold_percent/100) else 1000
    
    cwp <- CentWaveParam(ppm=50, peakwidth=c(5, 60), snthresh=5, noise=abs_noise, prefilter=c(3, abs_noise))
    mse <- findChromPeaks(mse, param=cwp, msLevel=target_lvl)
    
    # Fallback if strict parameters fail
    if(nrow(chromPeaks(mse)) == 0) {
        cwp_loose <- CentWaveParam(ppm=50, peakwidth=c(5,60), snthresh=3, noise=0, prefilter=c(3, 100))
        mse <- findChromPeaks(mse, param=cwp_loose, msLevel=target_lvl)
    }
    
    # Force metadata consistency
    if(nrow(chromPeaks(mse)) > 0) {
        cpd <- chromPeakData(mse); cpd$ms_level <- target_lvl; chromPeakData(mse) <- cpd
    }
    
    # 4. ALIGNMENT (Obiwarp)
    if(nrow(chromPeaks(mse)) > 0) {
        log_message("Performing Obiwarp alignment...")
        
        # A. Initial Grouping
        pdp <- PeakDensityParam(sampleGroups = pd$sample_group, minFraction=0.5, bw=30)
        mse <- groupChromPeaks(mse, param=pdp, msLevel=target_lvl)
        
        # B. Retention Time Correction (Obiwarp)
        # binSize=0.6 matches the documentation example
        mse <- adjustRtime(mse, param=ObiwarpParam(binSize=0.6), msLevel=target_lvl)
        
        # C. Re-group and Fill peaks (Standard workflow)
        mse <- groupChromPeaks(mse, param=pdp, msLevel=target_lvl)
        mse <- fillChromPeaks(mse, param=ChromPeakAreaParam(), msLevel=target_lvl)
    }

    # 5. EXTRACT "AFTER" DATA & DEVIATION
    log_message("Extracting Post-Alignment Chromatograms and Deviations...")
    
    # --- A. DEVIATION DATA (Bottom Plot) ---
    # Deviation = Adjusted - Raw
    rt_adj <- rtime(mse, adjusted = TRUE)
    rt_raw <- rtime(mse, adjusted = FALSE)
    file_indices <- fromFile(mse)
    
    output$results$retention_time_deviation <- lapply(1:length(file_paths), function(i) {
        mask <- file_indices == i
        # Calculate diff in seconds
        diffs <- rt_adj[mask] - rt_raw[mask]
        times <- rt_raw[mask] / 60 # X-axis in minutes
        
        # Create DF
        df <- data.frame(rt_min = times, deviation_seconds = diffs)
        
        # Downsample for UI performance (keep 1 every 20 points, the curve is smooth anyway)
        # This prevents sending megabytes of redundant line data to the browser
        if(nrow(df) > 1000) {
             df <- df[seq(1, nrow(df), 20), ]
        }
        
        list(filename = basename(file_paths[i]), data = df)
    })

    # --- B. AFTER ALIGNMENT CHROMATOGRAMS (Middle Plot) ---
    # We regenerate the chromatogram object, this time it will use adjusted times automatically
    bpc_adj_obj <- chromatogram(mse, aggregationFun = "max", msLevel = target_lvl)
    
    # Note: We must explicitly pull the adjusted times for the plotting data
    # The 'chromatogram' object might still hold raw times in some backend versions, 
    # so we manually stitch the BPC intensity with the adjusted time vector.
    
    # Get all adjusted times split by file
    adj_rt_list <- split(rtime(mse, adjusted=TRUE), fromFile(mse))
    
    output$results$chromatograms_after_alignment <- lapply(1:length(file_paths), function(i) {
        # Intensity from BPC
        ints <- intensity(bpc_adj_obj[1, i])
        # Correct Time from adjusted vector
        rts  <- adj_rt_list[[i]]
        
        min_len <- min(length(rts), length(ints))
        
        df <- data.frame(
            rt_min = rts[1:min_len] / 60,
            intensity = ints[1:min_len]
        ) %>% filter(!is.na(intensity))
        
        list(filename = basename(file_paths[i]), data = df)
    })
    
    # 6. STANDARD OUTPUTS (Mapping File 1 to standard fields for backward compatibility)
    # This ensures your existing Tables and Single-Plot logic still works
    output$results$raw_chromatogram_data <- output$results$chromatograms_after_alignment[[1]]$data
    output$results$smoothed_chromatogram_data <- output$results$chromatograms_after_alignment[[1]]$data

    peaks <- as.data.frame(chromPeaks(mse)) %>% filter(sample == 1)
    if(nrow(peaks) > 0) {
        peaks <- peaks %>% mutate(peak_number = row_number())
        total_area <- sum(peaks$into)
        output$results$quantitative_report <- peaks %>%
            mutate(rt_minutes = rt/60, peak_area = into, area_percent = (into/total_area)*100) %>%
            select(peak_number, rt_minutes, peak_area, area_percent)
        output$results$integrated_peaks_details <- peaks %>%
            mutate(rt_apex_min = rt/60, rt_start_min = rtmin/60, rt_end_min = rtmax/60, peak_height=maxo, peak_area=into)
            
        # Top 50 Spectra Logic (File 1)
        mse_sub <- mse[1]
        sps_sub <- filterMsLevel(spectra(mse_sub), msLevel=target_lvl)
        top_peaks <- peaks %>% arrange(desc(maxo)) %>% head(50)
        output$results$top_spectra_data <- lapply(1:nrow(top_peaks), function(i) {
             p <- top_peaks[i,]
             target_rt <- p$rt
             spec_win <- filterRt(sps_sub, rt=c(target_rt-2, target_rt+2))
             df <- data.frame(mz=numeric(), relative_intensity=numeric())
             if(length(spec_win) > 0) {
                 ii <- which.min(abs(rtime(spec_win) - target_rt))
                 d <- peaksData(spec_win)[[ii]]
                 if(nrow(d) > 0) {
                     max_i <- max(d[,2])
                     df <- data.frame(mz=d[,1], intensity=d[,2]) %>%
                         mutate(relative_intensity = (intensity/max_i)*100) %>%
                         filter(relative_intensity >= 1) %>% select(mz, relative_intensity)
                 }
             }
             list(peak_number=p$peak_number, spectrum_data=df)
        })
    } else {
        output$results$quantitative_report <- list()
        output$results$integrated_peaks_details <- list()
        output$results$top_spectra_data <- list()
    }
    
    return(output)
}

# --- ENTRY POINT ---
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("Usage: Rscript xcms_analysis.R <files> <output> [noise]")
  
  files <- trimws(strsplit(args[1], ",")[[1]])
  out_path <- args[2]
  noise <- if (length(args) >= 3) as.numeric(args[3]) else 2.5
  
  output_data <- list(status = "processing", results = list(), error = NULL)
  
  tryCatch({
      if (length(files) == 1) {
          res <- process_single_file(files[1], noise)
      } else {
          res <- process_multi_file(files, noise)
      }
      
      output_data$results <- res$results
      output_data$results$package_citations <- get_package_citations(c("Spectra", "xcms", "MsExperiment"))
      output_data$status <- "success"
      write(toJSON(output_data, auto_unbox = TRUE, pretty = TRUE), file = out_path)
  }, error = function(e) {
      write(toJSON(list(status="error", error=e$message), auto_unbox=TRUE), file=out_path)
      quit(status=1)
  })
}
# --- ENTRY POINT ---
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("Usage: Rscript xcms_analysis.R <files> <output> [noise]")
  
  files <- trimws(strsplit(args[1], ",")[[1]])
  out_path <- args[2]
  noise <- if (length(args) >= 3) as.numeric(args[3]) else 2.5
  
  output_data <- list(status = "processing", results = list(), error = NULL)
  
  tryCatch({
      if (length(files) == 1) {
          res <- process_single_file(files[1], noise)
      } else {
          res <- process_multi_file(files, noise)
      }
      
      output_data$results <- res$results
      output_data$results$package_citations <- get_package_citations(c("Spectra", "xcms", "MsExperiment"))
      output_data$status <- "success"
      write(toJSON(output_data, auto_unbox = TRUE, pretty = TRUE), file = out_path)
  }, error = function(e) {
      write(toJSON(list(status="error", error=e$message), auto_unbox=TRUE), file=out_path)
      quit(status=1)
  })
}