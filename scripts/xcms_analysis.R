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
  library(BiocParallel) # Added for memory control
})

# --- MEMORY OPTIMIZATION: ENFORCE SERIAL EXECUTION ---
# This prevents R from forking processes, which doubles memory usage per core.
register(SerialParam())

# --- HELPERS ---
source("./scripts/helpers.R")

# --- MEMORY OPTIMIZATION: DOWNSAMPLER ---
# Reduces 50,000+ points to ~2,000 for UI plotting.
# This massively reduces the JSON size and memory required to generate it.
downsample_chrom <- function(rt, int, n_out = 2000) {
  if (length(rt) <= n_out) {
    return(data.frame(rt_min = rt, intensity = int))
  }
  
  # Create bins
  breaks <- seq(1, length(rt), length.out = n_out + 1)
  idx <- floor(breaks)
  
  # fast binning
  rt_bin <- rt[idx[-length(idx)]] # Take start time of bin
  
  # Calculate max intensity per bin loop (vectorized aggregation is memory heavy)
  # We use a simple indexing approach to avoid creating large matrices
  int_bin <- vapply(1:n_out, function(i) {
    s <- idx[i]
    e <- idx[i+1]
    # Protect against index out of bounds
    if (s >= e) return(int[s]) 
    max(int[s:e], na.rm = TRUE)
  }, numeric(1))
  
  # Scale RT to minutes here to save compute later
  data.frame(rt_min = rt_bin / 60, intensity = int_bin)
}

# ==============================================================================
# LOGIC A: SINGLE FILE
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
    
    # Filter NULLs before rbind to save memory
    peak_list <- peak_list[!sapply(peak_list, is.null)]
    valid_peaks <- do.call(rbind, peak_list)
    
    if (is.null(valid_peaks) || nrow(valid_peaks) == 0) return(empty_peaks_df())
    distinct(valid_peaks, rt_start, rt_end, .keep_all = TRUE)
}

process_single_file <- function(mzMl_file_path, noise_threshold_percent) {
    log_message("Single file detected.")
    output <- list(results = list())
    
    # Use on-disk backend
    exp_spectra <- Spectra(mzMl_file_path, backend = MsBackendMzR())
    
    # Extract TIC (low memory footprint compared to full spectra)
    raw_rt <- rtime(exp_spectra)
    raw_int <- tic(exp_spectra)
    
    # Smooth
    sg_coeffs <- MsCoreUtils::coefSG(hws = 5L, k = 3L)
    smooth_int <- MsCoreUtils::smooth(raw_int, cf = sg_coeffs)
    
    # Peak Detection
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
        
        # Optimize spectrum extraction
        output$results$top_spectra_data <- lapply(1:nrow(top_peaks), function(i) {
             p <- top_peaks[i,]
             # Find index closest to apex
             idx <- which.min(abs(raw_rt - p$rt_apex))
             
             # Extract ONLY the specific scan needed
             spec_subset <- exp_spectra[idx]
             vals <- peaksData(spec_subset)[[1]]
             
             df <- data.frame(mz=numeric(), relative_intensity=numeric())
             if(!is.null(vals) && nrow(vals) > 0) {
                 max_i <- max(vals[,2], na.rm=TRUE)
                 if(max_i > 0) {
                     # Filter strictly to save JSON size
                     df <- data.frame(mz = vals[,1], intensity = vals[,2]) %>% 
                           mutate(relative_intensity = (intensity/max_i)*100) %>% 
                           filter(relative_intensity >= 1.0) %>% 
                           select(mz, relative_intensity)
                 }
             }
             list(peak_number = p$peak_number, spectrum_data = df)
        })
    } else {
        output$results$quantitative_report <- list()
        output$results$integrated_peaks_details <- list()
        output$results$top_spectra_data <- list()
    }
    
    # DOWNSAMPLE BEFORE EXPORT
    output$results$raw_chromatogram_data <- downsample_chrom(raw_rt, raw_int)
    output$results$smoothed_chromatogram_data <- downsample_chrom(raw_rt, smooth_int)
    
    gc() # Force cleanup
    return(output)
}

# ==============================================================================
# LOGIC B: MULTI-FILE
# ==============================================================================

process_multi_file <- function(file_paths, noise_threshold_percent) {
    log_message("Multiple files detected. Using XCMS...")
    output <- list(results = list())
    
    pd <- data.frame(sample_name = paste0("File_", 1:length(file_paths)), sample_group = "Study")
    
    # MsExperiment with On-Disk backend implicitly
    mse <- readMsExperiment(spectraFiles = file_paths, sampleData = pd)
    
    lvls <- unique(msLevel(spectra(mse)))
    target_lvl <- if(1L %in% lvls) 1L else 2L
    
    # 2. EXTRACT "BEFORE" DATA
    log_message("Extracting Pre-Alignment Chromatograms...")
    
    # BPC extraction
    bpc_raw_obj <- chromatogram(mse, aggregationFun = "max", msLevel = target_lvl)
    
    # Helper with Downsampling
    extract_chrom_data <- function(chrom_obj) {
        lapply(1:length(file_paths), function(i) {
            chr <- chrom_obj[1, i]
            # Use downsampler here
            df <- downsample_chrom(rtime(chr), intensity(chr))
            # Clean NAs
            df <- df %>% filter(!is.na(intensity))
            list(filename = basename(file_paths[i]), data = df)
        })
    }
    
    output$results$chromatograms_before_alignment <- extract_chrom_data(bpc_raw_obj)
    
    # Cleanup early
    rm(bpc_raw_obj)
    gc()

    # 3. PEAK DETECTION
    # Re-calculate max val roughly from existing extracted data to save memory
    all_ints <- unlist(lapply(output$results$chromatograms_before_alignment, function(x) x$data$intensity))
    max_val <- if(length(all_ints) > 0) max(all_ints, na.rm=TRUE) else 1000
    rm(all_ints) # Cleanup
    
    abs_noise <- if(is.finite(max_val)) (max_val/5)*(noise_threshold_percent/100) else 1000
    
    cwp <- CentWaveParam(ppm=50, peakwidth=c(5, 60), snthresh=5, noise=abs_noise, prefilter=c(3, abs_noise))
    mse <- findChromPeaks(mse, param=cwp, msLevel=target_lvl)
    
    if(nrow(chromPeaks(mse)) == 0) {
        cwp_loose <- CentWaveParam(ppm=50, peakwidth=c(5,60), snthresh=3, noise=0, prefilter=c(3, 100))
        mse <- findChromPeaks(mse, param=cwp_loose, msLevel=target_lvl)
    }
    
    if(nrow(chromPeaks(mse)) > 0) {
        cpd <- chromPeakData(mse); cpd$ms_level <- target_lvl; chromPeakData(mse) <- cpd
    }
    
    gc() # Critical GC after peak detection

    # 4. ALIGNMENT
    if(nrow(chromPeaks(mse)) > 0) {
        log_message("Performing Obiwarp alignment...")
        pdp <- PeakDensityParam(sampleGroups = pd$sample_group, minFraction=0.5, bw=30)
        mse <- groupChromPeaks(mse, param=pdp, msLevel=target_lvl)
        mse <- adjustRtime(mse, param=ObiwarpParam(binSize=0.6), msLevel=target_lvl)
        mse <- groupChromPeaks(mse, param=pdp, msLevel=target_lvl)
        mse <- fillChromPeaks(mse, param=ChromPeakAreaParam(), msLevel=target_lvl)
    }
    
    gc() # Critical GC after alignment

    # 5. EXTRACT "AFTER" DATA
    log_message("Extracting Post-Alignment Chromatograms...")
    
    # Deviation Data
    rt_adj <- rtime(mse, adjusted = TRUE)
    rt_raw <- rtime(mse, adjusted = FALSE)
    file_indices <- fromFile(mse)
    
    output$results$retention_time_deviation <- lapply(1:length(file_paths), function(i) {
        mask <- file_indices == i
        diffs <- rt_adj[mask] - rt_raw[mask]
        times <- rt_raw[mask]
        
        # DOWNSAMPLE DEVIATION DATA
        # We manually bin this because it's not a chromatogram
        n_pts <- length(times)
        if(n_pts > 2000) {
            idx <- seq(1, n_pts, length.out = 2000)
            times <- times[idx]
            diffs <- diffs[idx]
        }
        
        df <- data.frame(rt_min = times / 60, deviation_seconds = diffs)
        list(filename = basename(file_paths[i]), data = df)
    })
    
    # Cleanup large vectors
    rm(rt_adj, rt_raw, file_indices)
    gc()

    # Chromatograms After Alignment
    bpc_adj_obj <- chromatogram(mse, aggregationFun = "max", msLevel = target_lvl)
    
    # Adjusted times are heavy to pull all at once, so we do it per file inside the loop
    output$results$chromatograms_after_alignment <- lapply(1:length(file_paths), function(i) {
        chr <- bpc_adj_obj[1, i]
        # Use downsampler
        df <- downsample_chrom(rtime(chr), intensity(chr))
        df <- df %>% filter(!is.na(intensity))
        list(filename = basename(file_paths[i]), data = df)
    })
    
    rm(bpc_adj_obj)
    gc()

    # 6. STANDARD OUTPUTS (Mapped to File 1)
    # Reuse the downsampled data we already created
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
            
        # Top 50 Spectra Logic
        mse_sub <- mse[1]
        sps_sub <- filterMsLevel(spectra(mse_sub), msLevel=target_lvl)
        top_peaks <- peaks %>% arrange(desc(maxo)) %>% head(50)
        
        output$results$top_spectra_data <- lapply(1:nrow(top_peaks), function(i) {
             p <- top_peaks[i,]
             target_rt <- p$rt
             # Extract small window
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
      
      # Clean up R environment before massive JSON string allocation
      rm(res)
      gc()
      
      write(toJSON(output_data, auto_unbox = TRUE, pretty = TRUE), file = out_path)
  }, error = function(e) {
      write(toJSON(list(status="error", error=e$message), auto_unbox=TRUE), file=out_path)
      quit(status=1)
  })
}