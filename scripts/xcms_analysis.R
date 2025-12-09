# Suppress package startup messages
suppressPackageStartupMessages({
  library(jsonlite)
  library(Spectra)
  library(mzR)
  library(dplyr)
  library(tidyr)
  library(MsCoreUtils)
})

source('./scripts/helpers.R')

# --- 1. MEMORY HELPERS ---

# Downsampler: Reduces massive vectors to ~2000 points for UI plotting
# This is crucial to keep JSON size and RAM usage low.
downsample_chrom <- function(rt, int, n_out = 2000) {
  if (length(rt) <= n_out) {
    return(data.frame(rt_min = rt / 60, intensity = int))
  }
  
  idx <- floor(seq(1, length(rt), length.out = n_out + 1))
  res_rt <- numeric(n_out)
  res_int <- numeric(n_out)
  
  for (i in 1:n_out) {
    s <- idx[i]; e <- idx[i+1]
    if (s < e) {
      res_rt[i] <- rt[s]
      val <- max(int[s:e], na.rm = TRUE)
      res_int[i] <- if(is.infinite(val)) 0 else val
    }
  }
  data.frame(rt_min = res_rt / 60, intensity = res_int)
}

# --- 2. PEAK DETECTION LOGIC (From your Version 1) ---

empty_peaks_df <- function() {
  data.frame(
    rt_apex = numeric(), rt_start = numeric(), rt_end = numeric(),
    peak_height = numeric(), peak_area = numeric(), stringsAsFactors = FALSE
  )
}

detect_and_integrate_peaks <- function(rt, intensity, noise_thresh,
                                       min_width_sec = 3, min_prominence = NULL) {
  stopifnot(length(rt) == length(intensity))
  
  # Basic derivative-based peak detection
  is_peak <- which(diff(sign(diff(intensity))) == -2) + 1
  apexes <- is_peak[intensity[is_peak] > noise_thresh]
  if (length(apexes) == 0) return(empty_peaks_df())
  
  get_valley <- function(idx, dir){
    j <- idx
    while (j > 1 && j < length(intensity) &&
           ((dir < 0 && intensity[j-1] <= intensity[j]) ||
            (dir > 0 && intensity[j+1] <= intensity[j]))) {
      j <- j + dir
    }
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
    
    data.frame(
      rt_apex = rt[apex_idx],
      rt_start = rt[left_idx],
      rt_end = rt[right_idx],
      peak_height = intensity[apex_idx],
      peak_area = peak_area,
      stringsAsFactors = FALSE
    )
  })
  
  # Clean up NULLs and bind
  peak_list <- peak_list[!sapply(peak_list, is.null)]
  valid_peaks <- do.call(rbind, peak_list)
  if (is.null(valid_peaks) || nrow(valid_peaks) == 0) return(empty_peaks_df())
  dplyr::distinct(valid_peaks, rt_start, rt_end, .keep_all = TRUE)
}

# --- 3. MAIN PROCESSING FUNCTION ---

perform_feature_analysis <- function(file_paths, noise_threshold_percent = 2.5) {
  
  output_data <- list(status = "processing", results = list(), error = NULL)
  
  chrom_list <- list()      # Stores downsampled chromatograms for all files
  deviation_list <- list()  # Stores alignment info
  
  # Reference data (File 1) used for alignment calculation
  ref_rt <- NULL
  ref_int <- NULL
  
  # --- LOOP: Process Each File Sequentially ---
  for (i in seq_along(file_paths)) {
    fpath <- file_paths[i]
    fname <- basename(fpath)
    
    tryCatch({
      # A. Load Data (Spectra object - OnDisk)
      # This is memory efficient as it reads metadata only initially
      exp_spectra <- Spectra(fpath, backend = MsBackendMzR())
      
      # B. Extract TIC
      # tic() sums all ions in a scan. Works for MS1 or MS2 automatically.
      raw_rt <- rtime(exp_spectra)
      raw_int <- tic(exp_spectra)
      
      # C. Smooth
      sg_coeffs <- MsCoreUtils::coefSG(hws = 5L, k = 3L)
      smooth_int <- MsCoreUtils::smooth(raw_int, cf = sg_coeffs)
      
      # D. Alignment Calculation (Linear Shift vs File 1)
      shift_sec <- 0
      
      if (i == 1) {
        # Store Reference
        ref_rt <- raw_rt
        ref_int <- smooth_int
      } else {
        # Align File 'i' to File 1 using Cross-Correlation
        # This is much lighter than XCMS Obiwarp but effective for visualization
        
        # Interpolate to common grid for CCF (needed if scan rates differ)
        common_time <- seq(min(ref_rt), max(ref_rt), by = 0.5) # 0.5s bins
        
        # Linear approximation
        y1 <- approx(ref_rt, ref_int, xout = common_time, rule = 2)$y
        y2 <- approx(raw_rt, smooth_int, xout = common_time, rule = 2)$y
        
        # Cross-correlation
        cc <- ccf(y1, y2, plot = FALSE, lag.max = 200) # Max lag +/- 100 steps (50s)
        best_lag_idx <- which.max(cc$acf)
        lag_val <- cc$lag[best_lag_idx]
        
        # Calculate shift in seconds
        shift_sec <- lag_val * 0.5
      }
      
      # Store deviation info for UI
      # We create a simple line data frame to show "Constant Shift"
      dev_df <- data.frame(
        rt_min = c(min(raw_rt), max(raw_rt)) / 60, 
        deviation_seconds = c(shift_sec, shift_sec)
      )
      deviation_list[[i]] <- list(filename = fname, data = dev_df)
      
      # Apply shift to chromatogram for overlay plotting
      aligned_rt <- raw_rt + shift_sec
      
      # E. Downsample Chromatogram (Save Memory)
      df_small <- downsample_chrom(aligned_rt, smooth_int)
      chrom_list[[i]] <- list(filename = fname, data = df_small)
      
      # --- SPECIAL PROCESSING FOR FILE 1 (THE SAMPLE) ---
      # We only do detailed Peak Detection & Spectra Extraction on File 1
      # to keep the JSON light and the report focused.
      if (i == 1) {
        max_val <- max(smooth_int, na.rm = TRUE)
        abs_noise <- max_val * (noise_threshold_percent / 100)
        
        # 1. Detect Peaks
        peaks_df <- detect_and_integrate_peaks(raw_rt, smooth_int, abs_noise)
        
        if (nrow(peaks_df) > 0) {
          # Number peaks
          peaks_df <- peaks_df %>% arrange(rt_apex) %>% mutate(peak_number = row_number())
          
          # Quantitative Report
          total_area <- sum(peaks_df$peak_area)
          output_data$results$quantitative_report <- peaks_df %>%
            mutate(
              area_percent = if (total_area > 0) (peak_area / total_area) * 100 else 0,
              rt_minutes = round(rt_apex / 60, 3),
              peak_area = round(peak_area, 0),
              area_percent = round(area_percent, 2)
            ) %>%
            select(peak_number, rt_minutes, peak_area, area_percent)
            
          # Integrated Peak Details (for Plot shading)
          output_data$results$integrated_peaks_details <- peaks_df %>%
            mutate(rt_apex_min = rt_apex/60, rt_start_min = rt_start/60, rt_end_min = rt_end/60)
          
          # 2. Extract Top 50 Spectra
          top_50 <- peaks_df %>% arrange(desc(peak_height)) %>% head(50)
          spectra_res <- list()
          
          for(k in 1:nrow(top_50)) {
             p <- top_50[k,]
             
             # Find scan index closest to apex
             idx <- which.min(abs(raw_rt - p$rt_apex))
             
             # Extract scan data
             scan_obj <- exp_spectra[idx]
             pdata <- peaksData(scan_obj)[[1]]
             
             spec_df <- data.frame(mz=numeric(), relative_intensity=numeric())
             
             if (!is.null(pdata) && nrow(pdata) > 0) {
                 mz_vals <- pdata[,1]
                 int_vals <- pdata[,2]
                 max_i <- max(int_vals)
                 
                 # Filter noise (<1%) to keep JSON small
                 keep_mask <- int_vals > (max_i * 0.01)
                 
                 spec_df <- data.frame(
                    mz = mz_vals[keep_mask],
                    relative_intensity = (int_vals[keep_mask] / max_i) * 100
                 )
             }
             
             spectra_res[[k]] <- list(
                peak_number = p$peak_number,
                spectrum_data = spec_df
             )
          }
          output_data$results$top_spectra_data <- spectra_res
          
        } else {
          # No peaks found
          output_data$results$quantitative_report <- list()
          output_data$results$top_spectra_data <- list()
          output_data$results$integrated_peaks_details <- list()
        }
        
        # Save standard single-file plotting data for backward compatibility
        output_data$results$raw_chromatogram_data <- chrom_list[[1]]$data
        output_data$results$smoothed_chromatogram_data <- chrom_list[[1]]$data
      }
      
      # Clean up memory immediately after processing file
      rm(exp_spectra, raw_rt, raw_int, smooth_int)
      gc()
      
    }, error = function(e) {
      log_message(paste("Error processing file", fname, ":", e$message))
    })
  }
  
  # --- 4. FINALIZE ---
  
  # Add the multi-file lists to the result
  output_data$results$chromatograms_after_alignment <- chrom_list
  output_data$results$retention_time_deviation <- deviation_list
  
  if (is.null(output_data$error)) { output_data$status <- "success" }
  return(output_data)
}

# --- SCRIPT EXECUTION ---
if (!interactive()) {
  tryCatch({
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
      stop("Usage: Rscript xcms_analysis.R <files> <output> [noise]")
    }
    
    # Split comma-separated files
    file_list <- trimws(strsplit(args[1], ",")[[1]])
    output_path <- args[2]
    noise_thresh <- if (length(args) >= 3) as.numeric(args[3]) else 2.5
    
    result <- perform_feature_analysis(file_paths = file_list, noise_threshold_percent = noise_thresh)
    
    # Use write to stream to file (safer than generating massive string)
    json_output <- toJSON(result, auto_unbox = TRUE, pretty = TRUE)
    write(json_output, file = output_path)
    
  }, error = function(e) {
    error_json <- toJSON(list(status = "error", error = e$message), auto_unbox = TRUE)
    try(write(error_json, file = commandArgs(trailingOnly=TRUE)[2]), silent = TRUE)
    quit(status = 1)
  })
}