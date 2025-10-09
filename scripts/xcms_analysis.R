# Suppress package startup messages for cleaner logs
suppressPackageStartupMessages({
  library(jsonlite)
  library(Spectra)
  library(mzR)
  library(dplyr)
  library(tidyr)
  library(MsCoreUtils)
})

source('./scripts/helpers.R')

perform_feature_analysis <- function(mzMl_file_path, noise_threshold_percent = 2.5) {
  
  output_data <- list(
    status = "processing",
    results = list(),
    error = NULL
  )
  
  # Helper: empty peak table with required columns
  empty_peaks_df <- function() {
    data.frame(
      rt_apex = numeric(),
      rt_start = numeric(),
      rt_end = numeric(),
      peak_height = numeric(),
      peak_area = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  # Advanced peak detection and integration function
  detect_and_integrate_peaks <- function(rt, intensity, noise_thresh,
                                         min_width_sec = 3, min_prominence = NULL) {
    stopifnot(length(rt) == length(intensity))
    
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
      
      width_sec <- rt[right_idx] - rt[left_idx]
      if (!is.null(min_width_sec) && width_sec < min_width_sec) return(NULL)
      
      left_min  <- min(intensity[left_idx:apex_idx], na.rm = TRUE)
      right_min <- min(intensity[apex_idx:right_idx], na.rm = TRUE)
      prom <- intensity[apex_idx] - max(left_min, right_min)
      if (!is.null(min_prominence) && prom < min_prominence) return(NULL)
      
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
    
    valid_peaks <- do.call(rbind, peak_list)
    if (is.null(valid_peaks) || nrow(valid_peaks) == 0) return(empty_peaks_df())
    dplyr::distinct(valid_peaks, rt_start, rt_end, .keep_all = TRUE)
  }
  
  # --- 1. DATA LOADING, SMOOTHING, AND PEAK PROCESSING ---
  tryCatch({
    log_message("Loading MS data and chromatogram...")
    exp_spectra <- Spectra(mzMl_file_path, backend = MsBackendMzR())
    raw_chrom_data <- data.frame(rt_sec = rtime(exp_spectra), intensity = tic(exp_spectra))
    
    log_message("Smoothing chromatogram data...")
    sg_coeffs <- MsCoreUtils::coefSG(hws = 5L, k = 3L)
    smoothed_intensity <- MsCoreUtils::smooth(raw_chrom_data$intensity, cf = sg_coeffs)
    smoothed_chrom_data <- raw_chrom_data
    smoothed_chrom_data$intensity <- smoothed_intensity
    
    max_intensity_smoothed <- max(smoothed_intensity, na.rm = TRUE)
    absolute_noise_threshold <- max_intensity_smoothed * (noise_threshold_percent / 100)
    
    log_message(paste("Starting peak detection with", noise_threshold_percent, "% noise threshold..."))
    all_peaks_numbered <- detect_and_integrate_peaks(
      smoothed_chrom_data$rt_sec, smoothed_chrom_data$intensity, absolute_noise_threshold
    )
    if (nrow(all_peaks_numbered) > 0) {
        all_peaks_numbered <- all_peaks_numbered %>% arrange(rt_apex) %>% mutate(peak_number = dplyr::row_number())
    }
    
    log_message(paste("Successfully integrated and numbered", nrow(all_peaks_numbered), "peaks."))
    output_data$results$integrated_peaks <- all_peaks_numbered
    output_data$results$smoothed_chrom_data_internal <- smoothed_chrom_data
    
  }, error = function(e) {
    output_data$status <<- "error"; output_data$error <<- e$message; return(output_data)
  })
  if (!is.null(output_data$error)) return(output_data)
  
  # --- 2. QUANTITATIVE REPORT & SPECTRA DATA EXTRACTION ---
  tryCatch({
    log_message("Generating quantitative report...")
    all_peaks_numbered <- output_data$results$integrated_peaks
    
    if (nrow(all_peaks_numbered) == 0) {
      log_message("No peaks were found. Skipping report and spectra extraction.")
      output_data$results$quantitative_report <- data.frame()
      output_data$results$top_spectra_data <- list()
    } else {
      total_area <- sum(all_peaks_numbered$peak_area, na.rm = TRUE)
      output_data$results$quantitative_report <- all_peaks_numbered %>%
        mutate(
          area_percent = if (total_area > 0) (peak_area / total_area) * 100 else 0,
          rt_minutes = round(rt_apex / 60, 3),
          peak_area = round(peak_area, 0),
          area_percent = round(area_percent, 2)
        ) %>%
        select(peak_number, rt_minutes, peak_area, area_percent)
      
      log_message("Identifying and extracting spectra for the top 50 most intense peaks...")
      top_50_peaks_data <- all_peaks_numbered %>% arrange(desc(peak_height)) %>% head(50)
      
      if (nrow(top_50_peaks_data) > 0) {
        spec_rt <- rtime(exp_spectra)
        top_features <- top_50_peaks_data %>%
          rowwise() %>%
          mutate(spectrum_index = which.min(abs(spec_rt - rt_apex))) %>%
          ungroup() %>%
          select(peak_number, spectrum_index)
        output_data$results$top_features <- top_features

        top_spectra_objects <- exp_spectra[top_features$spectrum_index]
        peaks_list <- peaksData(top_spectra_objects)
        
        spectra_results_list <- lapply(1:nrow(top_features), function(i) {
          spec_df <- as.data.frame(peaks_list[[i]])
          if (ncol(spec_df) >= 2) {
              colnames(spec_df)[1:2] <- c("mz", "intensity")
              if(nrow(spec_df) > 0 && max(spec_df$intensity, na.rm=TRUE) > 0) {
                spec_df <- spec_df %>%
                    mutate(relative_intensity = (intensity / max(intensity, na.rm = TRUE)) * 100) %>%
                    filter(relative_intensity >= 1.0) %>%
                    select(mz, relative_intensity)
              } else {
                spec_df <- data.frame(mz=numeric(), relative_intensity=numeric())
              }
          } else {
             spec_df <- data.frame(mz=numeric(), relative_intensity=numeric())
          }
          list(
            peak_number = top_features$peak_number[i],
            spectrum_data = spec_df
          )
        })
        
        # <<< --- THIS IS THE FIX --- >>>
        # Use the name 'top_spectra_data' that the UI is expecting.
        output_data$results$top_spectra_data <- spectra_results_list

      } else {
        output_data$results$top_spectra_data <- list()
      }
    }
  }, error = function(e) {
    output_data$status <<- "error"; output_data$error <<- e$message; return(output_data)
  })
  if (!is.null(output_data$error)) return(output_data)
  
  # --- 3. FINALIZE AND RETURN ---
  output_data$results$raw_chromatogram_data <- raw_chrom_data %>% mutate(rt_min = rt_sec/60) %>% select(rt_min, intensity)
  output_data$results$smoothed_chromatogram_data <- output_data$results$smoothed_chrom_data_internal %>% mutate(rt_min = rt_sec/60) %>% select(rt_min, intensity)
  output_data$results$integrated_peaks_details <- output_data$results$integrated_peaks %>%
      mutate(rt_apex_min = rt_apex/60, rt_start_min = rt_start/60, rt_end_min = rt_end/60, peak_height = peak_height)
  
  # Clean internals
  output_data$results$integrated_peaks <- NULL
  output_data$results$smoothed_chrom_data_internal <- NULL
  
  if (is.null(output_data$error)) { output_data$status <- "success" }
  return(output_data)
}

# --- SCRIPT EXECUTION ---
if (!interactive()) {
  tryCatch({
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2 || length(args) > 3) {
      stop("Usage: Rscript xcms_analysis.R <input.mzML> <output.json> [noise_threshold_percent]")
    }
    input_path <- args[1]; output_path <- args[2]
    noise_thresh <- if (length(args) == 3) as.numeric(args[3]) else 3
    
    result <- perform_feature_analysis(mzMl_file_path = input_path, noise_threshold_percent = noise_thresh)
    json_output <- toJSON(result, auto_unbox = TRUE, pretty = TRUE)
    write(json_output, file = output_path)
    
  }, error = function(e) {
    error_json <- toJSON(list(status = "error", error = e$message), auto_unbox = TRUE)
    try(write(error_json, file = commandArgs(trailingOnly=TRUE)[2]), silent = TRUE)
    quit(status = 1)
  })
}