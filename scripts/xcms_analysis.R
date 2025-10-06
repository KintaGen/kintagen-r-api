# Suppress package startup messages for cleaner logs
suppressPackageStartupMessages({
  library(jsonlite)
  library(Spectra)
  library(mzR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("./scripts/helpers.R")

perform_gcms_analysis <- function(mzMl_file_path) {
  
  output_data <- list(
    status = "processing",
    results = list(),
    error = NULL
  )

  # Helper function is unchanged
  detect_and_integrate_peaks <- function(rt, intensity, noise_threshold_percent = 1.0) {
    # ... (function code is identical)
    if (length(rt) != length(intensity)) stop("Vectors must have the same length.")
    max_intensity <- max(intensity, na.rm = TRUE)
    if (max_intensity == 0) return(data.frame())
    noise_threshold <- max_intensity * (noise_threshold_percent / 100)
    is_peak <- which(diff(sign(diff(intensity))) == -2) + 1
    apexes <- is_peak[intensity[is_peak] > noise_threshold]
    if (length(apexes) == 0) return(data.frame())
    peak_list <- lapply(apexes, function(apex_idx) {
      left_idx <- apex_idx
      while (left_idx > 1 && intensity[left_idx - 1] < intensity[left_idx] && intensity[left_idx - 1] >= noise_threshold) {
        left_idx <- left_idx - 1
      }
      right_idx <- apex_idx
      while (right_idx < length(intensity) && intensity[right_idx + 1] < intensity[right_idx] && intensity[right_idx + 1] >= noise_threshold) {
        right_idx <- right_idx + 1
      }
      if (left_idx >= right_idx) return(NULL)
      peak_indices <- left_idx:right_idx
      peak_rt <- rt[peak_indices]; peak_int <- intensity[peak_indices]
      if(length(peak_rt) < 2) return(NULL)
      peak_area <- sum(diff(peak_rt) * (peak_int[-1] + peak_int[-length(peak_int)]) / 2)
      return(data.frame(rt_apex=rt[apex_idx], rt_start=rt[left_idx], rt_end=rt[right_idx], peak_height=intensity[apex_idx], peak_area=peak_area, integration_indices=I(list(peak_indices))))
    })
    valid_peaks <- do.call(rbind, peak_list)
    if (is.null(valid_peaks) || nrow(valid_peaks) == 0) return(data.frame())
    return(valid_peaks %>% distinct(rt_start, rt_end, .keep_all = TRUE))
  }

  # --- 1. DATA LOADING AND PEAK PROCESSING ---
  tryCatch({
    log_message("Loading MS data and chromatogram...")
    exp_spectra <- Spectra(mzMl_file_path, backend = MsBackendMzR())
    
    chrom_data <- data.frame(rt_sec = rtime(exp_spectra), intensity = tic(exp_spectra))

    log_message("Starting automated peak detection and integration...")
    unprocessed_peaks_df <- detect_and_integrate_peaks(chrom_data$rt_sec, chrom_data$intensity, noise_threshold_percent = 1.0)
    
    if (nrow(unprocessed_peaks_df) == 0) stop("No peaks detected.")
    
    # <<< --- KEY CHANGE 1: NUMBER ALL PEAKS BY RETENTION TIME FIRST --- >>>
    # This creates the definitive peak list that everything else will be based on.
    all_peaks_numbered <- unprocessed_peaks_df %>%
      arrange(rt_apex) %>%
      mutate(peak_number = row_number())
      
    log_message(paste("Successfully integrated and numbered", nrow(all_peaks_numbered), "peaks."))
    output_data$results$integrated_peaks <- all_peaks_numbered # Store this complete table

  }, error = function(e) {
    output_data$status <<- "error"; output_data$error <<- e$message; return(output_data)
  })
  if (!is.null(output_data$error)) return(output_data)

  # --- 2. QUANTITATIVE REPORT & TOP FEATURE SELECTION ---
  tryCatch({
    log_message("Generating quantitative report...")
    all_peaks_numbered <- output_data$results$integrated_peaks
    total_area <- sum(all_peaks_numbered$peak_area, na.rm = TRUE)
    if (total_area == 0) stop("Total integrated area is zero.")
    
    # The quant report is now directly created from the pre-numbered table
    output_data$results$quantitative_report <- all_peaks_numbered %>%
      mutate(
        area_percent = (peak_area / total_area) * 100,
        rt_minutes = round(rt_apex / 60, 3),
        peak_area = round(peak_area, 0),
        area_percent = round(area_percent, 2)
      ) %>%
      select(peak_number, rt_minutes, peak_area, area_percent)
    
    # <<< --- KEY CHANGE 2: SELECT TOP FEATURES BUT KEEP THEIR REAL PEAK NUMBER --- >>>
    log_message("Identifying spectra for the top 5 most intense peaks...")
    top_5_peaks_data <- all_peaks_numbered %>%
      arrange(desc(peak_height)) %>%
      head(5)

    # Now, build the `top_features` object for the UI from this selection
    all_runtimes_sec <- rtime(exp_spectra)
    all_tics <- tic(exp_spectra)
    max_tic_overall <- max(all_tics, na.rm = TRUE)
    
    top_features <- top_5_peaks_data %>%
      rowwise() %>%
      mutate(
        spectrum_index = which.min(abs(all_runtimes_sec - rt_apex)),
        retention_time_sec = all_runtimes_sec[spectrum_index],
        absolute_tic = all_tics[spectrum_index],
        relative_tic_percent = if(max_tic_overall > 0) (absolute_tic / max_tic_overall) * 100 else 0
      ) %>%
      ungroup() %>%
      select(peak_number, spectrum_index, retention_time_sec, relative_tic_percent) # Keep peak_number!
    
    output_data$results$top_features <- top_features

  }, error = function(e) {
    output_data$status <<- "error"; output_data$error <<- e$message; return(output_data)
  })
  if (!is.null(output_data$error)) return(output_data)

  # --- 3. GENERATE TOP 5 SPECTRA PLOT ---
  tryCatch({
    log_message("Generating annotated static plot of top 5 spectra...")
    
    top_features_to_plot <- output_data$results$top_features
    if (nrow(top_features_to_plot) == 0) stop("No features found to plot.")
    
    top_indices <- top_features_to_plot$spectrum_index
    top_spectra <- exp_spectra[top_indices]
    
    peaks_list <- peaksData(top_spectra)
    results_list <- lapply(1:nrow(top_features_to_plot), function(i) {
      df <- as.data.frame(peaks_list[[i]]); colnames(df) <- c("mz", "intensity")
      
      if (nrow(df) > 0 && max(df$intensity, na.rm = TRUE) > 0) {
        df <- df %>% mutate(relative_intensity = (intensity / max(df$intensity, na.rm = TRUE)) * 100) %>% filter(relative_intensity >= 1.0)
      } else { df$relative_intensity <- numeric(0) }
      
      # <<< --- KEY CHANGE 3: USE THE CORRECT PEAK NUMBER FOR THE LABEL --- >>>
      # Instead of using the loop index `i`, we use the actual `peak_number` from our data.
      correct_peak_number <- top_features_to_plot$peak_number[i]
      retention_time_min <- round(top_features_to_plot$retention_time_sec[i] / 60, 2)
      df$label <- paste0("#", correct_peak_number, " @ ", retention_time_min, " min")
      
      # ... rest of peak labeling logic is unchanged
      peaks_to_label_df <- data.frame()
      if (nrow(df) > 0) {
        molecular_ion_peak <- df %>% arrange(desc(mz)) %>% head(1) %>% mutate(peak_type = "Molecular Ion")
        fragment_peaks <- df %>% filter(mz != molecular_ion_peak$mz) %>% arrange(desc(relative_intensity)) %>% head(3) %>% mutate(peak_type = "Fragment")
        peaks_to_label_df <- rbind(molecular_ion_peak, fragment_peaks) %>% mutate(mz_label = format(round(mz, 2), nsmall = 2))
      }
      return(list(main_data = df, label_data = peaks_to_label_df))
    })
    
    combined_plot_df <- do.call(rbind, lapply(results_list, `[[`, "main_data"))
    combined_label_df <- do.call(rbind, lapply(results_list, `[[`, "label_data"))
    
    p_static <- ggplot(combined_plot_df, aes(x=mz, xend=mz, y=0, yend=relative_intensity)) +
      # ... ggplot layers are unchanged ...
      geom_segment(linewidth=0.7, color="gray20") +
      geom_text(data=combined_label_df, aes(x=mz, y=relative_intensity, label=mz_label, color=peak_type), angle=45, hjust=-0.1, vjust=-0.2, size=2.8, fontface="bold", show.legend=FALSE) +
      scale_color_manual(values = c("Fragment"="goldenrod3", "Molecular Ion"="firebrick")) +
      facet_wrap(~label, ncol=1) + 
      labs(title="Spectra of Top 5 Most Intense Peaks", x="m/z", y="Relative Intensity (%)") +
      theme_bw() + coord_cartesian(ylim=c(0, 150), expand=TRUE) + 
      theme(plot.title=element_text(hjust=0.5, face="bold"))
      
    output_data$results$top_5_spectra_plot_b64 <- gg_to_base64(p_static, width = 8, height = 12)
    log_message("Static spectra plot generation complete.")

  }, error = function(e) {
    output_data$status <<- "error"; output_data$error <<- e$message
  })
  if (!is.null(output_data$error)) return(output_data)

  # --- 4. FINALIZE AND RETURN ---
  output_data$results$chromatogram_data <- chrom_data %>% mutate(rt_min = rt_sec/60) %>% select(rt_min, intensity)
  output_data$results$integrated_peaks_details <- output_data$results$integrated_peaks %>% mutate(rt_apex_min = rt_apex/60, rt_start_min = rt_start/60, rt_end_min = rt_end/60, peak_height = peak_height)
  output_data$results$integrated_peaks <- NULL
  if (is.null(output_data$error)) { output_data$status <- "success" }
  return(output_data)
}

# --- SCRIPT EXECUTION ---
if (!interactive()) {
  tryCatch({
      args <- commandArgs(trailingOnly = TRUE)
      if (length(args) != 2) stop("Usage: Rscript <script.R> <input.mzML> <output.json>")
      input_path <- args[1]; output_path <- args[2]
      result <- perform_gcms_analysis(mzMl_file_path = input_path)
      json_output <- toJSON(result, auto_unbox = TRUE, pretty = TRUE)
      write(json_output, file = output_path)
  }, error = function(e) {
      error_json <- toJSON(list(status = "error", error = e$message), auto_unbox = TRUE)
      try(write(error_json, file = args[2]), silent = TRUE)
      quit(status = 1)
  })
}