perform_ms_analysis <- function(mzMl_file_path) {
  
  # A simple helper function to find local maxima (peaks)
  find_peaks <- function(x) {
    # A point is a peak if it's higher than its immediate neighbors
    peaks <- which(diff(sign(diff(x))) == -2) + 1
    return(peaks)
  }

  # --- 2. FEATURE DISCOVERY (with BPC Normalization) ---
  tryCatch({
    log_message("Starting feature discovery...")
    exp_spectra <- Spectra(mzMl_file_path, backend = MsBackendMzR())
    
    all_tics <- tic(exp_spectra)
    max_tic <- max(all_tics, na.rm = TRUE)
    relative_tic_percent <- if (max_tic > 0) (all_tics / max_tic) * 100 else all_tics

    spectra_summary <- data.frame(
      spectrum_index = 1:length(exp_spectra),
      precursor_mz = precursorMz(exp_spectra),
      retention_time_sec = rtime(exp_spectra),
      relative_tic_percent = relative_tic_percent
    )
    
    top_10_features <- spectra_summary %>%
      filter(relative_tic_percent > 0) %>%
      arrange(desc(relative_tic_percent)) %>%
      head(10)
      
    output_data$results$top_features <- top_10_features
    
    bpc_threshold_percent <- 2.0
    filtered_bpc_data <- spectra_summary %>%
      filter(relative_tic_percent >= bpc_threshold_percent)
    
    final_bpc_data <- filtered_bpc_data %>%
        select(rt_sec = retention_time_sec, relative_tic = relative_tic_percent)
        
    output_data$results$tic_data <- final_bpc_data
        
    log_message(paste("Found", nrow(top_10_features), "top features and prepared BPC data."))
    
  }, error = function(e) {
    log_message(paste("Error during feature discovery:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
    return(output_data)
  })
  
  if (!is.null(output_data$error)) return(output_data)

  # <<< --- NEW SECTION: 3. GENERATE BPC PLOT (BASE64 ENCODED) --- >>>
  tryCatch({
    log_message("Generating Base Peak Chromatogram plot...")
    
    if (nrow(filtered_bpc_data) < 2) {
      stop("Not enough data points to generate a BPC plot after filtering.")
    }

    p_bpc <- ggplot(filtered_bpc_data, aes(x = retention_time_sec / 60, y = relative_tic_percent)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      labs(
        title = "Base Peak Chromatogram (BPC)",
        subtitle = paste0("Showing data points > ", bpc_threshold_percent, "% relative intensity"),
        x = "Retention Time (minutes)",
        y = "Relative Intensity (%)"
      ) +
      theme_bw() +
      coord_cartesian(ylim = c(0, 110)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    # Assuming gg_to_base64 is defined in your helpers.R
    output_data$results$bpc_plot_b64 <- gg_to_base64(p_bpc)
    log_message("BPC plot generation complete.")

  }, error = function(e) {
    log_message(paste("Error during BPC plot generation:", e$message))
    # This is a non-fatal error; the analysis can continue without this plot.
    # We can add a placeholder or just log the error.
  })
  # Note: We don't return on this error, allowing the next plot to still be generated.


  # --- 4. GENERATE TOP 5 SPECTRA PLOT (BASE64 ENCODED) ---
  tryCatch({
    log_message("Generating normalized static plot of top 5 spectra...")
    
    num_features_to_plot <- min(5, nrow(output_data$results$top_features))
    if (num_features_to_plot == 0) {
        stop("No features found to plot.")
    }
    top_indices <- output_data$results$top_features$spectrum_index[1:num_features_to_plot]
    top_spectra <- exp_spectra[top_indices]
    
    peaks_list <- peaksData(top_spectra)
    plot_df_list <- lapply(1:num_features_to_plot, function(i) {
      df <- as.data.frame(peaks_list[[i]])
      colnames(df) <- c("mz", "intensity")
      
      if (nrow(df) > 0 && max(df$intensity) > 0) {
        df <- df %>%
          mutate(relative_intensity = (intensity / max(intensity)) * 100) %>%
          filter(relative_intensity >= 1.0)
      } else {
        df$relative_intensity <- numeric(0)
      }

      df$label <- paste0(
        "#", i, " @ ", round(rtime(top_spectra[i])/60, 2), " min\n",
        "Precursor: ", round(precursorMz(top_spectra[i]), 4)
      )
      return(df)
    })
    combined_plot_df <- do.call(rbind, plot_df_list)
    
    p_static <- ggplot(combined_plot_df, aes(x = mz, xend = mz, y = 0, yend = relative_intensity)) +
      geom_segment(linewidth = 0.7) +
      facet_wrap(~label) + 
      labs(
        title = "Top 5 Most Intense Fragmentation Spectra", 
        x = "m/z", 
        y = "Relative Intensity (%)" 
      ) +
      theme_bw() +
      coord_cartesian(ylim = c(0, 110)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
    output_data$results$top_5_spectra_plot_b64 <- gg_to_base64(p_static)
    log_message("Static plot generation complete.")

  }, error = function(e) {
    log_message(paste("Error during static plot generation:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
  })

  if (!is.null(output_data$error)) return(output_data)

  # --- 5. FINALIZE AND RETURN ---
  if (is.null(output_data$error)) {
    output_data$status <- "success"
  }
  
  return(output_data)
}