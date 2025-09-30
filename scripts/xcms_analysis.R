suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(Spectra)
  library(dplyr)
  library(plotly)
  library(mzR)
})

perform_ms_analysis <- function(mzMl_file_path) {
  #--- 1. Initialize the final output list and helpers ---

  source("./scripts/helpers.R")
  
  # --- 2. FEATURE DISCOVERY ---
  tryCatch({
    log_message("Starting feature discovery...")
    exp_spectra <- Spectra(mzMl_file_path, backend = MsBackendMzR())

    spectra_summary <- data.frame(
      spectrum_index = 1:length(exp_spectra),
      precursor_mz = precursorMz(exp_spectra),
      retention_time_sec = rtime(exp_spectra),
      total_ion_count = tic(exp_spectra)
    )
    
    top_10_features <- spectra_summary %>%
      filter(total_ion_count > 0) %>%
      arrange(desc(total_ion_count)) %>%
      head(10)
      
    output_data$results$top_features <- top_10_features
    log_message(paste("Found", nrow(top_10_features), "top features."))
    
  }, error = function(e) {
    log_message(paste("Error during feature discovery:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
    return(output_data)
  })
  
  if (!is.null(output_data$error)) return(output_data)


  # --- 3. GENERATE STATIC PLOT (BASE64 ENCODED) ---
  tryCatch({
    log_message("Generating static plot of top 5 spectra...")
    
    # Get data for the top 5 spectra
    top_5_indices <- output_data$results$top_features$spectrum_index[1:5]
    top_5_spectra <- exp_spectra[top_5_indices]
    
    # Combine peak data into a single data frame for ggplot faceting
    peaks_list <- peaksData(top_5_spectra)
    plot_df_list <- lapply(1:5, function(i) {
      df <- as.data.frame(peaks_list[[i]])
      colnames(df) <- c("mz", "intensity")
      # Create a label for each plot panel
      df$label <- paste0(
        "#", i, " @ ", round(rtime(top_5_spectra[i])/60, 2), " min\n",
        "Precursor: ", round(precursorMz(top_5_spectra[i]), 4)
      )
      return(df)
    })
    combined_plot_df <- do.call(rbind, plot_df_list)
    
    # Create the faceted plot
    p_static <- ggplot(combined_plot_df, aes(x = mz, xend = mz, y = 0, yend = intensity)) +
      geom_segment(color = "grey50") +
      facet_wrap(~label, scales = "free_y") +
      labs(title = "Top 5 Most Intense Fragmentation Spectra", x = "m/z", y = "Intensity") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
    output_data$results$top_5_spectra_plot_b64 <- gg_to_base64(p_static)
    log_message("Static plot generation complete.")

  }, error = function(e) {
    log_message(paste("Error during static plot generation:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
  })

  if (!is.null(output_data$error)) return(output_data)


  # --- 4. GENERATE INTERACTIVE PLOT DATA (PLOTLY JSON) ---
  tryCatch({
    log_message("Generating interactive TIC plot data...")
    
    # Re-use the full spectra_summary from the feature discovery step
    spectra_summary <- data.frame(
      rt_sec = rtime(exp_spectra),
      tic = tic(exp_spectra)
    )
    
    # Build the ggplot object (it will be converted to JSON later)
    p_interactive <- ggplot(spectra_summary, aes(x = rt_sec / 60, y = tic)) +
      geom_line(color = "steelblue") +
      labs(
        title = "Interactive MS2 TIC",
        x = "Retention Time (minutes)",
        y = "Total Ion Count (TIC)"
      ) +
      theme_bw()

    # Convert to plotly JSON and store it
    # This JSON can be directly fed into Plotly.js in the frontend
    output_data$results$interactive_tic_plot_json <- plotly_json(p_interactive, pretty = FALSE)
    log_message("Interactive plot JSON created.")
    
  }, error = function(e) {
    log_message(paste("Error generating plotly JSON:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
  })

  # --- 5. FINALIZE AND RETURN ---
  if (is.null(output_data$error)) {
    output_data$status <- "success"
  }
  
  return(output_data)
}