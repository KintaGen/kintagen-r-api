# identifier.R
# This script takes the output from feature_finder.R and a library file,
# and performs only the spectral matching step.

suppressPackageStartupMessages({
  library(jsonlite)
  library(Spectra)
  library(dplyr)
  library(MetaboAnnotation)
  library(MsBackendMsp)
  library(ggplot2)    # Now needed for plotting
  library(patchwork)  # For combining plots
})

source("./scripts/helpers.R")

perform_identification <- function(input_json_path, library_path) {
  
  output_data <- list(status = "processing", results = list(), error = NULL)
  chunk_size <- 200
  
  tryCatch({
    # Step 1: Read the JSON data from the first script
    log_message("Reading initial analysis data from JSON...")
    initial_data <- fromJSON(input_json_path)
    mzml_path <- initial_data$mzMl_file_path
    top_features <- initial_data$results$top_features
    if (is.null(mzml_path) || is.null(top_features)) {
      stop("Input JSON is missing required 'mzMl_file_path' or 'top_features' data.")
    }
    
    # Step 2: Load the necessary data
    log_message("Loading original spectra from .mzML file...")
    exp_spectra <- Spectra(mzml_path, backend = MsBackendMzR())
    
    log_message(paste("Loading reference library from:", library_path))
    ref_lib <- Spectra(
      readMsp(f = library_path, skip = 0, n = chunk_size), 
      backend = MsBackendDataFrame()
    )
    
    # Step 3: Perform the matching
    log_message("Performing spectral matching...")
    initial_query_spectra <- exp_spectra[top_features$spectrum_index]
    initial_query_spectra$peak_number <- top_features$peak_number
    
    non_empty_mask <- initial_query_spectra$peaksCount > 0
    if (!any(non_empty_mask)) { stop("All selected spectra are empty.") }
    
    query_spectra_for_matching <- initial_query_spectra[non_empty_mask]
    top_features_for_matching <- top_features[non_empty_mask, ]
    
    match_params <- MatchForwardReverseParam(requirePrecursor = FALSE, tolerance = 0.5, ppm = 0)
    matches <- matchSpectra(query_spectra_for_matching, ref_lib, param = match_params)
    
    # --- START: FINAL, CORRECTED RESULT PROCESSING ---
    log_message("Processing match results...")
    
    all_hits_df <- matchedData(matches)
    all_hits_df <- as.data.frame(all_hits_df)
    # Check if there were any matches at all
    if(any(!is.na(all_hits_df$score))) {
      successful_hits_df <- all_hits_df %>% filter(!is.na(score))
      successful_hits_df$target_idx <- targetIndex(matches)
      best_hits_df <- successful_hits_df %>%
        group_by(peak_number) %>% slice_max(order_by = score, n = 1, with_ties = FALSE) %>% ungroup() %>%
        transmute(peak_number = as.integer(peak_number), match_name = target_Name, similarity_score = round(score, 3), target_idx = target_idx)
      final_matches <- top_features_for_matching %>% select(peak_number) %>%
        left_join(best_hits_df, by = "peak_number") %>%
        mutate(match_name = ifelse(is.na(match_name), "No match found", match_name), similarity_score = ifelse(is.na(similarity_score), 0, similarity_score))
    } else {
      final_matches <- top_features_for_matching %>% select(peak_number) %>% mutate(match_name = "No match found", similarity_score = 0)
    }
    output_data$results$library_matches <- final_matches
    
    # --- Step 5: Generate Stacked Plots for Confident Matches ---
    confident_hits <- final_matches %>% filter(similarity_score > 0.7)
    if (nrow(confident_hits) > 0) {
      log_message(paste("Generating", nrow(confident_hits), "stacked validation plots..."))
      plot_list <- list()
      
      for (i in 1:nrow(confident_hits)) {
        hit <- confident_hits[i, ]
        hit_details <- best_hits_df %>% filter(peak_number == hit$peak_number)
        
        if (nrow(hit_details) > 0 && !is.na(hit_details$target_idx)) {
          query_spec <- query_spectra_for_matching[query_spectra_for_matching$peak_number == hit$peak_number]
          target_spec_lazy <- ref_lib[hit_details$target_idx]
          target_spec <- setBackend(target_spec_lazy, MsBackendDataFrame())
          
          query_peaks_list <- peaksData(query_spec)
          target_peaks_list <- peaksData(target_spec)
          
          if (length(query_peaks_list) > 0 && nrow(query_peaks_list[[1]]) > 0 && 
              length(target_peaks_list) > 0 && nrow(target_peaks_list[[1]]) > 0) {
            
            # --- START: NEW PUBLICATION-QUALITY PLOT LOGIC ---
            all_query_peaks <- as.data.frame(query_peaks_list[[1]]) %>% mutate(intensity = 100 * intensity / max(intensity, na.rm = TRUE))
            all_target_peaks <- as.data.frame(target_peaks_list[[1]]) %>% mutate(intensity = 100 * intensity / max(intensity, na.rm = TRUE))
            
            matched_pairs <- joinPeaks(query_peaks_list[[1]], target_peaks_list[[1]], tolerance = 0.5, ppm = 0)
            
            matched_query_peaks <- data.frame(mz=numeric(), intensity=numeric())
            matched_target_peaks <- data.frame(mz=numeric(), intensity=numeric())
            
            if (is.matrix(matched_pairs) && nrow(matched_pairs) > 0) {
              matched_query_peaks <- data.frame(mz = matched_pairs[, "x.mz"], intensity = 100 * matched_pairs[, "x.intensity"] / max(query_peaks_list[[1]][,2], na.rm=TRUE))
              matched_target_peaks <- data.frame(mz = matched_pairs[, "y.mz"], intensity = 100 * matched_pairs[, "y.intensity"] / max(target_peaks_list[[1]][,2], na.rm=TRUE))
            }
            
            # --- Plot 1: Query (Experimental) Spectrum ---
            p_query <- ggplot() +
              # Draw all peaks in gray first as a background
              geom_segment(data = all_query_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "grey80", linewidth = 0.8) +
              # Draw the matched peaks on top in blue
              geom_segment(data = matched_query_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "steelblue", linewidth = 1.2) +
              # Add text labels for matched peaks
              geom_text(data = matched_query_peaks, aes(x = mz, y = intensity, label = round(mz, 1)), angle = 45, hjust = -0.1, vjust = -0.2, size = 2.8) +
              labs(title = paste("Experimental Spectrum (Peak #", hit$peak_number, ")"), x = NULL, y = "Rel. Intensity (%)") +
              coord_cartesian(ylim = c(0, 25), expand = FALSE) + # More headroom for labels
              theme_bw() +
              theme(plot.title = element_text(face = "bold"))
            
            # --- Plot 2: Library (Reference) Spectrum ---
            p_library <- ggplot() +
              geom_segment(data = all_target_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "grey80", linewidth = 0.8) +
              geom_segment(data = matched_target_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "darkred", linewidth = 1.2) +
              geom_text(data = matched_target_peaks, aes(x = mz, y = intensity, label = round(mz, 1)), angle = 45, hjust = -0.1, vjust = -0.2, size = 2.8) +
              labs(title = paste("Library Reference:", hit$match_name), x = "m/z", y = "Rel. Intensity (%)") +
              coord_cartesian(ylim = c(0, 10), expand = FALSE) +
              theme_bw() +
              theme(plot.title = element_text(face = "bold"))
            
            # --- Combine Plots with Patchwork ---
            combined_plot_for_hit <- p_query / p_library + 
              plot_annotation(
                title = "Spectral Match Validation",
                subtitle = paste("Similarity Score:", hit$similarity_score)
              ) & theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
            
            plot_list[[i]] <- combined_plot_for_hit
            # --- END: NEW PUBLICATION-QUALITY PLOT LOGIC ---
          }
        }
      }
      
      if (length(plot_list) > 0) {
        final_plot <- wrap_plots(plot_list, ncol = 1)
        output_data$results$match_plots_b64 <- gg_to_base64(final_plot, height = 7 * length(plot_list), width = 10)
      }
    } else {
      log_message("No confident matches to plot.")
    }

    log_message("Library matching complete.")
    
    output_data$results$library_matches <- final_matches
    # --- END: FINAL, CORRECTED RESULT PROCESSING ---
    
    # --- Step 5: Generate Mirror Plots for Confident Matches ---
    confident_hits <- final_matches %>% filter(similarity_score > 0.7)
    
    for (i in 1:nrow(confident_hits)) {
      hit <- confident_hits[i, ]
      hit_details <- best_hits_df %>% filter(peak_number == hit$peak_number)
      
      if (nrow(confident_hits) > 0) {
        log_message(paste("Generating", nrow(confident_hits), "overlaid plots..."))
        plot_list <- list()
        
        for (i in 1:nrow(confident_hits)) {
          hit <- confident_hits[i, ]
          hit_details <- best_hits_df %>% filter(peak_number == hit$peak_number)
          
          if (nrow(hit_details) > 0 && !is.na(hit_details$target_idx)) {
            query_spec <- query_spectra_for_matching[query_spectra_for_matching$peak_number == hit$peak_number]
            target_spec_lazy <- ref_lib[hit_details$target_idx]
            target_spec <- setBackend(target_spec_lazy, MsBackendDataFrame())
            
            query_peaks_list <- peaksData(query_spec)
            target_peaks_list <- peaksData(target_spec)
            
            if (length(query_peaks_list) > 0 && nrow(query_peaks_list[[1]]) > 0 && 
                length(target_peaks_list) > 0 && nrow(target_peaks_list[[1]]) > 0) {
              
              # --- START: NEW OVERLAID PLOT LOGIC ---
              query_peaks <- as.data.frame(query_peaks_list[[1]]) %>%
                mutate(intensity = 100 * intensity / max(intensity, na.rm = TRUE), source = "Query (Experiment)")
              
              # Both spectra now have POSITIVE intensity
              target_peaks <- as.data.frame(target_peaks_list[[1]]) %>%
                mutate(intensity = 100 * intensity / max(intensity, na.rm = TRUE), source = "Library (Reference)")
              
              combined_peaks <- rbind(query_peaks, target_peaks)
              
              # --- Plot 1: Query Spectrum ---
              p_query <- ggplot() +
                geom_segment(data = all_query_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "gray80", linewidth = 0.8) +
                geom_segment(data = matched_query_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "steelblue", linewidth = 1.2) +
                geom_text(data = matched_query_peaks, aes(x = mz, y = intensity, label = round(mz, 2)), angle = 60, hjust = 0, vjust = -0.2, size = 2.5) +
                labs(title = "Query Spectrum (Experiment)", x = NULL, y = "Rel. Intensity (%)") +
                coord_cartesian(ylim = c(0, 120)) +
                theme_bw()
              
              # --- Plot 2: Library Spectrum ---
              p_library <- ggplot() +
                geom_segment(data = all_target_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "gray80", linewidth = 0.8) +
                geom_segment(data = matched_target_peaks, aes(x = mz, y = 0, xend = mz, yend = intensity), color = "darkred", linewidth = 1.2) +
                geom_text(data = matched_target_peaks, aes(x = mz, y = intensity, label = round(mz, 2)), angle = 60, hjust = 0, vjust = -0.2, size = 2.5) +
                labs(title = "Library Spectrum (Reference)", x = "m/z", y = "Rel. Intensity (%)") +
                coord_cartesian(ylim = c(0, 120)) +
                theme_bw()
              
              # --- Combine Plots with Patchwork ---
              combined_plot_for_hit <- p_query / p_library + 
                plot_annotation(
                  title = paste("Peak #", hit$peak_number, ": ", hit$match_name),
                  subtitle = paste("Similarity Score:", hit$similarity_score)
                )
              
              plot_list[[i]] <- combined_plot_for_hit
              # --- END: NEW OVERLAID PLOT LOGIC ---
            }
          }
        } # end for loop
        
        if (length(plot_list) > 0) {
          combined_plot <- wrap_plots(plot_list, ncol = 1)
          output_data$results$match_plots_b64 <- gg_to_base64(combined_plot, height = 4 * length(plot_list), width = 8)
        }
      } else {
        log_message("No confident matches to plot.")
      }
    } # end for loop
    
    output_data$status <- "success"
    log_message("Identification complete.")
    
  }, error = function(e) {
    output_data$status <<- "error"; output_data$error <<- e$message
  })
  
  # Corrected the return variable name
  return(output_data)
}

# --- SCRIPT EXECUTION ---
if (!interactive()) {
  tryCatch({
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 3) {
      stop("Usage: Rscript identifier.R <input_results.json> <library_file> <output_matches.json>")
    }
    input_json <- args[1]; library_file <- args[2]; output_json <- args[3]
    result <- perform_identification(input_json_path = input_json, library_path = library_file)
    json_output <- toJSON(result, auto_unbox = TRUE, pretty = TRUE)
    write(json_output, file = output_json)
  }, error = function(e) {
    error_json <- toJSON(list(status = "error", error = e$message), auto_unbox = TRUE)
    try(write(error_json, file = commandArgs(trailingOnly=TRUE)[3]), silent = TRUE)
    quit(status = 1)
  })
}