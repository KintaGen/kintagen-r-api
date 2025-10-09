# identifier.R
# This script is STATELESS. It takes the JSON output from the feature finder
# (which includes embedded raw spectra) and a library file, and performs
# spectral matching and generates validation plots.

# Suppress package startup messages for cleaner logs
suppressPackageStartupMessages({
  library(jsonlite)
  library(Spectra)
  library(dplyr)
  library(MetaboAnnotation)
  library(MsBackendMsp) 
  library(httr)
})

# Make the source path robust to where the script is called from
source('./scripts/helpers.R')

fetch_candidates_from_mona <- function(mz, tolerance = 0.5) {
  mz_low <- mz - tolerance
  mz_high <- mz + tolerance
  
  # This is the smarter query string. It filters by mass.
  query_string <- paste0(
    "exists(compound.metaData.name:'total exact mass' and compound.metaData.value>:", mz_low,
    " and compound.metaData.value<:", mz_high, ")"
  )
  
  encoded_query <- URLencode(query_string)
  
  # --- STEP 1: Get the total count of candidates ---
  count_api_url <- paste0(
    "https://mona.fiehnlab.ucdavis.edu/rest/spectra/search/count?query=",
    encoded_query
  )
  
  log_message(paste("Querying MoNA API for candidate count for m/z ~", round(mz, 2)))
  
  count_response <- GET(count_api_url)
  if (http_status(count_response)$category != "Success") {
    log_message("MoNA count API request failed.")
    return(NULL)
  }
  
  total_candidates <- as.integer(content(count_response, "text", encoding = "UTF-8"))
  
  if (is.na(total_candidates) || total_candidates == 0) {
    log_message("MoNA API returned 0 candidates for this mass.")
    return(NULL)
  }
  
  log_message(paste("Found", total_candidates, "total candidates. Fetching all of them in pages..."))
  
  # --- STEP 2: Paginate through all results ---
  # Define a page size. 500 is a reasonable number to be efficient without overloading the server.
  page_size <- 500 
  from_index <- 0
  master_list_of_dfs <- list() # To store results from all pages
  total_pages <- ceiling(total_candidates / page_size)
  
  while (from_index < total_candidates) {
    current_page <- (from_index / page_size) + 1
    log_message(paste("Fetching page", current_page, "of", total_pages, "... (records", from_index + 1, "to", min(from_index + page_size, total_candidates), ")"))

    # Construct the API URL for the current page using 'from' and 'size'
    api_url <- paste0(
      "https://mona.fiehnlab.ucdavis.edu/rest/spectra/search?query=",
      encoded_query,
      "&size=", page_size,
      "&from=", from_index
    )
    
    response <- GET(api_url)
    if (http_status(response)$category != "Success") {
      log_message(paste("API request failed for page", current_page, ". Aborting fetch."))
      # We return NULL here because the data would be incomplete and potentially misleading.
      return(NULL)
    }
    
    api_results_df <- fromJSON(content(response, "text", encoding = "UTF-8"))
    
    if (is.null(api_results_df) || nrow(api_results_df) == 0) {
      log_message(paste("Received an empty page (page", current_page, ") unexpectedly. Stopping fetch."))
      break # Exit the loop if a page is empty
    }
    
    # --- STEP 3: Process the results from the current page ---
    for (i in 1:nrow(api_results_df)) {
      hit <- api_results_df[i, ]
      
      if (is.null(hit$spectrum) || !is.character(hit$spectrum) || nchar(hit$spectrum) == 0) {
        next 
      }
      
      peaks <- strsplit(hit$spectrum, " ")[[1]]
      pks_matrix <- do.call(rbind, lapply(peaks, function(p) as.numeric(strsplit(p, ":")[[1]])))
      if (is.null(pks_matrix) || nrow(pks_matrix) == 0) next
      colnames(pks_matrix) <- c("mz", "intensity")
      
      compound_name <- "Unknown"
      if (length(hit$compound[[1]]) > 0 && length(hit$compound[[1]]$names[[1]]) > 0) {
        compound_name <- hit$compound[[1]]$names[[1]]$name[1]
      }
      
      # Add the processed DataFrame to our master list
      master_list_of_dfs[[length(master_list_of_dfs) + 1]] <- DataFrame(msLevel = 2L, name = compound_name, pks = list(pks_matrix))
    }
    
    # Increment 'from' for the next page
    from_index <- from_index + page_size
  } # End of while loop
  
  
  # --- STEP 4: Combine all pages and return the final Spectra object ---
  if (length(master_list_of_dfs) == 0) {
    log_message("Processing completed, but no valid spectrum data was found across all pages.")
    return(NULL)
  }
  
  log_message(paste("Successfully fetched and processed", length(master_list_of_dfs), "valid spectra from MoNA."))
  combined_df <- do.call(rbind, master_list_of_dfs)
  return(Spectra(combined_df, backend = MsBackendDataFrame()))
}


perform_identification <- function(input_json_path) {
  
  output_data <- list(status = "processing", results = list(), error = NULL)
  
  tryCatch({
    log_message("Reading single spectrum data from JSON...")
    # The input JSON is the simple object for one peak
    spec_info <- fromJSON(input_json_path, simplifyDataFrame = TRUE)
    
    if (is.null(spec_info$peak_number) || is.null(spec_info$spectrum_data)) {
      stop("Input JSON is missing required 'peak_number' or 'spectrum_data'.")
    }

    peak_num <- spec_info$peak_number 
    

    query_pks_matrix <- as.matrix(as.data.frame(spec_info$spectrum_data))
    
    if (nrow(query_pks_matrix) == 0) {
      final_match <- data.frame(peak_number = peak_num, match_name = "No Peaks in Query", similarity_score = 0)
    } else {
      # The UI sends 'relative_intensity', but matching functions need 'intensity'.
      colnames(query_pks_matrix) <- c("mz", "intensity")
      
      molecular_ion_mz <- max(query_pks_matrix[, "mz"])
      candidate_lib <- fetch_candidates_from_mona(molecular_ion_mz)
      
      if (is.null(candidate_lib)) {
        final_match <- data.frame(peak_number = peak_num, match_name = "No API Candidates", similarity_score = 0)
      } else {
        query_spec <- Spectra(S4Vectors::DataFrame(pks = list(query_pks_matrix)), backend = MsBackendDataFrame())
        match_params <- MatchForwardReverseParam(requirePrecursor = FALSE, tolerance = 0.5, ppm = 0)
        matches <- matchSpectra(query_spec, candidate_lib, param = match_params)
        all_hits_df <- as.data.frame(matchedData(matches))
        
        if (any(!is.na(all_hits_df$score))) {
          best_hit <- all_hits_df %>% filter(!is.na(score)) %>% slice_max(order_by = score, n = 1)
          final_match <- data.frame(
            peak_number = peak_num,
            match_name = best_hit$target_name,
            similarity_score = round(best_hit$score, 3)
          )
        } else {
          final_match <- data.frame(peak_number = peak_num, match_name = "No Match Found", similarity_score = 0)
        }
      }
    }
    
    output_data$results$library_matches <- final_match
    output_data$status <- "success"
    log_message("Single peak identification complete.")

  }, error = function(e) {
    output_data$status <<- "error"; output_data$error <<- e$message
  })
  
  return(output_data)
}

# --- SCRIPT EXECUTION ---
if (!interactive()) {
  tryCatch({
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2) {
      stop("Usage: Rscript identifier.R <input_results.json> <output_matches.json>")
    }
    input_json <- args[1]; output_json <- args[2]
    cat(input_json)
    result <- perform_identification(input_json_path = input_json)
    json_output <- toJSON(result, auto_unbox = TRUE, pretty = TRUE)
    write(json_output, file = output_json)
  }, error = function(e) {
    error_json <- toJSON(list(status = "error", error = e), auto_unbox = TRUE)
    try(write(error_json, file = commandArgs(trailingOnly=TRUE)[3]), silent = TRUE)
    quit(status = 1)
  })
}