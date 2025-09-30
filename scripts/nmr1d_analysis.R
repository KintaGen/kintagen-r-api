perform_nmr_analysis <- function(zip_file_path) {

  # --- Initialize the final output list and helpers ---
  source("./scripts/helpers.R")

  # --- Local variables for cleanup ---
  temp_zip_path <- NULL
  temp_unzip_dir <- NULL

  # --- 1. HANDLE INPUT & DECODE/UNZIP DATA ---
  tryCatch({
    if (is.null(zip_file_path) || !file.exists(zip_file_path)) {
      stop("No input ZIP file was provided or found at the temporary path.")
    }
    cat(zip_file_path)
    temp_unzip_dir <- tempfile(pattern = "nmr_unzipped_")
    dir.create(temp_unzip_dir)
    log_message(paste("Unzipping data from", zip_file_path, "into", temp_unzip_dir))
    # Use the robust archive_extract function.
    archive::archive_extract(archive = zip_file_path, dir = temp_unzip_dir)
    
    all_contents <- list.files(temp_unzip_dir, full.names = TRUE)
    
    if (file.exists(file.path(temp_unzip_dir, "fid"))) {
      varian_sample_path <- temp_unzip_dir
      log_message("Found Varian files at the root of the ZIP.")
    } else if (length(all_contents) == 1 && dir.exists(all_contents[1])) {
      varian_sample_path <- all_contents[1]
      log_message(paste("Found Varian data in a single sub-directory:", basename(varian_sample_path)))
    } else {
      stop("Could not determine the Varian data path in the ZIP file.")
    }
    cat(all_contents)
  }, error = function(e) {
    output_data$error <<- paste("Input Error:", e$message)
  })
  
  # --- 2. NMR DATA PROCESSING (Only run if input was successful) ---
  if (is.null(output_data$error)) {
    tryCatch({
      log_message("Processing raw Varian data using Rnmr1D...")
      procParams <- Spec1rProcpar
      procParams$VENDOR <- 'varian'
      procParams$INPUT_SIGNAL <- 'fid'
      procParams$LB <- 0.3
      procParams$ZEROFILLING <- TRUE
      procParams$ZFFAC <- 2
      procParams$OPTPHC1 <- TRUE
      procParams$TSP <- TRUE
      procParams$LOGFILE <- ""
      
      spec_object <- Spec1rDoProc(Input = varian_sample_path, param = procParams)
      
      spectrum_df <- data.frame(PPM = spec_object$ppm, Intensity = spec_object$int)
      output_data$results$spectrum_data <- spectrum_df
      log_message("NMR processing complete.")
      
    }, error = function(e) {
      output_data$error <<- paste("NMR Processing Error:", e$message)
    })
  }

  # --- 3. GENERATE PLOT (Only run if processing was successful) ---
  if (is.null(output_data$error)) {
    tryCatch({
      log_message("Generating plot...")
      p_nmr <- ggplot(data = spectrum_df, aes(x = PPM, y = Intensity)) +
        geom_line(color = "#00529B", linewidth = 0.7) +
        scale_x_reverse(name = "Chemical Shift (ppm)") +
        labs(title = "Processed 1H NMR Spectrum", y = "Intensity") +
        theme_bw()
      
      output_data$results$plot_b64 <- gg_to_base64(p_nmr)
      log_message("Plot generation complete.")
      
    }, error = function(e) {
      output_data$error <<- paste("Plot Generation Error:", e$message)
    })
  }

  # --- 4. FINALIZE AND RETURN ---
  if (!is.null(temp_zip_path)) unlink(temp_zip_path)
  if (!is.null(temp_unzip_dir)) unlink(temp_unzip_dir, recursive = TRUE)
  
  output_data$status <- ifelse(is.null(output_data$error), "success", "error")
  return(output_data)
}