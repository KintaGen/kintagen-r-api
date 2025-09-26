# This is the main API router file.
# Its only job is to define endpoints and delegate work to other scripts.

# 1. Load libraries required for the API itself
suppressPackageStartupMessages({
  library(plumber)
  library(archive)
  library(drc)
  library(jsonlite)
  library(ggplot2)
  library(base64enc)
  library(Rnmr1D)
})

perform_drc_analysis <- function(input_csv_string) {

  # --- Initialize the final output list ---
  output_data <- list(
    status = "processing",
    error = NULL,
    log = c(),
    results = list()
  )
  
  # --- Helper function for logging ---
  log_message <- function(msg) {
    output_data$log <<- c(output_data$log, msg)
  }
  
  # --- Helper function for plot encoding ---
  gg_to_base64 <- function(gg, width = 7, height = 5) {
      temp_file <- tempfile(fileext = ".png")
      ggsave(temp_file, plot = gg, width = width, height = height, dpi = 150)
      base64_string <- base64enc::base64encode(temp_file)
      unlink(temp_file)
      return(paste0("data:image/png;base64,", base64_string))
  }

  # --- 1. HANDLE INPUT DATA ---
  tryCatch({
    if (length(input_csv_string) == 0 || nchar(input_csv_string) == 0) {
      log_message("No input data provided. Using internal sample data.")
      data <- data.frame(
        dose = c(0.1, 0.5, 1, 5, 10, 20),
        total = rep(50, 6),
        response = c(1, 5, 10, 25, 40, 48)
      )
    } else {
      log_message("Reading data from provided CSV text.")
      data <- read.csv(text = input_csv_string)
    }
  }, error = function(e) {
    log_message(paste("Error reading input data:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
    return(output_data) # Exit function early on error
  })
  
  # If there was an error reading data, the error field will be set.
  # We should stop here.
  if (!is.null(output_data$error)) {
    return(output_data)
  }

  # --- 2. DOSE-RESPONSE ANALYSIS ---
  tryCatch({
    log_message("Performing dose-response modeling...")
    model <- drm(response / total ~ dose, weights = total, data = data, fct = LL.2(), type = "binomial")
    ed_results <- ED(model, 50, interval = "delta", level = 0.95, display = FALSE)
    model_summary_obj <- summary(model)
    
    output_data$results$ld50_estimate <- ed_results[1]
    output_data$results$standard_error <- ed_results[2]
    output_data$results$confidence_interval_lower <- ed_results[3]
    output_data$results$confidence_interval_upper <- ed_results[4]
    output_data$results$model_coefficients <- coef(model_summary_obj)
    
    log_message("Dose-response analysis complete.")
  }, error = function(e) {
    log_message(paste("Error during DRC modeling:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
    return(output_data) # Exit function early
  })
  
  if (!is.null(output_data$error)) {
    return(output_data)
  }

  # --- 3. GENERATE PLOT ---
  tryCatch({
    log_message("Generating plot...")
    plot_data <- data.frame(dose = data$dose, proportion = data$response / data$total)
    min_dose_nonzero <- min(plot_data$dose[plot_data$dose > 0], na.rm = TRUE)
    max_dose <- max(plot_data$dose, na.rm = TRUE)
    curve_data <- data.frame(dose = exp(seq(log(min_dose_nonzero), log(max_dose), length.out = 100)))
    curve_data$p <- predict(model, newdata = curve_data)
    ld50_val <- output_data$results$ld50_estimate
    
    p_ld50 <- ggplot(plot_data, aes(x = dose, y = proportion)) +
        geom_point(size = 3, shape = 16) +
        geom_line(data = curve_data, aes(x = dose, y = p), color = "blue", linewidth = 1) +
        annotate("point", x = ld50_val, y = 0.5, color = "red", size = 4, shape = 18) +
        annotate("segment", x = ld50_val, y = 0, xend = ld50_val, yend = 0.5, linetype = "dashed", color = "darkgrey") +
        annotate("segment", x = min_dose_nonzero, y = 0.5, xend = ld50_val, yend = 0.5, linetype = "dashed", color = "darkgrey") +
        annotate("label", x = ld50_val, y = 0.1, label = sprintf("LD50 = %.3f", ld50_val), hjust = 0, nudge_x = 0.05, fontface = "bold") +
        scale_x_log10(name = "Dose (log scale)", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(title = "Dose-Response Curve with LD50 Estimate", y = "Response Proportion") +
        annotation_logticks(sides = "b") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
    output_data$results$plot_b64 <- gg_to_base64(p_ld50)
    log_message("Plot generation complete.")
  }, error = function(e) {
    log_message(paste("Error during plot generation:", e$message))
    output_data$status <<- "error"
    output_data$error <<- e$message
  })

  # --- 4. FINALIZE AND RETURN ---
  if (is.null(output_data$error)) {
    output_data$status <- "success"
  }
  
  return(output_data)
}

# This script contains the core 1D NMR spectrum processing logic.
# It is designed to be sourced by another script (like plumber.R).

perform_nmr_analysis <- function(zip_file_path) {

  # --- Initialize the final output list and helpers ---
  output_data <- list(status = "processing", error = NULL, log = c(), results = list())
  log_message <- function(msg) { output_data$log <<- c(output_data$log, msg) }
  
  gg_to_base64 <- function(gg, width = 10, height = 6) {
    temp_file <- tempfile(fileext = ".png")
    ggsave(temp_file, plot = gg, width = width, height = height, dpi = 150)
    base64_string <- base64enc::base64encode(temp_file)
    unlink(temp_file)
    return(paste0("data:image/png;base64,", base64_string))
  }

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

# =========================================================================
#  3. CORS FILTER
#  This filter runs before any endpoint and adds the necessary headers
#  to tell the browser that requests from your UI are allowed.
# =========================================================================
#* @filter cors
function(req, res) {
  # Set the permission slip header
  res$setHeader("Access-Control-Allow-Origin", "*")

  # Handle preflight OPTIONS requests
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$setHeader("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
    res$setHeader("Access-Control-Allow-Headers", req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS)
    res$status <- 204 # No Content status
    return(list())
  }
  
  # Forward the request to the endpoint
  plumber::forward()
}

#* Liveness check to confirm the API is running
#* @get /healthcheck
function() {
  list(status = "OK", timestamp = Sys.time(),wd = getwd())
}

#* Perform Dose-Response (LD50/ED50) Analysis
#* This endpoint receives the data and passes it to the analysis function.
#* @post /analyze/drc
#* @param req The HTTP request object containing the raw CSV data in its body.
#* @serializer unboxedJSON
function(req) {
  # Get the raw CSV string from the request body
  input_data <- req$body
  cat(input_data)
  # Call the dedicated analysis function from our other script
  results <- perform_drc_analysis(input_data)
  
  # Return the results. Plumber handles the JSON conversion.
  return(results)
}

#* Perform 1D NMR Spectrum Processing
#* @parser multi
#* @post /analyze/nmr
#* @serializer unboxedJSON
function(req) {
  # Your log shows the file data is in req$body$file
  file_object <- req$body$file
  
  # A robust check
  if (is.null(file_object)) {
    return(list(status = "error", error = "No file was uploaded in a form field named 'file'."))
  }
  
  # Raw zip data
  raw_zip_content <- file_object$value
  # We manually write these raw bytes to a temporary file.
  temp_zip_file_path <- tempfile(fileext = ".zip")
  writeBin(raw_zip_content, temp_zip_file_path)
  # Now we have what our analysis function needs: a path to the ZIP file.
  results <- perform_nmr_analysis(temp_zip_file_path)
  
  # Clean up the temporary file we created.
  unlink(temp_zip_file_path)
  
  return(results)
}

