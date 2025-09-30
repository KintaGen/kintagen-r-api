# This is the main API router file.
# Its only job is to define endpoints and delegate work to other scripts.

# 1. Load libraries required for the API itself
suppressPackageStartupMessages({
  library(plumber)
})

source("./scripts/drc_analysis.R")

source("./scripts/xcms_analysis.R")

source("./scripts/nmr1d_analysis.R")




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
  suppressPackageStartupMessages({
    library(drc)
    library(jsonlite)
    library(ggplot2)
  })
  # Get the raw CSV string from the request body
  input_data <- req$body
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

  suppressPackageStartupMessages({
    library(drc)
    library(jsonlite)
    library(ggplot2)
    library(Rnmr1D)
    library(archive)
  })
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


#* XC-MS Analysis Endpoint
#* An Endpoint to analyze a user-uploaded mzML file.

#* @parser multi
#* @post /analyze/xcms
#* @serializer unboxedJSON
function(req) {
  # Plumber receives the file and saves it to a temporary location.
  # The 'mzml_file' argument is a list containing info about the upload.
  # The most important part is '$value', which is the path to the temp file.
  suppressPackageStartupMessages({
    library(jsonlite)
    library(ggplot2)
    library(Spectra)
    library(dplyr)
    library(plotly)
    library(mzR)
  })
  mzml_file <- req$body$file
  if (is.null(mzml_file)) {
    return(list(error = "No file was uploaded. Please include a file named 'mzml_file'."))
  }
  temp_mzml_file_path <- tempfile(fileext = ".mzML")
  
  writeBin(mzml_file$value, temp_mzml_file_path)  # Call the main analysis function with the path to the temporary file
  results <- perform_ms_analysis(temp_mzml_file_path)
  
  # Plumber will automatically clean up the temporary file after the request is finished.
  return(results)
}

