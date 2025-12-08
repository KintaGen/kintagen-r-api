suppressPackageStartupMessages({
  library(base64enc)
})

# --- 1. INITIALIZE OUTPUT AND HELPERS ---
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

# --- Function to Get Package Citations ---
get_package_citations <- function(packages) {
  
  # Helper function to format a single citation
  format_citation <- function(pkg) {
    citation_text <- tryCatch({
      # Capture the output of the print command for the citation
      paste(capture.output(print(citation(pkg), style = "text")), collapse = "\n")
    }, error = function(e) {
      paste("Could not retrieve citation for package:", pkg)
    })
    
    version_text <- tryCatch({
      as.character(packageVersion(pkg))
    }, error = function(e) {
      "Version not found"
    })
    
    list(
      package = pkg,
      version = version_text,
      citation = citation_text
    )
  }
  
  # Apply the formatting to each package
  citation_list <- lapply(packages, format_citation)
  return(citation_list)
}