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