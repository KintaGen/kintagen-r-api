suppressPackageStartupMessages({
  library(drc)
  library(jsonlite)
  library(ggplot2)
})
perform_drc_analysis <- function(input_csv_string) {
  # --- 1. Initialize the final output list and helpers ---

  source("./scripts/helpers.R")


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