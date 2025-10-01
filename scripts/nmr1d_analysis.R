perform_nmr_analysis <- function(zip_file_path, gg_to_base64_func = NULL) {
  # --------------------------- Setup & guards ---------------------------------
  output_data <- list(status = "processing", error = NULL, results = list())
  temp_unzip_dir <- NULL

  log_message <- function(msg) cat(sprintf("%s - %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))

  need_pkgs <- c("Rnmr1D", "ggplot2", "dplyr", "tidyr")
  for (p in need_pkgs) if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is required.", p))
  `%>%` <- dplyr::`%>%`

  if (is.null(gg_to_base64_func)) {
    if (!requireNamespace("base64enc", quietly = TRUE)) stop("Provide gg_to_base64_func or install 'base64enc'.")
    gg_to_base64_func <- function(p) {
      tf <- tempfile(fileext = ".png")
      ggplot2::ggsave(tf, plot = p, width = 11, height = 6, dpi = 150)
      uri <- base64enc::dataURI(file = tf, mime = "image/png")
      unlink(tf)
      uri
    }
  }

  # ---------------------- Solvent reference (1H residuals) --------------------
  solvent_reference <- data.frame(
    stringsAsFactors = FALSE,
    solvent_name   = c(
      "DMSO-d6","CDCl3","Acetone-d6","CD3CN (Acetonitrile-d3)","CD3OD (Methanol-d4)",
      "D2O","C6D6 (Benzene-d6)","THF-d8","Toluene-d8","CD2Cl2 (Dichloromethane-d2)","DMF-d7"
    ),
    reference_ppm  = c(2.50, 7.26, 2.05, 1.94, 3.31, 4.79, 7.16, 1.72, 2.09, 5.32, 2.75),
    secondary_ppm  = c(3.33, NA,   NA,   NA,   4.87, NA,   NA,   3.58, NA,   NA,   2.92),
    tertiary_ppm   = c(NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   8.03)
  )

  # ------------------------------ Helpers -------------------------------------
  find_local_maxima <- function(y) { which(diff(sign(diff(y))) < 0) + 1L }
  peak_indices_prominent <- function(y, k = 6) {
    m <- stats::median(y, na.rm = TRUE)
    mad <- stats::mad(y, center = m, constant = 1.4826, na.rm = TRUE)
    idx <- find_local_maxima(y)
    idx[y[idx] > (m + k * mad)]
  }
  
  refine_peak_position <- function(full_spectrum_df, approx_ppm) {
    half_width_ppm <- 0.05
    peak_region <- full_spectrum_df %>%
      dplyr::filter(PPM_uncorrected > (approx_ppm - half_width_ppm) & 
                    PPM_uncorrected < (approx_ppm + half_width_ppm))
    if (nrow(peak_region) < 3) return(approx_ppm)
    stats::weighted.mean(peak_region$PPM_uncorrected, peak_region$Intensity)
  }

  find_tms_candidate <- function(significant_peaks_df) {
    ISOLATION_PPM_THRESHOLD <- 0.20
    rightmost_peak <- significant_peaks_df %>% dplyr::slice_min(PPM, n = 1, with_ties = FALSE)
    if (nrow(rightmost_peak) == 0) return(NULL)
    other_peaks <- significant_peaks_df %>% dplyr::filter(PPM != rightmost_peak$PPM)
    if (nrow(other_peaks) > 0) {
      min_dist <- min(abs(other_peaks$PPM - rightmost_peak$PPM))
      if (min_dist < ISOLATION_PPM_THRESHOLD) {
        log_message(sprintf("TMS Check: Rightmost significant peak at %.2f ppm is not isolated.", rightmost_peak$PPM))
        return(NULL)
      }
    }
    log_message(sprintf("TMS Check: Found plausible candidate at %.4f ppm.", rightmost_peak$PPM))
    return(rightmost_peak)
  }

  find_best_solvent_match <- function(candidate_peaks_df, solvent_ref_df) {
    best_match <- list(score = -1)
    ref_peaks <- tidyr::pivot_longer(
      solvent_ref_df,
      cols = c("reference_ppm", "secondary_ppm", "tertiary_ppm"),
      names_to = "type", values_to = "ppm"
    ) %>%
      dplyr::filter(!is.na(ppm)) %>%
      dplyr::mutate(weight = dplyr::case_when(
        type == "reference_ppm" ~ 1.0, type == "secondary_ppm" ~ 0.5, TRUE ~ 0.4
      ))
    for (i in 1:nrow(candidate_peaks_df)) {
      for (j in 1:nrow(ref_peaks)) {
        obs_peak <- candidate_peaks_df[i, ]; ref_peak <- ref_peaks[j, ]
        current_offset <- obs_peak$PPM - ref_peak$ppm
        score <- 0
        solvent_peaks_to_check <- ref_peaks %>% dplyr::filter(solvent_name == ref_peak$solvent_name)
        for (k in 1:nrow(solvent_peaks_to_check)) {
          sp <- solvent_peaks_to_check[k, ]
          expected_pos <- sp$ppm + current_offset
          match_candidates <- candidate_peaks_df %>%
            dplyr::mutate(dist = abs(PPM - expected_pos)) %>% dplyr::filter(dist < 0.05)
          if (nrow(match_candidates) > 0) {
            best_candidate <- match_candidates %>% dplyr::slice_min(dist, n = 1)
            score <- score + sp$weight * (1 - best_candidate$dist / 0.05)
          }
        }
        if (score > best_match$score) {
          best_match <- list(score = score, offset = current_offset, solvent = ref_peak$solvent_name)
        }
      }
    }
    if (best_match$score < 0.5) return(NULL)
    primary_ref_ppm <- solvent_ref_df$reference_ppm[solvent_ref_df$solvent_name == best_match$solvent]
    best_match$result <- list(
      detected_solvent = best_match$solvent,
      expected_primary_ppm = primary_ref_ppm,
      found_primary_ppm = primary_ref_ppm + best_match$offset,
      ppm_deviation = best_match$offset
    )
    return(best_match)
  }
  
  detect_solvent_from_calibrated_spectrum <- function(calibrated_peaks_df, solvent_ref_df) {
    best_solvent <- list(name = "Unknown", score = 0)
    for (i in 1:nrow(solvent_ref_df)) {
      s <- solvent_ref_df[i, ]
      target_ppm <- s$reference_ppm
      match_candidates <- calibrated_peaks_df %>%
        dplyr::filter(PPM > target_ppm - 0.05 & PPM < target_ppm + 0.05)
      if(nrow(match_candidates) > 0) {
        score <- max(match_candidates$Intensity)
        if (score > best_solvent$score) {
          best_solvent <- list(name = s$solvent_name, score = score)
        }
      }
    }
    if (best_solvent$score > 0) return(best_solvent$name)
    return("Unknown")
  }
  
  detect_vendor <- function(root_dir) {
    all <- list.files(root_dir, full.names = TRUE, recursive = TRUE)
    if (any(basename(all) == "fid") && any(basename(all) == "procpar"))
      return(list(vendor="varian", path=dirname(all[basename(all)=="fid"][1])))
    fid_candidates <- all[basename(all) == "fid"]
    if (length(fid_candidates)) {
      for (f in fid_candidates) {
        if (file.exists(file.path(dirname(f), "acqus")) || file.exists(file.path(dirname(dirname(f)), "acqus")))
          return(list(vendor="bruker", path=dirname(f)))
      }
      return(list(vendor="varian", path=dirname(fid_candidates[1])))
    }
    NULL
  }
  summarize_bands <- function(ppm, intensity) {
    bands <- list(
      "Aliphatic (0.5–1.8 ppm)" = c(0.5, 1.8),
      "Benzylic/Allylic/α-heteroatom (2–3 ppm)" = c(2.0, 3.0),
      "C–O/C–N (3.2–4.2 ppm)" = c(3.2, 4.2),
      "Water/HOD or vinylic (4.5–5.0 ppm)" = c(4.5, 5.0),
      "Aromatic (6–8.5 ppm)" = c(6.0, 8.5),
      "Aldehyde (9–10.5 ppm)" = c(9.0, 10.5)
    )
    idx <- peak_indices_prominent(intensity, k = 6)
    pk_ppm <- ppm[idx]
    rows <- lapply(names(bands), function(nm) {
      lo <- bands[[nm]][1]; hi <- bands[[nm]][2]
      mask <- pk_ppm >= lo & pk_ppm <= hi
      ex <- if (any(mask)) paste(head(sort(round(pk_ppm[mask], 2)), 6), collapse = ", ") else ""
      data.frame(Band = nm, ppm_min = lo, ppm_max = hi, count = sum(mask), examples_ppm = ex, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)[order(-sapply(rows, `[[`, "count")), , drop = FALSE]
  }

  # --- Tunable Parameters ---
  top_n_peaks_to_consider <- 30
  MIN_PEAK_INTENSITY_RATIO <- 0.01

  # ------------------------------ Main block ----------------------------------
  tryCatch({
    # ---- Unzip, Vendor Detection, Rnmr1D Processing ----
    if (is.null(zip_file_path) || !file.exists(zip_file_path)) stop("No input ZIP file provided or found.")
    temp_unzip_dir <- tempfile("nmr_unzipped_"); dir.create(temp_unzip_dir, showWarnings = FALSE)
    if (requireNamespace("archive", quietly = TRUE)) archive::archive_extract(archive = zip_file_path, dir = temp_unzip_dir)
    else utils::unzip(zip_file_path, exdir = temp_unzip_dir)
    log_message(sprintf("Unzipped into: %s", temp_unzip_dir))

    v <- detect_vendor(temp_unzip_dir)
    if (is.null(v)) stop("Could not detect vendor or find a 'fid' file in the archive.")
    log_message(sprintf("Detected vendor: %s | sample path: %s", v$vendor, v$path))

    log_message("Processing raw 1D data with Rnmr1D...")
    prm <- Rnmr1D::Spec1rProcpar
    prm$VENDOR <- v$vendor; prm$INPUT_SIGNAL <- "fid"
    prm$LB <- 0.3; prm$ZEROFILLING <- TRUE; prm$ZFFAC <- 2; prm$OPTPHC1 <- TRUE
    prm$TSP <- FALSE; prm$LOGFILE <- ""
    spec <- Rnmr1D::Spec1rDoProc(Input = v$path, param = prm)
    spectrum_df <- data.frame(PPM_uncorrected = spec$ppm, Intensity = spec$int)
    if (!nrow(spectrum_df) || all(!is.finite(spectrum_df$Intensity))) stop("Empty or invalid spectrum after processing.")
    log_message("Initial NMR processing complete.")

    # ---- Calibration and Solvent Detection Block -----------------------------
    referencing_successful <- FALSE
    ppm_correction <- 0.0
    calibration_standard <- "None"
    detected_solvent <- "Unknown"
    referencing_info_out <- list(status = "failed", message = "Could not find a high-confidence reference peak.")

    tryCatch({
      # --- Step 1: Find all significant peaks for analysis ---
      log_message("Finding significant peaks for analysis...")
      idx_pk <- peak_indices_prominent(spectrum_df$Intensity, k = 6)
      if (length(idx_pk) < 1) stop("No prominent peaks found.")
      
      all_prominent_peaks_df <- data.frame(
          PPM = spectrum_df$PPM_uncorrected[idx_pk],
          Intensity = spectrum_df$Intensity[idx_pk]
      )
      
      max_intensity <- max(all_prominent_peaks_df$Intensity, na.rm = TRUE)
      significant_peaks_df <- all_prominent_peaks_df %>%
        dplyr::filter(Intensity >= max_intensity * MIN_PEAK_INTENSITY_RATIO)
      
      if (nrow(significant_peaks_df) < 1) stop("No visually significant peaks found after filtering.")
      log_message(sprintf("Found %d visually significant peaks.", nrow(significant_peaks_df)))

      # --- Step 2: Determine Calibration Standard (TMS first) ---
      calibration_ref <- NULL
      tms_cand <- find_tms_candidate(significant_peaks_df)
      
      if (!is.null(tms_cand)) {
        log_message("Using plausible TMS candidate for calibration.")
        refined_ppm <- refine_peak_position(spectrum_df, tms_cand$PPM)
        calibration_standard <- "TMS"
        ppm_correction <- -refined_ppm
        calibration_ref <- list(expected = 0.0, found = refined_ppm)
      } else {
        log_message("TMS check failed. Using solvent matching for calibration.")
        solvent_match <- find_best_solvent_match(significant_peaks_df, solvent_reference)
        if (!is.null(solvent_match)) {
          refined_ppm <- refine_peak_position(spectrum_df, solvent_match$result$found_primary_ppm)
          calibration_standard <- solvent_match$result$detected_solvent
          ppm_correction <- -(refined_ppm - solvent_match$result$expected_primary_ppm)
          calibration_ref <- list(expected = solvent_match$result$expected_primary_ppm, found = refined_ppm)
          detected_solvent <- calibration_standard # If we calibrate on a solvent, we've also detected it.
        } else {
          stop("Could not determine a reliable calibration standard (neither TMS nor solvent).")
        }
      }
      
      referencing_successful <- TRUE
      log_message(sprintf("Calibration successful using %s. Applying correction of %.4f ppm.", calibration_standard, ppm_correction))
      
      # --- Step 3: Apply calibration and then detect bulk solvent (if needed) ---
      spectrum_df$PPM <- spectrum_df$PPM_uncorrected + ppm_correction
      
      if (calibration_standard == "TMS") {
        log_message("Calibrated with TMS. Now detecting bulk solvent...")
        calibrated_peaks_df <- significant_peaks_df
        calibrated_peaks_df$PPM <- calibrated_peaks_df$PPM + ppm_correction
        detected_solvent <- detect_solvent_from_calibrated_spectrum(calibrated_peaks_df, solvent_reference)
        log_message(sprintf("Detected bulk solvent: %s", detected_solvent))
      }
      
      referencing_info_out <- list(
        status = "success",
        calibration_standard = calibration_standard,
        expected_ppm = calibration_ref$expected,
        found_peak_at_ppm = calibration_ref$found,
        ppm_correction_applied = ppm_correction,
        detected_solvent = detected_solvent
      )
      
    }, error = function(e) {
      log_message(paste("WARNING: Automatic referencing failed:", e$message))
      referencing_info_out <- list(status = "failed", message = e$message)
    })
    
    # ---- Apply Correction and Generate All Outputs ----
    if(!("PPM" %in% names(spectrum_df))) spectrum_df$PPM <- spectrum_df$PPM_uncorrected
    
    output_data$results$referencing_info <- referencing_info_out
    output_data$results$spectrum_data <- spectrum_df %>% dplyr::select(PPM, Intensity)

    plot_subtitle <- if (referencing_successful) {
      sprintf("Calibrated to %s at %.2f ppm (Δ=%.4f ppm). Detected Solvent: %s.",
              referencing_info_out$calibration_standard, referencing_info_out$expected_ppm, 
              referencing_info_out$ppm_correction_applied, referencing_info_out$detected_solvent)
    } else {
      "WARNING: Automatic referencing failed. PPM scale may be inaccurate."
    }
    
    p_main <- ggplot2::ggplot(spectrum_df, ggplot2::aes(PPM, Intensity)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::scale_x_reverse(name = "Chemical Shift (ppm)") +
      ggplot2::labs(title = "Processed 1H NMR Spectrum", subtitle = plot_subtitle, y = "Intensity") +
      ggplot2::theme_bw()
    output_data$results$plot_b64 <- gg_to_base64_func(p_main)

    if (referencing_successful) {
      zoom_half <- 0.5
      expected_ppm <- referencing_info_out$expected_ppm
      df_zoom <- subset(spectrum_df, PPM >= (expected_ppm - zoom_half) & PPM <= (expected_ppm + zoom_half))
      if(nrow(df_zoom) > 2) {
          p_zoom <- ggplot2::ggplot(df_zoom, ggplot2::aes(PPM, Intensity)) +
            ggplot2::geom_line(linewidth = 0.7) +
            ggplot2::geom_vline(xintercept = expected_ppm, linetype = 2) +
            ggplot2::scale_x_reverse(name = sprintf("Region around %.2f ppm", expected_ppm)) +
            ggplot2::labs(title = sprintf("Calibration Standard Alignment: %s", referencing_info_out$calibration_standard), y = "Intensity") +
            ggplot2::theme_bw()
          output_data$results$residual_zoom_plot_b64 <- gg_to_base64_func(p_zoom)
      }
    }

    final_peaks_idx <- which(spectrum_df$PPM_uncorrected %in% significant_peaks_df$PPM)
    peaks_tbl <- data.frame(PPM = round(spectrum_df$PPM[final_peaks_idx], 4),
                            Intensity = spectrum_df$Intensity[final_peaks_idx], stringsAsFactors = FALSE)
    output_data$results$peaks <- head(peaks_tbl[order(-peaks_tbl$Intensity), , drop = FALSE], 200)
    
    bands_df <- summarize_bands(spectrum_df$PPM, spectrum_df$Intensity)
    output_data$results$summary_table <- bands_df

    summary_text_ref <- if (referencing_successful) {
      sprintf("Spectrum calibrated using %s.\nExpected @ %.2f ppm; observed @ %.4f ppm → correction Δ = %.4f ppm (applied).\nDetected bulk solvent: %s.",
              referencing_info_out$calibration_standard, referencing_info_out$expected_ppm,
              referencing_info_out$found_peak_at_ppm, referencing_info_out$ppm_correction_applied,
              referencing_info_out$detected_solvent)
    } else {
      "WARNING: Automatic referencing failed. The chemical shifts (ppm) shown below are likely inaccurate."
    }
    
    present <- subset(bands_df, count > 0)
    bullets <- if (nrow(present)) {
      paste0("• ", present$Band, " — n=", present$count,
             ifelse(nzchar(present$examples_ppm), paste0(" (e.g., ", present$examples_ppm, ")"), ""),
             collapse = "\n")
    } else "• No clear band enrichments found."
    
    output_data$results$summary_text <- sprintf(
      "%s\n\nFunctional-group hints (by 1H band counts):\n%s",
      summary_text_ref, bullets
    )

    output_data$status <- "success"

  }, error = function(e) {
    log_message(paste("FATAL ERROR:", e$message))
    output_data$status <- "error"
    output_data$error <- e$message
  }, finally = {
    if (!is.null(temp_unzip_dir) && dir.exists(temp_unzip_dir)) {
      log_message(paste("Cleaning up:", temp_unzip_dir))
      unlink(temp_unzip_dir, recursive = TRUE, force = TRUE)
    }
  })

  return(output_data)
}