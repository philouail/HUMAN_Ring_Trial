# Compute median BPC/TIC ratio per RT bin
compute_bin_ratio <- function(rt, int_tic, int_bpc, n_bins = 30) {
    bins <- cut(rt, breaks = n_bins)
    tapply(seq_along(rt), bins, function(idx) {
        median(int_bpc[idx] / int_tic[idx], na.rm = TRUE)
    })
}

# Compute number of data points per RT bin
compute_bin_count <- function(rt, n_bins = 30) {
    bins <- cut(rt, breaks = n_bins)
    tapply(seq_along(rt), bins, length)
}


## load data
#' @return spectra object

load_data <- function(lab, study_group, return = "sp"){
    ## Load MS1 level data.
    dr <- file.path("..", lab)
    mse <- readMsObject(XcmsExperiment(),
                        AlabasterParam(path = file.path(dr, "results", study_group, "mse")),
                        spectraPath = file.path(dr, "HE_mzml")
    )

    sampleData(mse)$mixture <- sub(".*_", "", sampleData(mse)$Sample.Name)
    sampleData(mse)$mixture <- gsub("\\.", "_", sampleData(mse)$mixture)

    spectra(mse)$lab <- lab
    spectra(mse)$mixture  <- sampleData(mse)[match(spectra(mse)$dataOrigin,
                                             sampleData(mse)$spectraOrigin),
                                             "mixture"]
    sp <- spectra(mse)
    if (length(unique(sp$msLevel)) > 1) sp <- filterMsLevel(sp, 1L)
    if (return == "sp")
        return(sp)
    else if (length(unique(spectra(mse)$msLevel)) > 1)
            mse <- filterMsLevel(mse, 1L)
        else return(mse)
}

##create_detect object
#' @return spectra object
detect_signal <- function(study_group, annotated = TRUE, lab = character(),
                          bpparam)  {
    ## lab should be character of length 1.
    a <- load_data(lab = lab, study_group = study_group, return = "mse")
    dr <- file.path("..", lab)
    if (annotated) {
        res <- read.csv(file.path(dr, "results",
                                  study_group, "ring_trial_library_HE.csv"))
        cpks <- chromPeaks(a)[res$X,
                              c("rtmin", "rtmax", "mzmin", "mzmax", "sample")]

    }
    else cpks <- chromPeaks(a)[, c("rtmin", "rtmax", "mzmin", "mzmax", "sample")]
    spectra(a) <- setBackend(spectra(a), MsBackendMemory())
    cpk_split <- split(as.data.frame(cpks), cpks[, "sample"])
    cpk_split <- lapply(cpk_split, function(df) as.matrix(df[,
                                                             c("rtmin",
                                                               "rtmax",
                                                               "mzmin",
                                                               "mzmax")]))
    bg <- bpmapply(FUN = function(s, pks) {
        s <- Spectra::filterPeaksRanges(s, mz = pks[, c("mzmin", "mzmax")],
                                        rtime = pks[, c("rtmin", "rtmax")],
                                        keep = TRUE)
        Spectra::applyProcessing(s)
    }, split(spectra(a), spectraSampleIndex(a)), cpk_split, BPPARAM = bpparam)
    bg_full <- concatenateSpectra(bg)
    bg_full <- filterEmptySpectra(bg_full)
    return(bg_full)
}

## For Savitzky-Golay Smoothing
setClass("SavitzkyGolayParam",
         slots = c(window = "integer",
                   polynomial = "integer"),
         prototype = list(window = 11L,
                          polynomial = 2L),
         validity = function(object) {
             if (object@window %% 2 == 0 || object@window < 3) {
                 return("Smoothing 'window' must be an odd integer >= 3.")
             }
             if (object@polynomial >= object@window) {
                 return("'polynomial' order must be less than the 'window' size.")
             }
             return(TRUE)
         })

#' Constructor for SavitzkyGolayParam
SavitzkyGolayParam <- function(window = 11, polynomial = 2) {
    new("SavitzkyGolayParam", window = as.integer(window), polynomial = as.integer(polynomial))
}


## For Isolating the Main Peak
setClass("IsolatePeakParam",
         slots = c(frac = "numeric", tail_point = "numeric"),
         prototype = list(frac = 0.05, tail_point = 1),
         validity = function(object) {
             if (object@frac <= 0 || object@frac >= 1) {
                 return("'frac' must be between 0 and 1.")
             }
             return(TRUE)
         })
#' Constructor for IsolatePeakParam
IsolatePeakParam <- function(frac = 0.05, tail_point = 1) {
    new("IsolatePeakParam", frac = frac, tail_point = tail_point)
}

## -----------------------------------------------------------------------------
## Generic and Methods for processEICs
## -----------------------------------------------------------------------------

# Define the generic function
setGeneric("processEICs", function(object, param) {
    standardGeneric("processEICs")
})

# Method for Smoothing
setMethod("processEICs", signature(object = "Chromatograms", param = "SavitzkyGolayParam"),
          function(object, param) {
              df_list <- peaksData(object)

              res <- lapply(seq_along(df_list), function(i) {
                  df <- df_list[[i]]
                  if (nrow(df) >= param@window) {
                      smoothed <- signal::sgolayfilt(df$intensity, p = param@polynomial, n = param@window)
                      df$intensity <- pmax(0, smoothed)
                  } else {
                      warning(sprintf("EIC %d: Skipped smoothing (n=%d < window=%d)", i, nrow(df), param@window), call. = FALSE)
                  }
                  return(df)
              })

              peaksData(object) <- res
              object
          })

# Method for Peak Isolation
setMethod("processEICs", signature(object = "Chromatograms", param = "IsolatePeakParam"),
          function(object, param) {
              df_list <- peaksData(object)

              res <- lapply(df_list, function(df) {
                  if (nrow(df) < 3) return(df)

                  peaks <- pracma::findpeaks(df$intensity, sortstr = TRUE)
                  if (is.null(peaks)) {
                      warning("No peaks found in chromatograms; returning empty data frame.", call. = FALSE)
                      return(df[F,]) # Return empty if no peaks
                  }

                  main_center <- peaks[1, 2]
                  apex <- df$intensity[main_center]
                  threshold <- apex * param@frac
                  n <- nrow(df)

                  left <- main_center
                  while (left > 1 && df$intensity[left - 1] < df$intensity[left] && df$intensity[left - 1] > threshold) {
                      left <- left - 1
                  }
                  right <- main_center
                  while (right < n && df$intensity[right + 1] < df$intensity[right] && df$intensity[right + 1] > threshold) {
                      right <- right + 1
                  }
                  final_left <- max(1, left - param@tail_point)
                  final_right <- min(n, right + param@tail_point)
                  # Return the subset of the original dataframe
                  return(df[final_left:final_right, ])
              })
              peaksData(object) <- res
              object
          })

#' Flag Potentially Problematic EICs
#'
#' This function evaluates a list of EICs based on several quality metrics
#' and returns a data frame of those that fail the checks.
#'
#' @param df_list A list of data frames, where each has 'rtime' and 'intensity'.
#' @param min_points The minimum number of data points an EIC must have.
#' @param min_snr The minimum signal-to-noise ratio required. Noise is calculated
#'   as the Median Absolute Deviation (MAD), a robust metric.
#' @param max_peaks The maximum number of local maxima (peaks) an EIC can have.
#'   More peaks suggest a noisy, spiky signal rather than a clean one.
#' @param min_cv The minimum coefficient of variation (sd / mean). Helps flag
#'   unnaturally flat signals.
#'
#' @return A data frame with 'index' and 'reason' columns for each flagged EIC.
#'   Returns an empty data frame if all EICs pass.

flagBadEICs <- function(df_list,
                        min_points = 10,
                        min_snr = 3,
                        max_peaks = 8,
                        min_cv = 0.1) {

    # Ensure the 'pracma' package is available for findpeaks()
    if (!requireNamespace("pracma", quietly = TRUE)) {
        stop("Package 'pracma' is required. Please install it with install.packages('pracma')")
    }

    flagged_results <- list()

    # Loop through each EIC with its index
    for (i in seq_along(df_list)) {
        df <- df_list[[i]]
        reasons <- c() # Store failure reasons for the current EIC

        # --- Perform Checks ---

        # 1. Check for sufficient data points
        if (nrow(df) < min_points) {
            reasons <- c(reasons, "too_few_points")
        }

        # Skip further checks if the data is empty or too small
        if (nrow(df) > 2 && sum(df$intensity, na.rm = TRUE) > 0) {

            # 2. Check Signal-to-Noise Ratio (S/N)
            signal <- max(df$intensity, na.rm = TRUE)
            noise <- mad(df$intensity, na.rm = TRUE) # MAD is robust to outliers
            snr <- if (noise > 0) signal / noise else Inf
            if (snr < min_snr) {
                reasons <- c(reasons, "low_snr")
            }

            # 3. Check for excessive noise/spikes by counting local peaks
            num_peaks <- pracma::findpeaks(df$intensity)
            n_peaks <- if (is.null(num_peaks)) 0 else nrow(num_peaks)
            if (n_peaks > max_peaks) {
                reasons <- c(reasons, "too_spiky")
            }

            # 4. Check for unnaturally flat signals using Coefficient of Variation (CV)
            mean_int <- mean(df$intensity, na.rm = TRUE)
            sd_int <- sd(df$intensity, na.rm = TRUE)
            cv <- if (mean_int > 0) sd_int / mean_int else 0
            if (cv < min_cv) {
                reasons <- c(reasons, "too_flat")
            }
        } else if (nrow(df) > 0) {
            reasons <- c(reasons, "no_signal") # If not empty but sum is 0
        }


        # If any check failed, record the index and the reasons
        if (length(reasons) > 0) {
            flagged_results[[length(flagged_results) + 1]] <- data.frame(
                index = i,
                reason = paste(reasons, collapse = "; ")
            )
        }
    }

    # Combine list of data frames into a single one
    if (length(flagged_results) == 0) {
        return(data.frame(index = integer(0), reason = character(0)))
    } else {
        return(do.call(rbind, flagged_results))
    }
}

#' @title S4 Parameter Class for ALS Baseline Correction
#' @description Holds parameters for baseline correction using the ptw package.
#' @slot lambda numeric. The smoothness parameter (default is 1e5).
#' @slot p numeric. The asymmetry parameter (default is 0.01).
setClass("BaselineParam",
         slots = c(lambda = "numeric",
                   p = "numeric"),
         prototype = list(lambda = 1e5,
                          p = 0.01),
         validity = function(object) {
             if (object@lambda <= 0) {
                 return("'lambda' must be a positive number.")
             }
             if (object@p <= 0 || object@p >= 1) {
                 return("'p' must be between 0 and 1.")
             }
             return(TRUE)
         })

#' @title Constructor for BaselineParam
#' @param lambda Smoothness parameter.
#' @param p Asymmetry parameter.
#' @return A BaselineParam object.
BaselineParam <- function(lambda = 1e5, p = 0.01) {
    new("BaselineParam", lambda = lambda, p = p)
}


#' @title processEICs Method for Baseline Correction
#' @description Applies Asymmetric Least Squares baseline correction to each EIC.
setMethod("processEICs", signature(object = "Chromatograms", param = "BaselineParam"),
          function(object, param) {

              # Ensure the required package is installed
              if (!requireNamespace("ptw", quietly = TRUE)) {
                  stop("Package 'ptw' is required for baseline correction. Please install it.")
              }

              df_list <- peaksData(object)

              res <- lapply(df_list, function(df) {
                  # Proceed only if there's data to process
                  if (nrow(df) > 0 && sum(df$intensity, na.rm = TRUE) > 0) {

                      # Perform the baseline correction
                      corrected_intensity <- ptw::baseline.corr(df$intensity,
                                                                lambda = param@lambda,
                                                                p = param@p)

                      # Set any values below zero to zero
                      df$intensity <- pmax(0, corrected_intensity)
                  }
                  return(df)
              })

              peaksData(object) <- res
              object
          })


#' Calculate Comprehensive Metrics for a List of EICs
#'
#' This function processes a list of EIC data frames (with 'rtime' and
#' 'intensity' columns) and computes several key chromatographic peak metrics.
#'
#' @param df_list A list of EIC data frames.
#'
#' @return A data frame where each row corresponds to an EIC and each column
#'   is a calculated metric.
#' Calculate Comprehensive Metrics for a List of EICs
#'
#' This function processes a list of EIC data frames (with 'rtime' and
#' 'intensity' columns) and computes several key chromatographic peak metrics.
#'
#' @param df_list A list of EIC data frames.
#'
#' @return A data frame where each row corresponds to an EIC and each column
#'   is a calculated metric.
#' Calculate Comprehensive Metrics for a List of EICs
#'
#' This function processes a list of EIC data frames (with 'rtime' and
#' 'intensity' columns) and computes several key chromatographic peak metrics.
#'
#' @param df_list A list of EIC data frames.
#'
#' @return A data frame where each row corresponds to an EIC and each column
#'   is a calculated metric.
calculatePeakMetrics <- function(df_list) {

    # Check for the 'pracma' package, needed for AUC calculation
    if (!requireNamespace("pracma", quietly = TRUE)) {
        stop("Package 'pracma' is required. Please install it with install.packages('pracma')")
    }

    # Use lapply to iterate over each EIC and calculate its metrics
    metrics_list <- lapply(seq_along(df_list), function(i) {
        df <- df_list[[i]]

        # Handle empty or very small EICs by returning a row of NAs
        if (is.null(df) || nrow(df) < 3) {
            return(data.frame(eic_index = i,
                              rt_apex = NA_real_,
                              apex_intensity = NA_real_,
                              num_points = nrow(df),
                              fwhm = NA_real_,
                              tailing_factor = NA_real_,
                              gaussian_similarity = NA_real_, # <-- NEW
                              rt_width = NA_real_,
                              auc = NA_real_,
                              baseline = NA_real_,
                              bs_apex_ratio = NA_real_,
                              entropy = NA_real_))
        }

        # 1. Apex Intensity and Retention Time
        apex_intensity <- max(df$intensity, na.rm = TRUE)
        apex_idx <- which.max(df$intensity)
        rt_apex <- df$rtime[apex_idx]

        # 2. Number of Data Points
        num_points <- nrow(df)

        # Define peak sides for FWHM and Tailing Factor
        rising_side <- df[1:apex_idx, ]
        falling_side <- df[apex_idx:num_points, ]

        # 3. Full Width at Half Maximum (FWHM)
        fwhm <- tryCatch({
            half_max <- apex_intensity / 2
            if (min(rising_side$intensity, na.rm = TRUE) <= half_max &&
                min(falling_side$intensity, na.rm = TRUE) <= half_max) {
                rt_left <- approx(rising_side$intensity, rising_side$rtime, xout = half_max)$y
                rt_right <- approx(falling_side$intensity, falling_side$rtime, xout = half_max)$y
                rt_right - rt_left
            } else {
                NA_real_
            }
        }, error = function(e) NA_real_)

        # 4. Tailing Factor (at 10% height)
        tailing_factor <- tryCatch({
            height_10pct <- apex_intensity * 0.10
            if (min(rising_side$intensity, na.rm = TRUE) <= height_10pct &&
                min(falling_side$intensity, na.rm = TRUE) <= height_10pct) {
                rt_left_10 <- approx(rising_side$intensity, rising_side$rtime, xout = height_10pct)$y
                rt_right_10 <- approx(falling_side$intensity, falling_side$rtime, xout = height_10pct)$y
                if (!is.na(rt_left_10) && !is.na(rt_right_10)) {
                    a <- rt_apex - rt_left_10
                    b <- rt_right_10 - rt_apex
                    if (a > 0) b / a else NA_real_
                } else {
                    NA_real_
                }
            } else {
                NA_real_
            }
        }, error = function(e) NA_real_)

        # 5. Gaussian Similarity (R-squared) <-- NEW
        gaussian_similarity <- tryCatch({
            # This metric depends on a valid FWHM
            if (is.na(fwhm) || fwhm == 0) {
                NA_real_
            } else {
                # Convert FWHM to sigma (std. deviation)
                sigma <- fwhm / (2 * sqrt(2 * log(2)))

                # Generate the ideal Gaussian curve
                y_ideal <- apex_intensity * exp(-((df$rtime - rt_apex)^2) / (2 * sigma^2))

                # Calculate R^2 (coefficient of determination)
                cor(df$intensity, y_ideal)^2
            }
        }, error = function(e) NA_real_)

        # 6. Retention Time Width (at the base)
        rt_width <- max(df$rtime) - min(df$rtime)

        # 7. Area Under the Curve (AUC)
        auc <- pracma::trapz(df$rtime, df$intensity)

        # 8. Value at Peak Ground (Baseline)
        baseline <- mean(c(df$intensity[1], df$intensity[num_points]), na.rm = TRUE)

        # 9. Baseline-to-Apex Ratio
        bs_apex_ratio <- ifelse(apex_intensity > 0, baseline / apex_intensity, NA_real_)

        # 10. Entropy (as a measure of peak shape complexity)
        total_intensity <- sum(df$intensity, na.rm = TRUE)
        entropy <- if (total_intensity > 0) {
            p <- df$intensity / total_intensity
            p_nonzero <- p[p > 0]
            -sum(p_nonzero * log2(p_nonzero))
        } else {
            NA_real_
        }

        # Combine all metrics into a single-row data frame
        data.frame(eic_index = i,
                   rt_apex = rt_apex,
                   apex_intensity = apex_intensity,
                   num_points = num_points,
                   fwhm = fwhm,
                   tailing_factor = tailing_factor,
                   gaussian_similarity = gaussian_similarity, # <-- ADDED
                   rt_width = rt_width,
                   auc = auc,
                   baseline = baseline,
                   bs_apex_ratio = bs_apex_ratio,
                   entropy = entropy)
    })

    # Combine the list of results into a single output data frame
    return(do.call(rbind, metrics_list))
}


#' Plot Signal Comparison Violin Plot
#'
#' Generates a grouped violin plot comparing the signal area (TIC or BPC)
#' from Full, Detected, and Annotated data across labs for a specific polarity.
#'
#' @param signal_to_plot Character string, either "TIC" or "BPC".
#' @param polarity_to_plot Numeric, either 1 (for positive) or 0 (for negative).
#'
#' @return A ggplot object.
plot_signal_comparison <- function(signal_to_plot = "TIC", polarity_to_plot = 1) {

    # --- A. Load Data Objects ---
    # Dynamically build file names based on the parameter
    file_suffix <- tolower(signal_to_plot)

    tryCatch({
        load(paste0("object/", file_suffix, "_full.RData"))
        load(paste0("object/", file_suffix, "_detect.RData"))
        load(paste0("object/", file_suffix, "_ann.RData"))
    }, error = function(e) {
        stop("Could not load data files. Make sure objects like 'tic_full.RData' or 'bpc_full.RData' exist in the 'object/' folder.")
    })

    # Use generic names for the loaded objects
    sig_full <- get(paste0(file_suffix, "_full"))
    sig_detect <- get(paste0(file_suffix, "_detect"))
    sig_ann <- get(paste0(file_suffix, "_ann"))

    # --- B. Filter for Polarity ---
    sig_full <- sig_full[sig_full$polarity == polarity_to_plot]
    sig_detect <- sig_detect[sig_detect$polarity == polarity_to_plot]
    sig_ann <- sig_ann[sig_ann$polarity == polarity_to_plot]

    pol_label <- ifelse(polarity_to_plot == 1, "Positive", "Negative")

    # --- C. Process and Combine Data ---

    # Process "Full" data
    s_full <- vapply(intensity(sig_full), sum, numeric(1), na.rm = TRUE)
    data_full <- data.frame(
        Lab = sig_full$lab,
        Mixture = sig_full$mixture,
        Area = s_full,
        Signal_Type = "Full"
    )

    # Process "Detected" data
    s_det <- vapply(intensity(sig_detect), sum, numeric(1), na.rm = TRUE)
    data_det <- data.frame(
        Lab = sig_detect$lab,
        Mixture = sig_detect$mixture,
        Area = s_det,
        Signal_Type = "Detected"
    )

    # Process "Annotated" data
    s_ann <- vapply(intensity(sig_ann), sum, numeric(1), na.rm = TRUE)
    data_ann <- data.frame(
        Lab = sig_ann$lab,
        Mixture = sig_ann$mixture,
        Area = s_ann,
        Signal_Type = "Annotated"
    )

    # Combine all three into one data frame
    all_data <- bind_rows(data_full, data_det, data_ann)

    # --- D. Prepare for Plotting ---

    # Convert Area to log2 scale for better visualization
    all_data <- all_data %>%
        mutate(log2_Area = log2(Area))

    # Set the order of the Signal_Type for the plot
    all_data$Signal_Type <- factor(all_data$Signal_Type,
                                   levels = c("Full", "Detected", "Annotated"))

    # --- E. Generate the Plot ---
    comparison_plot <- ggplot(all_data,
                              aes(x = Lab, y = log2_Area, fill = Signal_Type)) +

        # Create the violin plots, dodged side-by-side
        geom_violin(
            position = position_dodge(width = 0.9),
            alpha = 0.7,
            trim = TRUE,
            scale = "width" # Makes all violins have the same width
        ) +

        # Use a color-blind friendly palette
        scale_fill_brewer(palette = "Set1") +

        # Dynamic labels
        labs(
            title = paste(signal_to_plot,
                          "Area Comparison by Lab and Signal Type (", pol_label, "Polarity)"),
            subtitle = "Full (Noise+Signal) vs. Detected (Signal) vs. Annotated (Standards)",
            x = "Laboratory",
            y = paste("log2(Total", signal_to_plot, "Area)"),
            fill = "Signal Type"
        ) +

        theme_minimal() +
        theme(
            legend.position = "bottom",
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    # Return the plot object
    return(comparison_plot)
}

