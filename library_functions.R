## function for library building

#' @title Peak matching
#' @param xcms_obj An xcms object containing the chromatographic data.
#' @param meta A data frame containing the metadata for the samples.
#' @param adducts A character vector of adducts to match against. Default is
#'  c("[M+H]+", "[M+Na]+", "[M+H-H2O]+").
#' @param ppm A numeric value representing the parts per million tolerance for
#' matching. Default is 10.
#' @param tol A numeric value representing the tolerance for matching. Default
#' is 0.05.
match_chrom_peaks <- function(xcms_obj, meta, adducts = c("[M+H]+", "[M+Na]+",
                                                          "[M+H-H2O]+"),
                              ppm = 10, tol = 0.05) {
    cpks <- as.data.frame(chromPeaks(xcms_obj))
    cpks$chrom_peak_id <- rownames(cpks)
    #' match chrom peaks to possible ions
    matched <- matchValues(cpks, meta, Mass2MzParam(adducts, tolerance = tol, ppm = ppm), mzColname = "mz", massColname = "M")
    matched <- matched[whichQuery(matched)]
    #' extract and process the results
    res <- matchedData(matched, c("chrom_peak_id", "rtmin", "rtmax", "into",
                                  "target_ChEBI name", "target_ChEBI",
                                  "target_InChIKey", "target_formula",
                                  "target_M", "adduct", "score", "ppm_error")) |>
        as.data.frame()
    #' Add the m/z of the adduct
    res$adduct_mz <- mapply(res$target_M, res$adduct, FUN = mass2mz)
    #' Calculate the chemical formula of the adduct
    af <- mapply(res$target_formula, res$adduct, FUN = adductFormula)
    res$adduct_formula <- gsub("^\\[|\\].*$", "", af)
    res
}

#' @title MS1 Adduct analysis
#' @description This function extracts MS1 adduct information from the
#' chromatographic data and check for theoretical adduct that could be
#' present in the MS1 spectrum.
#' @param ms1_spectra An xcms object containing the MS1 spectra.
#' @param match_df A data frame containing the matched chromatographic data.
#'                 This will be the output of the previous step.
#' @return A list containing the updated match_df with additional columns for
#' the MS1 adduct information and the MS1 spectra with the adduct peaks
#' identified. THe later can be used for plotting later on.
extract_ms1_adduct_info <- function(ms1_spectra, match_df) {
    ms1_spectra$exactmass <- match_df$target_M
    ms1_spectra$adduct_mz <- match_df$adduct_mz
    ms1_spectra <- setBackend(backend = MsBackendMemory(), object = ms1_spectra)

    ## determine which other adduct peaks would be present in the MS1 spectrum.
    adducts_detected <- spectrapply(ms1_spectra, function(z) {
        mzs <- sort(mass2mz(z$exactmass, adductNames("positive"))[1, ])
        z <- filterMzValues(z, mzs, ppm = 10)
        idx <- MsCoreUtils::closest(mz(z)[[1L]], mzs)
        z$peak_adduct_name <- list(names(mzs)[idx])
        applyProcessing(z)
    })
    adducts_combined <- concatenateSpectra(adducts_detected)

    adduct_count <- lengths(adducts_combined)
    ## In addition we add the number of adducts that have an intensity >= 0.5 the
    ## intensity of the chrom peak's intensity
    adduct_05 <- spectrapply(adducts_combined, function(z) {
        adct <- match_df[z$chrom_peak_id, "adduct"]
        idx <- which(z$peak_adduct_name[[1L]] == adct)
        sum(intensity(z)[[1L]] >= 0.5 * intensity(z)[[1L]][idx])
    }) |> unlist()

    match_df$ms1_adduct_count <- adduct_count
    match_df$ms1_adduct_05_count <- adduct_05
    list(data = match_df, spectra = adducts_combined)
}


#' @title Isotope pattern matching
#' @description This function extracts the isotope pattern information from
#' the chromatographic data and check for theoretical isotope pattern that
#' could be present in the MS1 spectrum.
#' @param ms1_spectra An xcms object containing the MS1 spectra.
#' @param match_df A data frame containing the matched chromatographic data.
#'                This will be the output of the previous step.
#' @return A list containing the updated match_df with additional columns for
#' the isotope pattern information and the MS1 spectra with the isotope
#' peaks identified.
calculate_isotope_similarity <- function(ms1_spectra, match_df) {
    library(enviPat)
    data(isotopes)
    #' For each MS1 spectrum, extract the peaks matching potential isotope peaks.
    #' determine which other adduct peaks would be present in the MS1 spectrum.
    iso_spectra <- spectrapply(ms1_spectra, function(z) {
        idx <- isotopologues(peaksData(z)[[1L]], ppm = 20, seedMz = z$adduct_mz)
        if (length(idx) == 1L) {
            addProcessing(z, function(x, i = idx[[1L]], ...) x[i, , drop = FALSE]) |>
                applyProcessing()
        } else {
            filterMzValues(z, mz = z$adduct_mz, ppm = 20, tolerance = 0) |>
                applyProcessing()
        }
    }) |>
        concatenateSpectra() |>
        scalePeaks()
    #' Create theoretical isotope pattern for all chrom peaks/adducts
    ip <- isopattern(
        isotopes, check_chemform(isotopes, match_df$adduct_formula)$new_formula,
        threshold = 0.001, charge = adductCharge(match_df$adduct), rel_to = 2
    )

    theoretical_spectra <- isopattern_to_spectra(ip)
    match_df$isopeak_count <- lengths(iso_spectra)
    match_df$isopeak_sim <- diag(compareSpectra(iso_spectra, theoretical_spectra, ppm = 20))
    match_df

    list(data = match_df, spectra = iso_spectra)
}


## helper function for isotope pattern matching
adductCharge <- function(x) {
    MetaboCoreUtils:::.process_adduct_arg(x, "charge")
}
isopattern_to_spectra <- function(x) {
    df <- data.frame(msLevel = 1L, formula = names(x))
    df$mz <- lapply(x, function(z) z[, 1L])
    df$intensity <- lapply(x, function(z) z[, 2L])
    Spectra(df)
}


## Plotting function


## other ?
