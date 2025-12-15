## check for update - run this at least one at the beginning of the week
#BiocManager::install()

## Package load
library(readxl)
library(S4Vectors)
library(MsExperiment)
library(xcms) ## will also be devel for incorporating a few fixes
library(Spectra)
library(Biobase)
library(pheatmap)
library(alabaster.base)
library(MsIO)
library(RColorBrewer)
library(MetaboCoreUtils) ## devel version at the moment
library(MetaboAnnotation)
library(MsCoreUtils)
library(RSQLite)
library(MsBackendSql)
library(MsFeatures)
library(pander)

meta <- read_xlsx(
  file.path("1_preprocessing", "standards.xlsx"),
  col_names = TRUE,
  skip = 1,
  .name_repair = "minimal"
) |>
  as.data.frame(check.names = FALSE)

meta$M <- as.numeric(meta$M)

#' #' logp standard - not running this anymore
#' library(rcdk)
#' parsed <- parse.smiles(meta$SMILES)
#' logP <- vapply(parsed, get.xlogp, FUN.VALUE = numeric(1))
#' meta$theo_logP <- logP

## function for preprocessing

#' Helper function to extract a retention time matrix from the `matchedData`
#' `data.frame`. The function simply iterates over the sample and extracts the
#' retention time for each chrom peak matching a NAPS. If multiple chrom peaks
#' match a NAPS, the one with the highest intensity (`"maxo"`) is selected.
#' The function throws an error if the extracted retention times in one sample
#' are not ordered increasingly.
#'
#' @param x `data.frame` with at least columns `"rt"`, `"maxo"`, `"sample"`
#'     and `"target_Name"`
#'
#' @param naps_info `data.frame` with NAPS information. Needs to have a column
#'     `"Name"` with the name of the NAPS and with the data odered increasingly
#'     by retention time.
#'
#' @return `matrix` with retention times of the NAPS in the individual samples.
#'     Columns are samples and rows NAPS.
naps_rt_matrix <- function(x, naps_info) {
  smpls <- unique(x[, "sample"])
  res <- matrix(NA_real_, ncol = length(smpls), nrow = nrow(naps_info))
  rownames(res) <- naps_info$Name
  for (i in smpls) {
    z <- x[x$sample == i, ]
    rti <- split(z, factor(z$target_Name, naps_info$Name))
    rti <- do.call(
      rbind,
      lapply(rti, function(z) z[which.max(z$maxo), , drop = FALSE])
    )
    if (is.unsorted(rti$rt)) {
      warning("Retention times are not increasingly ordered for sample ", i)
      rti$rt <- force_sorted(rti$rt) # this can be improved.
    }
    res[match(rownames(rti), rownames(res)), i] <- rti$rt
  }
  res
}

## set up the sample data for future RTI
get_closest_index <- function(
  x,
  idx,
  method = c("next", "previous", "closest")
) {
  method <- match.arg(method)
  switch(
    method,
    `next` = {
      nxt <- idx > x
      if (any(nxt)) {
        idx[nxt][1]
      } else {
        idx[!nxt][sum(!nxt)]
      }
    },
    `previous` = {
      prv <- idx < x
      if (any(prv)) {
        idx[prv][sum(prv)]
      } else {
        idx[!prv][1]
      }
    },
    closest = {
      dst <- abs(idx - x)
      idx[which.min(dst)]
    }
  )
}


## function for library building

#' @title Peak matching
#' @param mse An xcms object containing the chromatographic data.
#' @param meta A data frame containing the metadata for the samples.
#' @param adducts A character vector of adducts to match against. Default is
#'  c("[M+H]+", "[M+Na]+", "[M+H-H2O]+").
#' @param ppm A numeric value representing the parts per million tolerance for
#' matching. Default is 10.
#' @param tol A numeric value representing the tolerance for matching. Default
#' is 0.05.
automate_matching <- function(meta, mse) {
  all_matches <- list()
  mixtures <- unique(meta$Mixture)

  for (x in mixtures) {
    # Subset for the current mixture
    x_meta <- split(meta, meta$Mixture)[[x]]
    x_mse <- mse[sampleData(mse)$mixture == x]

    cpks <- as.data.frame(chromPeaks(x_mse))
    cpks$chrom_peak_id <- rownames(cpks)

    # Positive mode
    adds_pos <- c("[M+H]+", "[M+Na]+")
    match_pos <- matchValues(
      cpks[cpks$sample == 1, ],
      x_meta,
      Mass2MzParam(adds_pos, tolerance = 0.05, ppm = 10),
      mzColname = "mz",
      massColname = "M"
    )
    match_pos <- match_pos[whichQuery(match_pos)]

    # Negative mode
    adds_neg <- c("[M-H]-", "[M+CHO2]-")
    match_neg <- matchValues(
      cpks[cpks$sample == 2, ],
      x_meta,
      Mass2MzParam(adds_neg, tolerance = 0.05, ppm = 10),
      mzColname = "mz",
      massColname = "M"
    )
    match_neg <- match_neg[whichQuery(match_neg)]

    # Extract results
    match_res_pos <- matchedData(
      match_pos,
      c(
        "chrom_peak_id",
        "mz",
        "rt",
        "mzmin",
        "mzmax",
        "rtmin",
        "rtmax",
        "into",
        "target_ChEBI name",
        "target_ChEBI",
        "target_InChIKey",
        "target_formula",
        "target_M",
        "adduct",
        "score",
        "ppm_error"
      )
    ) |>
      as.data.frame()
    match_res_pos$polarity <- "pos"

    match_res_neg <- matchedData(
      match_neg,
      c(
        "chrom_peak_id",
        "mz",
        "rt",
        "mzmin",
        "mzmax",
        "rtmin",
        "rtmax",
        "into",
        "target_ChEBI name",
        "target_ChEBI",
        "target_InChIKey",
        "target_formula",
        "target_M",
        "adduct",
        "score",
        "ppm_error"
      )
    ) |>
      as.data.frame()
    match_res_neg$polarity <- "neg"

    # Combine pos/neg
    match_res <- rbind(match_res_pos, match_res_neg)
    match_res$mixture <- x

    # Add adduct m/z
    match_res$adduct_mz <- mapply(
      match_res$target_M,
      match_res$adduct,
      FUN = mass2mz
    )

    # Add adduct formula
    af <- mapply(
      match_res$target_formula,
      match_res$adduct,
      FUN = adductFormula
    )
    match_res$adduct_formula <- gsub("^\\[|\\].*$", "", af)

    # Store results
    all_matches[[x]] <- match_res
  }

  # Combine all mixtures
  final_results <- do.call(rbind, all_matches)
  return(final_results)
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
    if (z$polarity == 1) {
      mzs <- sort(mass2mz(z$exactmass, adductNames("positive"))[1, ])
    } else {
      mzs <- sort(mass2mz(z$exactmass, adductNames("negative"))[1, ])
    }
    z <- filterMzValues(z, mzs, ppm = 10)
    idx <- MsCoreUtils::closest(mz(z)[[1L]], mzs, duplicates = "closest")
    z$peak_adduct_name <- list(names(mzs)[idx])
    z <- applyProcessing(z)
    z
  })
  adducts_combined <- concatenateSpectra(adducts_detected)

  ## In addition we add the number of adducts that have an intensity >= 0.5 the
  ## intensity of the chrom peak's intensity
  adduct_05 <- spectrapply(adducts_combined, function(z) {
    adct <- match_df[z$chrom_peak_id, "adduct"]
    idx <- which(z$peak_adduct_name[[1L]] == adct)
    sum(intensity(z)[[1L]] >= 0.5 * intensity(z)[[1L]][idx])
  }) |>
    unlist()

  match_df$ms1_adduct_count <- lengths(adducts_combined)
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
    isotopes,
    check_chemform(isotopes, match_df$adduct_formula)$new_formula,
    threshold = 0.001,
    charge = adductCharge(match_df$adduct),
    rel_to = 2
  )

  theoretical_spectra <- isopattern_to_spectra(ip)
  match_df$isopeak_count <- lengths(iso_spectra)
  match_df$isopeak_sim <- diag(compareSpectra(
    iso_spectra,
    theoretical_spectra,
    ppm = 20
  ))
  match_df

  list(
    data = match_df,
    spectra = iso_spectra,
    theoretical = theoretical_spectra
  )
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

#' @title MS2 matching
#' @description This function matches MS2 spectra to a reference library
#'
match_ms2_to_library <- function(
  ms2_spectra,
  ref_spectra,
  meta_df,
  min_score = 0.6
) {
  prm <- CompareSpectraParam(
    ppm = 40,
    tolerance = 0.05,
    requirePrecursor = TRUE,
    THRESHFUN = function(x) which(x >= min_score)
  )
  matches <- matchSpectra(ms2_spectra, ref_spectra, param = prm)
  mtch <- matches[whichQuery(matches)]

  df <- matchedData(mtch)[, c(
    "chrom_peak_id",
    "target_compound_name",
    "target_inchikey",
    "score",
    ".original_query_index",
    "target_adduct"
  )]
  df$target_index <- targetIndex(mtch)
  df$query_index <- queryIndex(mtch)
  df <- df |>
    as.data.frame() |>
    split(f = paste(df$chrom_peak_id, df$target_inchikey)) |>
    lapply(function(z) {
      z[which.max(z$score), ]
    }) |>
    do.call(what = rbind.data.frame)

  ## Count number of MS2 spectra per chrom peak with match
  match_summary <- df[!duplicated(df$.original_query_index), ]
  match_table <- table(match_summary$chrom_peak_id)

  ## TRUE positive: inchikey match with known standard
  match_summary$short_inchikey <- sub("-.*", "", match_summary$target_inchikey)
  meta_keys <- sub("-.*", "", meta_df$InChIKey)
  true_matches <- match_summary[match_summary$short_inchikey %in% meta_keys, ]

  list(match_obj = mtch, true_matches = true_matches, match_table = match_table)
}

#' Compute and update RTI values in match_res for each polarity
#'
#' @param mixture ID of the mixture to process (e.g., "1.1")
#' @param match_res data.frame from above
#' @param mse XcmsExperiemnt object ## that's the thing that is a bit dumb..
#' @param naps_res data.frame of NAPS results, from preprocessing
#' @param polarities Character vector of polarities to process (default: c("pos", "neg"))
#'
#' @return Updated match_res with new "RTI" values per polarity
update_all_with_rti <- function(
  match_res,
  mse,
  naps_res,
  polarities = c("pos", "neg")
) {
  # Internal RTI calculator
  compute_rti <- function(polarity, mixture, match_res, mse, naps_res) {
    mix_rt <- match_res[
      match_res$mixture == mixture & match_res$polarity == polarity,
      "rt"
    ]
    idx <- sampleData(mse)$mixture == mixture &
      sampleData(mse)$polarity == polarity

    naps_before <- sampleData(mse)$NAPS_before[idx]
    naps_after <- sampleData(mse)$NAPS_after[idx]

    testnaps <- naps_res[, c(naps_before, naps_after), drop = FALSE]
    av_naps <- rowMeans(testnaps, na.rm = TRUE)

    indexRtime(
      x = mix_rt,
      y = data.frame(rtime = av_naps, rindex = naps_res$RTI)
    )
  }

  # Loop through unique mixtures
  for (mix in unique(match_res$mixture)) {
    for (pol in polarities) {
      rti <- compute_rti(pol, mix, match_res, mse, naps_res)
      rows_to_update <- match_res$mixture == mix & match_res$polarity == pol
      match_res$RTI[rows_to_update] <- rti
    }
  }

  return(match_res)
}


#' @title Automated MS2 matching
#' @description This function automates the matching of MS2 spectra to
#' chromatographic peaks based on the provided match results.
#' @param mse2 An `XcmsExperiment` object containing the MS2 spectra.
#' @param match_res A data frame containing the match results from the
#' `automate_matching` function.
#' @param toleranceMz A numeric value representing the m/z tolerance for
#' matching. Default is 0.1.
#' @param toleranceRt A numeric value representing the retention time
#' tolerance for matching. Default is 5 seconds.
#' @return A list of matched MS2 spectra for each mixture.
automate_matching_ms2 <- function(
  mse2,
  match_res,
  toleranceMz = 0.1,
  toleranceRt = 5,
  waters_data = FALSE
) {
  all_matches <- list()
  mixtures <- unique(match_res$mixture)
  ## fix precursorMz

  for (x in mixtures) {
    #positive mode
    target <- match_res[
      match_res$mixture == x & match_res$polarity == "pos",
      c("rtmin", "rtmax", "mzmin", "mzmax")
    ]

    target$rtmin <- target$rtmin - toleranceRt
    target$rtmax <- target$rtmax + toleranceRt
    target$mzmin <- target$mzmin - toleranceMz
    target$mzmax <- target$mzmax + toleranceMz

    x_mse2 <- mse2[
      sampleData(mse2)$mixture == x &
        sampleData(mse2)$polarity == "pos"
    ]

    sp <- spectra(x_mse2)
    if (waters_data) {
      sp$precursorMz <- estimatePrecursorMz(sp)
    }
    sp <- filterMsLevel(sp, 2)

    sp_filt <- apply(
      as.matrix(target),
      MARGIN = 1,
      FUN = filterRanges,
      object = sp,
      spectraVariables = c("rtime", "precursorMz")
    )
    l <- lengths(sp_filt)
    sp_filt <- concatenateSpectra(sp_filt)
    sp_filt$chrom_peak_id <- rep(names(l), l)

    ## negative
    target <- match_res[
      match_res$mixture == x & match_res$polarity == "neg",
      c("rtmin", "rtmax", "mzmin", "mzmax")
    ]
    target$rtmin <- target$rtmin - toleranceRt
    target$rtmax <- target$rtmax + toleranceRt
    target$mzmin <- target$mzmin - toleranceMz
    target$mzmax <- target$mzmax + toleranceMz
    x_mse2 <- mse2[
      sampleData(mse2)$mixture == x &
        sampleData(mse2)$polarity == "neg"
    ]
    sp <- spectra(x_mse2)
    if (waters_data) {
      sp$precursorMz <- estimatePrecursorMz(sp)
    }
    sp <- filterMsLevel(sp, 2)
    sp_filt_neg <- apply(
      as.matrix(target),
      MARGIN = 1,
      FUN = filterRanges,
      object = sp,
      spectraVariables = c("rtime", "precursorMz")
    )
    l <- lengths(sp_filt_neg)
    sp_filt_neg <- concatenateSpectra(sp_filt_neg)
    sp_filt_neg$chrom_peak_id <- rep(names(l), l)
    sp_filt <- c(sp_filt, sp_filt_neg)

    ## add to list

    all_matches[[x]] <- sp_filt
  }
  all_sp <- concatenateSpectra(all_matches)
  return(all_sp)
}
