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
    spectra(mse)$mixture  <- sampleData(mse)[match(spectra(mse)$dataOrigin, sampleData(mse)$spectraOrigin), "mixture"]
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
detect_signal <- function(study_group, annotated = TRUE, lab = character(), bpparam)  {
    ## lab shouldf be character of length 1.
    a <- load_data(lab = lab, study_group = study_group, return = "mse")
    dr <- file.path("..", lab)
    if (annotated) {
        res <- read.csv(file.path(dr, "results",
                                  study_group, "ring_trial_library_HE.csv"))
        cpks <- chromPeaks(a)[res$chrom_peak_id, c("rtmin", "rtmax", "mzmin", "mzmax", "sample")]

    }
    else cpks <- chromPeaks(a)[, c("rtmin", "rtmax", "mzmin", "mzmax", "sample")]

    spectra(a) <- setBackend(spectra(a), MsBackendMemory())
    cpk_split <- split(as.data.frame(cpks), cpks[, "sample"])
    cpk_split <- lapply(cpk_split, function(df) as.matrix(df[, c("rtmin", "rtmax", "mzmin", "mzmax")]))
    bg <- bpmapply(FUN = function(s, pks) {
        s <- Spectra::filterPeaksRanges(s, mz = pks[, c("mzmin", "mzmax")],
                                        rtime = pks[, c("rtmin", "rtmax")],
                                        keep = TRUE)
        Spectra::applyProcessing(s)
    }, split(spectra(a), spectraSampleIndex(a)), cpk_split, BPPARAM = bpparam)
    res <- concatenateSpectra(bg)
    res <- filterEmptySpectra(res)
    return(res)
}



