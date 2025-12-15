## gnps library to sqlite format

library(Spectra)
library(MsBackendMgf)
library(MsBackendSql)

#' Import and re-format GNPS MatchMS Cleaned MGF
mms <- c(
    rtime = "RTINSECONDS",
    acquisitionNum = "SCANS",
    precursorMz = "PRECURSOR_MZ",
    precursorIntensity = "PEPMASSINT",
    precursorCharge = "CHARGE",
    spectrum_id = "TITLE",
    compound_name = "COMPOUND_NAME",
    adduct = "ADDUCT",
    smiles = "SMILES",
    inchi = "INCHI",
    inchikey = "INCHIKEY",
    ms_mass_analyzer = "MS_MASS_ANALYZER",
    ionization = "MS_IONISATION",
    exactmass = "PARENT_MASS",
    collisionEnergy = "COLLISION_ENERGY",
    ms_manufacturer = "MS_MANUFACTURER",
    ms_dissociation_method = "MS_DISSOCIATION_METHOD",
    salt_ions = "SALT_IONS",
    formula = "FORMULA"
)

gnps <- Spectra("HUMAN_Ring_Trial/cleaned_spectra.mgf",
                source = MsBackendMgf(), mapping = mms)
#' Reformat variables
gnps$exactmass <- as.numeric(gnps$exactmass)
pol <- rep(-1L, length(gnps))
pol[gnps$IONMODE == "positive"] <- 1L
pol[gnps$IONMODE == "negative"] <- 0L
gnps$polarity <- pol

#' Drop redundant or not needed variables
keep <- c("msLevel", "acquisitionNum", "dataOrigin",
          "centroided", "polarity", "precursorMz", "precursorCharge",
          "collisionEnergy", "spectrum_id", "compound_name", "adduct",
          "smiles", "ms_mass_analyzer", "ionization", "inchikey",
          "inchi", "exactmass", "ms_manufacturer", "ms_dissociation_method",
          "salt_ions", "formula", "mz", "intensity")
dta <- spectraData(gnps@backend, keep)
dta$dataOrigin <- basename(dta$dataOrigin)

library(MsBackendSql)
library(RSQLite)
dbname <- "MsBackendSql.GNPS.matchms.cleaned.v1.sqlite"
con <- dbConnect(SQLite(), dbname)

be <- backendInitialize(MsBackendSql(), dbcon = con, data = dta)

s <- Spectra(be)

dbDisconnect(con)

