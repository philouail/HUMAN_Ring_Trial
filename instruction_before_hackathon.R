## Install Latest R version
## Install latest Rstudio version
## install Github desktop
## Make a github account

# Install BiocManager if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# List of all packages to install via BiocManager
bioc_packages <- c(
    "readxl",
    "S4Vectors",
    "MsExperiment",
    "Spectra",
    "Biobase",
    "pheatmap",
    "alabaster.base",
    "MsIO",
    "RColorBrewer",
    "MetaboAnnotation",
    "MsCoreUtils",
    "RSQLite",
    "MsBackendSql",
    "MsFeatures"
)

# Install all packages with dependencies and force reinstallation
BiocManager::install(bioc_packages, dependencies = TRUE, force = TRUE)

BiocManager::install("MetaboCoreUtils") ## can the devel be loaded from bioconductor now ? ask johannes.
BiocManager::install("xcms", ref = "devel")

## Preparing metadata file
## see the example seq sheet in the shared_data folder
## this was the one created by Michael and i added one column of info.
## if you r injection does nto look exactly like that
## please update is with the same format but matching exactly your
## injection order, filenames,...)


