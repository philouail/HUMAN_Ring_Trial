----- Install packages -----

## Install Latest R version 4.5:
## Install latest Rstudio version
## instruction can be found here: https://rstudio-pubs-static.s3.amazonaws.com/1215682_b2ed746543234e2c961ea4b67b4da290.html

# Install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.21")

# List of all packages to install via BiocManager
bioc_packages <- c(
    "readxl",
    "S4Vectors",
    "MsExperiment",
    "Spectra",
    "Biobase",
    "pheatmap",
    "alabaster.base",
    "RColorBrewer",
    "MetaboAnnotation",
    "MsCoreUtils",
    "RSQLite",
    "MsBackendSql",
    "MsFeatures",
    "MsBackendMgf",
    "remotes",
    "pander",
    "shinyjs" # that's for the blood mapping annotation.
)

# Install all packages with dependencies and force reinstallation
BiocManager::install(bioc_packages, dependencies = TRUE, force = TRUE)

BiocManager::install(c("rformassspectrometry/MetaboCoreUtils",
                       "sneumann/xcms", "rformassspectrometry/MsIO"))
## if the code above is bugging please tell me.

## Check the packages are okay - please check the output in the R console,
##  error message hide well. Do they one by one to be sure not to miss anything.
library(readxl)
library(pander)
library(S4Vectors)
library(MsExperiment)
library(xcms)
library(Spectra)
library(Biobase)
library(pheatmap)
library(alabaster.base)
library(MsIO)
library(RColorBrewer)
library(MetaboCoreUtils)
library(MetaboAnnotation)
library(MsCoreUtils)
library(RSQLite)
library(MsBackendSql)
library(MsFeatures)
library(MsBackendMgf)
## if R is telling you some other package is missing please install them as such:
## `BiocManager::install("alabaster.base")`, and replace alabaster.base with
## the package name it is asking for.

----- Metadata preparation -----

## Preparing metadata file:
## see the example seq sheets I sent you.
## this was the one created by Michael and i added 2 columns of info.
## if your injections does not look exactly like that
## please update is with the same format but matching EXACTLY your
## injection order, filenames,...) the ".mzml is not necessary"
##
## If you have duplicate/triplicate. please just choose ONE file per mixture.
## Maybe the one with the most chromatographic peaks / higher signals/ less noise..
##
## for people that have spearate MS1 and MS2 injections put both of them. in
## the right injection order.

----- Raw data preparation -----
## For ring trial:
## have your raw data centroided, mzml and NAPS and BLANKS, pos and neg.
##
## For blood mapping:
## have an .mgf file of the preprocessing output ( if you can filter to only
## have MS2 that is great, if not just tell me so i adjust that in the code
## later.)
## Bring you MS1 annotation results.


