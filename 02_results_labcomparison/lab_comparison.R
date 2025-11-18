# load required libraries
library(tidyverse)
library(readxl)
library(ggvenn)
library(matrixStats)

# read data from laboratories
afekta <- read_xlsx("01_results_afekta/results/HE/peak_evidence_rt_grouped_manual.xlsx")
cembio <- read_xlsx("01_results_cembio/results/HE/peak_evidence_rt_grouped_manual.xlsx")
hmgu <- read_xlsx("01_results_hmgu/results/HE/peak_evidence_rt_grouped_manual.xlsx")
icl <- read_xlsx("01_results_hmgu/results/HE/peak_evidence_rt_grouped_manual.xlsx")

# filter only true values
afekta_filter <- afekta %>% 
  filter(`T/F` == TRUE)

cembio_filter <- cembio %>% 
  filter(`T/F` == TRUE)

hmgu_filter <- hmgu %>% 
  filter(`T/F` == TRUE)

icl_filter <- icl %>% 
  filter(`T/F` == TRUE)

# calculate average RI
afekta_filter$ri <- unlist(lapply(afekta_filter$RTI, function(x) {
  mean(as.numeric(unlist(str_split(x, "\\|"))))
}))

cembio_filter$ri <- unlist(lapply(cembio_filter$RTI, function(x) {
  mean(as.numeric(unlist(str_split(x, "\\|"))))
}))

hmgu_filter$ri <- unlist(lapply(hmgu_filter$RTI, function(x) {
  mean(as.numeric(unlist(str_split(x, "\\|"))))
}))

icl_filter$ri <- unlist(lapply(icl_filter$RTI, function(x) {
  mean(as.numeric(unlist(str_split(x, "\\|"))))
}))

# make venn diagram for comparison
venn <- list(afekta = afekta_filter$target_ChEBI.name,
             cembio = cembio_filter$target_ChEBI.name,
             hmgu = hmgu_filter$target_ChEBI.name,
             icl = icl_filter$target_ChEBI.name)

ggvenn(venn)

# make new data frame with RTs and RI
afekta_select <- afekta_filter %>% 
  select(target_ChEBI.name, target_ChEBI, target_InChIKey, target_formula,
         rtmed, ri) %>% 
  filter(ri > 300)

cembio_select <- cembio_filter %>% 
  select(target_ChEBI.name, target_ChEBI, target_InChIKey, target_formula,
         rtmed, ri) %>% 
  filter(ri > 300)

hmgu_select <- hmgu_filter %>% 
  select(target_ChEBI.name, target_ChEBI, target_InChIKey, target_formula,
         rtmed, ri) %>% 
  filter(ri > 300)

icl_select <- icl_filter %>% 
  select(target_ChEBI.name, target_ChEBI, target_InChIKey, target_formula,
         rtmed, ri) %>% 
  filter(ri > 300)

full_data <- full_join(afekta_select,
                       cembio_select,
                       suffix = c(".afekta", ".cembio"),
                       by = join_by(target_ChEBI.name, target_ChEBI, target_InChIKey, target_formula))

full_data <- full_join(full_data,
                       hmgu_select,
                       by = join_by(target_ChEBI.name, target_ChEBI, target_InChIKey, target_formula)) %>% 
  rename("rtmed.hmgu" = "rtmed",
         "ri.hmgu" = "ri")

full_data <- full_join(full_data,
                       icl_select,
                       by = join_by(target_ChEBI.name, target_ChEBI, target_InChIKey, target_formula)) %>% 
  rename("rtmed.icl" = "rtmed",
         "ri.icl" = "ri")


# find outlier metabolites
afekata_outlier <- full_data %>% 
  mutate(ri_mean = rowMeans(select(., c("ri.cembio", "ri.hmgu", "ri.icl"))),
         ri_sd = rowSds(as.matrix(select(., c("ri.cembio", "ri.hmgu", "ri.icl"))))) %>% 
  filter(abs((ri.afekta - ri_mean) / ri_mean) > 0.10) 

cembio_outlier <- full_data %>% 
  mutate(ri_mean = rowMeans(select(., c("ri.afekta", "ri.hmgu", "ri.icl")))) %>% 
  filter(abs((ri.cembio - ri_mean) / ri_mean) > 0.10) 

hmgu_outlier <- full_data %>% 
  mutate(ri_mean = rowMeans(select(., c("ri.afekta", "ri.cembio", "ri.icl")))) %>% 
  filter(abs((ri.hmgu - ri_mean) / ri_mean) > 0.10)  

icl_outlier <- full_data %>% 
  mutate(ri_mean = rowMeans(select(., c("ri.afekta", "ri.cembio", "ri.hmgu")))) %>% 
  filter(abs((ri.icl - ri_mean) / ri_mean) > 0.10) 
