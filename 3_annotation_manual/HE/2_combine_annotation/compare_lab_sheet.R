library(tidyverse)
library(readxl)
library(openxlsx)

# 1. Setup paths
lab_data_folder <- "post_manual_curation"

# 2. Read and Combine Data
all_labs_data <- list.files(
  path = lab_data_folder,
  pattern = "\\.xlsx$",
  full.names = TRUE
) %>%
  set_names(basename(.)) %>%
  map_dfr(~ read_excel(.x), .id = "lab_file") %>%
  mutate(
    lab_id = lab_file %>%
      str_remove("\\.xlsx$") %>% # Remove extension
      str_remove("^peak_evidence_rt_grouped_manual_fixed_") # Remove prefix
  )

## Export table of all annotations (raw combined)
write.xlsx(
  all_labs_data,
  file = "combine_annotation/all_lab_annotations.xlsx"
)

df <- all_labs_data


# Define the universe of labs (to determine who is missing)
all_participating_labs <- unique(df$lab_id)

# A. Calculate General Stats
general_stats <- df %>%
  mutate(
    # Clean RTI for calculation
    RTI_clean = map_dbl(
      str_split(RTI, "\\|"),
      ~ mean(as.numeric(.x), na.rm = TRUE)
    )
  ) %>%
  group_by(target_ChEBI.name) %>%
  summarise(
    # --- LAB COUNT ---
    n_labs = n_distinct(lab_id),

    # --- DETECTION REPORTING ---
    detected_labs = list(unique(lab_id)),

    # RT string (Now clean: "afekta:397.33 | icl:399.10")
    rt_per_lab = paste(paste0(lab_id, ":", round(rtmed, 2)), collapse = " | "),

    # --- RT REPORTING ---
    rt_min = min(rtmed, na.rm = TRUE),
    rt_max = max(rtmed, na.rm = TRUE),
    rt_mean = mean(rtmed, na.rm = TRUE),
    rt_sd = sd(rtmed, na.rm = TRUE),

    # --- RTI COMPARISON ---
    rti_mean = mean(RTI_clean, na.rm = TRUE),
    rti_sd = sd(RTI_clean, na.rm = TRUE),

    # --- MS2 CONFIRMATION ---
    labs_with_ms2 = n_distinct(lab_id[ms2_true_count > 0]),

    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    labs_missing = paste(
      setdiff(all_participating_labs, detected_labs),
      collapse = ", "
    ),

    # Some flags
    status_flag = case_when(
      n_labs == 1 & labs_with_ms2 == 0 ~ "WARNING: Single Lab / No MS2",
      n_labs < length(all_participating_labs) ~ "Check: Missing in some labs",
      rt_sd > 10 ~ "Check: High RT Variance", ## Is that too high ?
      TRUE ~ "Consensus OK"
    )
  ) %>%
  select(-detected_labs) %>%
  ungroup()

mz_stats <- df %>%
  select(target_ChEBI.name, starts_with("mz_")) %>%
  pivot_longer(
    cols = starts_with("mz_"),
    names_to = "adduct_type",
    values_to = "mz_value",
    values_drop_na = TRUE
  ) %>%
  group_by(target_ChEBI.name, adduct_type) %>%
  summarise(
    mz_mean = mean(mz_value, na.rm = TRUE),
    mz_sd = sd(mz_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = adduct_type,
    values_from = c(mz_mean, mz_sd),
    names_glue = "{.value}_{adduct_type}"
  )

consensus_table <- left_join(
  general_stats,
  mz_stats,
  by = "target_ChEBI.name"
) %>%
  arrange(status_flag, desc(n_labs))

print(head(consensus_table))

write.xlsx(
  consensus_table,
  file = "combine_annotation/consensus_summary.xlsx"
)

## Some labs need to chekc this in their raw data.
## and update the manual curation sheet (peak_evidence_rt_grouped.xlsx) and
## remove peak id that are not valid.
## how can they add new ones, they cannot add peak ids, make a new sheet ?
## i  think that's the easiest.
