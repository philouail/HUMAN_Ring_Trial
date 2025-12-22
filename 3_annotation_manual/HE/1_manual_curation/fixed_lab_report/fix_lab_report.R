# ==============================================================================
# SCRIPT PURPOSE: Rebuild Grouped Data from Manual Validation
# ==============================================================================
#
# WHY THIS IS NEEDED:
# The manual validation file contains the 'chrom_peak_id's that the student/reviewer
# marked as TRUE (valid), but for some columns, the full aggregated metrics is
# now inaccurate. We need to go back to the raw 'peak_evidence' file to
# retrieve the data and re-summarize it.
#
# SPECIFIC LOGIC APPLIED:
# 1. Filter: Only process rows where Validation == TRUE.
# 2. Retrieve: Look up the component peaks in the raw 'peak_evidence.csv'.
# 3. Deduplicate Adducts (The Fix):
#    Sometimes the same adduct (e.g., [M+H]+) is assigned to multiple peaks
#    within the same group. To prevent double-counting or using poor data:
#    -> For every unique adduct, we keep ONLY the row with the lowest ppm_error.
# 4. Aggregation: Recalculate means (rt, intensity) and sums (counts).
#    Note: m/z is NOT averaged; specific m/z values are pivoted to new columns.
# 5. Output: Save the corrected grouped data to a new Excel file.
# ==============================================================================

library(openxlsx)
library(dplyr) # Required for bind_rows

process_lab_evidence <- function(lab_name) {
  # 1. Setup Paths
  pe_path <- file.path("2_annotation_auto", lab_name, "HE", "peak_evidence.csv")
  pg_path <- file.path(
    "3_annotation_manual",
    "HE",
    "1_manual_curation",
    "lab_report",
    paste0("annotation_", lab_name, ".xlsx")
  )

  # Output: Saves to current directory with lab name appended
  out_path <- file.path(
    "3_annotation_manual",
    "HE",
    "1_manual_curation",
    "fixed_lab_report",
    paste0("fixed_annotation_", lab_name, ".xlsx")
  )

  message(paste0("Processing: ", lab_name))

  # 2. Validation Checks
  if (!file.exists(pe_path) || !file.exists(pg_path)) {
    warning(paste("  -> SKIPPING: Files not found for", lab_name))
    return(NULL)
  }

  pe <- read.csv(pe_path)
  pg <- read.xlsx(pg_path)

  if ("Validation" %in% colnames(pg)) {
    pg$Validation[is.na(pg$Validation)] <- FALSE
    pg <- pg[pg$Validation == TRUE, ]
  } else {
    warning(paste("  -> SKIPPING: 'Validation' column missing in", lab_name))
    return(NULL)
  }

  if (nrow(pg) == 0) {
    message("  -> No validated rows found.")
    return(NULL)
  }

  # 3. Process Rows
  res_list <- lapply(1:nrow(pg), function(i) {
    # Extract IDs
    target_ids <- unlist(strsplit(pg$chrom_peak_id[i], "\\|"))
    x_res <- pe[pe$chrom_peak_id %in% target_ids, ]

    # --- FIX: Keep lowest PPM for duplicate adducts ---
    sub <- do.call(
      rbind,
      lapply(split(x_res, x_res$adduct), function(l) {
        l[l$ppm_error == min(l$ppm_error), ]
      })
    )

    # Remove columns to be recalculated or pivoted
    cols_to_remove <- c(
      "adduct_formula",
      "adduct_mz",
      "score",
      "ppm_error",
      "rt",
      "mz",
      "polarity",
      "mzmin",
      "mzmax"
    )
    res <- sub[1, !(colnames(sub) %in% cols_to_remove)]

    # --- PIVOT: Create columns for each mz per adduct ---
    for (r in 1:nrow(sub)) {
      col_name <- paste0("mz_", sub$adduct[r])
      res[[col_name]] <- sub$mz[r]
    }

    # Summarise Group Metrics
    res$chrom_peak_id <- paste(sub$chrom_peak_id, collapse = "|")
    res$adduct <- paste(sub$adduct, collapse = "|")
    res$rtmed <- mean(sub$rt)
    res$rtmin <- min(sub$rtmin)
    res$rtmax <- max(sub$rtmax)
    res$into <- mean(sub$into)
    res$ms1_adduct_count <- sum(sub$ms1_adduct_count)
    res$ms1_adduct_05_count <- sum(sub$ms1_adduct_05_count)
    res$isopeak_count <- sum(sub$isopeak_count)
    res$isopeak_sim <- mean(sub$isopeak_sim, na.rm = TRUE)
    res$ms2_count <- sum(sub$ms2_count, na.rm = TRUE)
    res$ms2_matched_count <- sum(sub$ms2_matched_count, na.rm = TRUE)
    res$ms2_true_count <- sum(sub$ms2_true_count, na.rm = TRUE)
    res$number_of_peaks <- nrow(sub)
    res$RTI <- paste(sub$RTI, collapse = "|")

    return(res)
  })

  # 4. Bind and Save
  # bind_rows aligns columns, filling NA where an adduct column is missing for a specific group
  res_all <- dplyr::bind_rows(res_list)

  write.xlsx(res_all, file = out_path)
  message(paste("  -> Done. Saved to:", out_path))
}

# ==============================================================================
# EXECUTION
# ==============================================================================

process_lab_evidence("hmgu")
process_lab_evidence("icl")
process_lab_evidence("afekta")
process_lab_evidence("cembio")
