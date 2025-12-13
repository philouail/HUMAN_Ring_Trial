# ==============================================================================
# SCRIPT PURPOSE: Rebuild Grouped Data from Manual Validation
# ==============================================================================
#
# WHY THIS IS NEEDED:
# The manual validation file contains the 'chrom_peak_id's that the student/reviewer
# marked as TRUE (valid), but the full aggregated metrics is now innacurate. 
# We need to go back to the raw 'peak_evidence' file to
# retrieve the data and re-summarize it.
#
# SPECIFIC LOGIC APPLIED:
# 1. Filter: Only process rows where Validation == TRUE.
# 2. Retrieve: Look up the component peaks in the raw 'peak_evidence.csv'.
# 3. Deduplicate Adducts (The Fix):
#    Sometimes the same adduct (e.g., [M+H]+) is assigned to multiple peaks
#    within the same group. To prevent double-counting or using poor data:
#    -> For every unique adduct, we keep ONLY the row with the lowest ppm_error.
# 4. Aggregation: Recalculate means (rt, mz, intensity) and sums (counts)
#    based on this cleaned subset of peaks.
# ==============================================================================

# Note to self: need to re-adjust the lib_gen file afterward.and adjust read me and overall workflow explanbation with this. 

library(openxlsx)

process_lab_evidence <- function(lab_name, manual_filename) {
  # 1. Construct file paths dynamically
  #    Using file.path handles slash direction automatically for OS compatibility
  base_dir <- file.path(lab_name, "results", "HE")
  pe_path <- file.path(base_dir, "peak_evidence.csv")
  pg_path <- file.path(base_dir, manual_filename)
  out_path <- file.path(base_dir, "peak_evidence_rt_grouped_manual_fixed.xlsx")

  message(paste0("Processing: ", lab_name))

  # 2. Check if files exist
  if (!file.exists(pe_path) || !file.exists(pg_path)) {
    warning(paste("  -> SKIPPING: Files not found for", lab_name))
    return(NULL)
  }

  # 3. Load Data
  pe <- read.csv(pe_path)
  pg <- read.xlsx(pg_path)

  # 4. Filter for Validated rows only
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

  # 5. Process each validated group
  res_all <- do.call(
    rbind,
    lapply(1:nrow(pg), function(i) {
      # Extract the IDs associated with this group (split by pipe if necessary)
      target_ids <- unlist(strsplit(pg$chrom_peak_id[i], "\\|"))
      x_res <- pe[pe$chrom_peak_id %in% target_ids, ]

      # --- CRITICAL STEP: Keep lowest PPM for duplicates ---
      sub <- do.call(
        rbind,
        lapply(split(x_res, x_res$adduct), function(l) {
          l[l$ppm_error == min(l$ppm_error), ]
        })
      )

      # Remove columns that will be re-calculated or aren't needed in the summary
      cols_to_remove <- c(
        "adduct_formula",
        "adduct_mz",
        "score",
        "ppm_error",
        "rt",
        "mz",
        "polarity"
      )
      res <- sub[1, !(colnames(sub) %in% cols_to_remove)]

      # Re-summarise the data based on the 'sub' dataframe
      res$chrom_peak_id <- paste(sub$chrom_peak_id, collapse = "|")
      res$adduct <- paste(sub$adduct, collapse = "|")
      res$mzmin <- min(sub$mzmin)
      res$mzmax <- max(sub$mzmax)
      res$mzmed <- mean(sub$mz)
      res$rtmin <- min(sub$rtmin)
      res$rtmax <- max(sub$rtmax)
      res$rtmed <- mean(sub$rt)
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
  ) 

  # 6. Save Result
  write.xlsx(res_all, file = out_path)
  message(paste("  -> Done. Saved to:", out_path))
}

# ==============================================================================
# APPLY TO LABS
# ==============================================================================

process_lab_evidence("hmgu", "peak_evidence_rt_grouped_manual_Pualine.xlsx")
process_lab_evidence("icl", "peak_evidence_rt_grouped_curated_ICL.xlsx")
process_lab_evidence("afekta", "peak_evidence_rt_grouped_he_Dennisse.xlsx")
process_lab_evidence("cembio", "peak_evidence_rt_groupedTF.xlsx")
