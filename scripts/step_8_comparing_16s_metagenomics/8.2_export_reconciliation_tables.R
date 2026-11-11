# Exports the taxonomy reconciliation tables from reconcile_phyloseq.R
# (database synonym mapping + genus-to-lineage corrections) as CSVs for the
# supplement, so the reconciliation applied in
# 8.1_metagenome_vs_16s_comparison.R is documented and auditable.
#
# supplementary_taxonomy_reconciliation_synonyms.csv: database-specific
# taxonomy strings (SILVA 16S rRNA, NCBI/Kaiju, GTDB) mapped to a single
# canonical name for cross-method comparison. Rows where the canonical name
# equals the database name are included for completeness.
#
# supplementary_taxonomy_reconciliation_lineage_corrections.csv: order and
# family assignments filled in or corrected from genus identity when a
# database left these ranks unclassified or misassigned (e.g., Kaiju/NCBI's
# sparse lineage for candidate division NC10). Applied by
# fix_lineage_from_genus() after name reconciliation.
#
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed.

# --- Resolve project paths relative to this script's own location ----------
.get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_flag) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_flag))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(ctx)) return(dirname(normalizePath(ctx)))
  }
  stop(
    "Could not determine this script's location automatically (needed to find ",
    "reconcile_phyloseq.R and the project's supplementary folder). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
# -----------------------------------------------------------------------------

library(dplyr)

# taxonomy_map / genus_to_correct_lineage
source(file.path(script_dir, "reconcile_phyloseq.R"))

supplementary_dir <- file.path(project_root, "supplementary")
dir.create(supplementary_dir, showWarnings = FALSE)

rank_order <- c("phylum", "order", "family", "genus")

synonyms_table <- taxonomy_map %>%
  arrange(factor(rank, levels = rank_order), canonical_name, db_source) %>%
  transmute(
    Rank             = rank,
    Database         = db_source,
    `Database name`  = db_name,
    `Canonical name` = canonical_name
  )

lineage_table <- genus_to_correct_lineage %>%
  arrange(canonical_genus) %>%
  transmute(
    `Canonical genus`  = canonical_genus,
    `Corrected order`  = correct_order,
    `Corrected family` = correct_family
  )

synonyms_path <- file.path(supplementary_dir, "TableS3_supplementary_taxonomy_reconciliation_synonyms.csv")
lineage_path  <- file.path(supplementary_dir, "TableS4_supplementary_taxonomy_reconciliation_lineage_corrections.csv")

write.csv(synonyms_table, synonyms_path, row.names = FALSE)
write.csv(lineage_table,  lineage_path,  row.names = FALSE)

message("Saved: ", synonyms_path)
message("Saved: ", lineage_path)
