# Reads raw DIAMOND (cycDB) output tsvs, counts gene hits per sample, maps
# genes to functional pathways, and caches the pathway x sample count matrix
# for scripts 8.2 and 8.3.
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
    "constants.R/helpers.R and the project's cycdb_diamond_output folder). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(stringr)
library(purrr)

source(file.path(script_dir, "constants.R"))
source(file.path(script_dir, "helpers.R"))

gene_map <- read_csv(gene_map_path, show_col_types = FALSE)
colnames(gene_map) <- c("ProteinID", "Gene", "database")

files <- list.files(diamond_dir, pattern = "_processed_diamond\\.tsv$", full.names = TRUE)

count_list <- lapply(files, function(file) {
  sample_name <- extract_sample_from_diamond_filename(file)
  if (is.na(sample_name)) {
    warning("Could not extract sample name from: ", basename(file), " -- skipping")
    return(NULL)
  }

  df <- read_tsv(file, col_names = diamond_column_names, show_col_types = FALSE)

  df_filtered <- df %>% filter(e_value <= e_value_threshold, identity >= identity_threshold)

  df_mapped <- df_filtered %>%
    inner_join(gene_map, by = c("protein" = "ProteinID"), relationship = "many-to-many") %>%
    select(Gene)

  df_mapped %>% dplyr::count(Gene, name = sample_name)
})
names(count_list) <- vapply(files, extract_sample_from_diamond_filename, character(1))
count_list <- count_list[!vapply(count_list, is.null, logical(1))]

# Gene x sample hit counts.
aggregated_counts <- purrr::reduce(count_list, full_join, by = "Gene") %>% replace(is.na(.), 0)

# Map each gene to its pathway, then sum gene counts up to the pathway level.
# Genes that don't belong to any tracked pathway (check_string() -> NA) are
# dropped by the inner_join against pathwaymap, same as the original script.
agg_counts_path <- aggregated_counts %>%
  mutate(pathway = vapply(Gene, check_string, character(1), lop = lop)) %>%
  select(pathway, everything(), -Gene) %>%
  inner_join(pathwaymap, by = "pathway") %>%
  select(-pathway) %>%
  rename(pathway = long_pathway)

agg_counts_path_summary <- agg_counts_path %>%
  group_by(pathway) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

write_csv(agg_counts_path_summary, file.path(diamond_dir, "pathway_raw_counts.csv"))
saveRDS(agg_counts_path_summary, pathway_counts_path)

agg_counts_path_summary
