# Top-15 genus stacked bar plots for 16S Bacteria and Archaea -- same format
# (plot_top_taxa_bar(), location annotation bar, phylum-labeled legend) as
# step_6's metagenomic (kaiju) version, so the two are directly comparable.
# Saves bac_rank_16s.png and ar_rank_16s.png.
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed. Requires 4.1_load_16s_data.R to have been run at least once.

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
    "constants.R, step_6's constants.R/helpers.R, and the project's ",
    "updated_16s_run/supplementary folders). Run this script with Rscript, or ",
    "open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
taxonomy_dir <- file.path(project_root, "scripts", "step_6_metagenomic_taxonomy")
# -----------------------------------------------------------------------------

library(dplyr)
library(phyloseq)
library(ggplot2)

source(file.path(taxonomy_dir, "constants.R"))  # figures_dir, sample_list
source(file.path(taxonomy_dir, "helpers.R"))    # get_top_taxa_by_rank, plot_top_taxa_bar, label_genus_with_phylum
source(file.path(script_dir, "constants.R"))    # phyloseq_16s_bac_path, phyloseq_16s_ar_path

physeq_16s_bac <- readRDS(phyloseq_16s_bac_path)
physeq_16s_ar  <- readRDS(phyloseq_16s_ar_path)

top15_bac_rank_16s <- get_top_taxa_by_rank(physeq_16s_bac, rank = "genus", top_n = 15)
top15_ar_rank_16s  <- get_top_taxa_by_rank(physeq_16s_ar,  rank = "genus", top_n = 15)

top15_bac_rank_16s <- label_genus_with_phylum(top15_bac_rank_16s)
top15_ar_rank_16s  <- label_genus_with_phylum(top15_ar_rank_16s)

p_bac_rank_16s <- plot_top_taxa_bar(top15_bac_rank_16s, sample_order = sample_list, rank = "genus", title = "")
p_ar_rank_16s  <- plot_top_taxa_bar(top15_ar_rank_16s,  sample_order = sample_list, rank = "genus", title = "")

print(p_bac_rank_16s)
print(p_ar_rank_16s)

ggsave(file.path(supplementary_dir, "FigureS2_bac_rank_16s.png"), plot = p_bac_rank_16s, scale = 1, width = 10, height = 7, dpi = 600)
ggsave(file.path(supplementary_dir, "Figure_S3_ar_rank_16s.png"),  plot = p_ar_rank_16s,  scale = 1, width = 10, height = 7, dpi = 600)
ggsave(file.path(supplementary_dir, "FigureS2_bac_rank_16s.pdf"), plot = p_bac_rank_16s, scale = 1, width = 10, height = 7, dpi = 600)
ggsave(file.path(supplementary_dir, "FigureS3_ar_rank_16s.pdf"),  plot = p_ar_rank_16s,  scale = 1, width = 10, height = 7, dpi = 600)
