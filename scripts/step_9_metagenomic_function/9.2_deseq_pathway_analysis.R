# DESeq2 differential abundance testing of functional pathway counts:
# full-dataset model (Year + Site + sediment + Methane) plus within-year
# (2022, 2023) depth + methane models. Writes deseq2_*.csv result tables and
# lollipop plots of significant pathways per contrast.
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed. Requires 9.1_build_pathway_counts.R to have been run at
# least once.

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
    "constants.R/helpers.R and the project's cycdb_diamond_output/deseq_output folders). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(DESeq2)
library(ggplot2)

source(file.path(script_dir, "constants.R"))
source(file.path(script_dir, "helpers.R"))

loaded  <- load_pathway_matrix_and_coldata(pathway_counts_path, exp_metadata_path, sample_id_map)
mat     <- loaded$mat
coldata <- loaded$coldata

# =============================================================================
# FULL DATASET: Year, Site, sediment depth, Methane
# (Oxygen not testable -- confounded with depth in 2023, not manipulated in 2022)
# =============================================================================
keep_full <- complete.cases(coldata[, c("Year", "Site", "sediment", "Methane")])
mat_full     <- mat[, keep_full]
coldata_full <- droplevels(coldata[keep_full, ])

dds_full <- DESeqDataSetFromMatrix(
  countData = mat_full,
  colData   = coldata_full,
  design    = ~ Year + Site + sediment + Methane
)
dds_full <- estimateSizeFactors(dds_full, type = "poscounts")
dds_full <- DESeq(dds_full)

res_sediment <- results(dds_full, name = "sediment_deep_vs_surface")
res_methane  <- results(dds_full, name = "Methane_Yes_vs_No")
res_site     <- results(dds_full, name = "Site_Jetty_vs_CCS")
res_year     <- results(dds_full, name = "Year_23_vs_22")

# 2023 within-year depth + methane effects
keep_23 <- coldata$Year == "23" & complete.cases(coldata[, c("sediment", "Methane")])
mat_23     <- mat[, keep_23]
coldata_23 <- droplevels(coldata[keep_23, ])

dds_23 <- DESeqDataSetFromMatrix(
  countData = mat_23,
  colData   = coldata_23,
  design    = ~ Methane + sediment  # can't add Oxygen due to confound
)
dds_23 <- estimateSizeFactors(dds_23, type = "poscounts")
dds_23 <- DESeq(dds_23)

res_23_sed  <- results(dds_23, name = "sediment_deep_vs_surface")
res_23_meth <- results(dds_23, name = "Methane_Yes_vs_No")

# 2022 within-year depth + methane effects
keep_22 <- coldata$Year == "22" & complete.cases(coldata[, c("sediment", "Methane")])
mat_22     <- mat[, keep_22]
coldata_22 <- droplevels(coldata[keep_22, ])

dds_22 <- DESeqDataSetFromMatrix(
  countData = mat_22,
  colData   = coldata_22,
  design    = ~ Methane + sediment
)
dds_22 <- estimateSizeFactors(dds_22, type = "poscounts")
dds_22 <- DESeq(dds_22)

res_22_sed  <- results(dds_22, name = "sediment_deep_vs_surface")
res_22_meth <- results(dds_22, name = "Methane_Yes_vs_No")

write_deseq_results(res_sediment, "full_sediment_effect", deseq_output_dir)
write_deseq_results(res_methane,  "full_methane_effect",  deseq_output_dir)
write_deseq_results(res_site,     "full_site_effect",     deseq_output_dir)
write_deseq_results(res_year,     "full_year_effect",     deseq_output_dir)
write_deseq_results(res_23_sed,   "sediment_effect_2023_only", deseq_output_dir)
write_deseq_results(res_23_meth,  "methane_effect_2023_only",  deseq_output_dir)
write_deseq_results(res_22_sed,   "sediment_effect_2022_only", deseq_output_dir)
write_deseq_results(res_22_meth,  "methane_effect_2022_only",  deseq_output_dir)

# =============================================================================
# Lollipop plots of significant pathways per contrast
# =============================================================================
df_sediment <- tidy_deseq_results(res_sediment, "Deep vs Surface")
df_methane  <- tidy_deseq_results(res_methane,  "+CH4 vs -CH4")
df_site     <- tidy_deseq_results(res_site,     "Jetty vs CCS")
df_year     <- tidy_deseq_results(res_year,     "2023 vs 2022")

p_loll_sed  <- lollipop_plot(df_sediment, "Significant pathways: Deep vs Surface")
p_loll_meth <- lollipop_plot(df_methane,  "Significant pathways: +CH4 vs -CH4")
p_loll_site <- lollipop_plot(df_site,     "Significant pathways: Jetty vs CCS")
p_loll_year <- lollipop_plot(df_year,     "Significant pathways: 2023 vs 2022")

ggsave(file.path(deseq_output_dir, "lollipop_year.png"),     p_loll_year, width = 8, height = 8, dpi = 300)
ggsave(file.path(deseq_output_dir, "lollipop_sediment.png"), p_loll_sed,  width = 8, height = 8, dpi = 300)
ggsave(file.path(deseq_output_dir, "lollipop_methane.png"),  p_loll_meth, width = 8, height = 6, dpi = 300)
ggsave(file.path(deseq_output_dir, "lollipop_site.png"),     p_loll_site, width = 8, height = 6, dpi = 300)
ggsave(file.path(deseq_output_dir, "lollipop_year.pdf"),     p_loll_year, width = 8, height = 8, dpi = 300)
ggsave(file.path(deseq_output_dir, "lollipop_sediment.pdf"), p_loll_sed,  width = 8, height = 8, dpi = 300)
ggsave(file.path(deseq_output_dir, "lollipop_methane.pdf"),  p_loll_meth, width = 8, height = 6, dpi = 300)
ggsave(file.path(deseq_output_dir, "lollipop_site.pdf"),     p_loll_site, width = 8, height = 6, dpi = 300)
