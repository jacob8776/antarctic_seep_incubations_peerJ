# VST-normalized heatmaps of functional pathway abundance across all
# samples (pheatmap + ComplexHeatmap versions), annotated by Year, Site,
# sediment, depth, Methane, and Oxygen.
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
    "constants.R/helpers.R and the project's cycdb_diamond_output/figures folders). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
# -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(tibble)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

source(file.path(script_dir, "constants.R"))
source(file.path(script_dir, "helpers.R"))

loaded  <- load_pathway_matrix_and_coldata(pathway_counts_path, exp_metadata_path, sample_id_map)
mat     <- loaded$mat
coldata <- loaded$coldata

dds_all <- DESeqDataSetFromMatrix(countData = mat, colData = coldata, design = ~1)
dds_all <- estimateSizeFactors(dds_all, type = "poscounts")
vsd_all <- varianceStabilizingTransformation(dds_all, blind = TRUE)
vst_mat <- assay(vsd_all)

col_annot <- coldata[, c("Year", "Site", "sediment", "depth", "Methane", "Oxygen"), drop = FALSE]

ann_colors <- list(
  Year     = c(`22` = "#b2df8a", `23` = "#33a02c"),
  Site     = c(CCS = "#6a3d9a", Jetty = "#cab2d6"),
  sediment = c(surface = "#fdbf6f", deep = "#ff7f00"),
  depth    = c(`0_3` = "#fee5d9", `0_4` = "#fcae91",
               `3_6` = "#fb6a4a", `4_8` = "#de2d26", `6_9` = "#a50f15"),
  Methane  = c(No = "#cccccc", Yes = "#000000"),
  Oxygen   = c(No = "#cccccc", Yes = "#1f78b4")
)

# T0 (time-zero) samples have no Methane/Oxygen treatment applied yet;
# recode NA -> "T0" so they get their own annotation color instead of
# rendering blank.
col_annot$Methane <- as.character(col_annot$Methane)
col_annot$Oxygen  <- as.character(col_annot$Oxygen)
col_annot$Methane[is.na(col_annot$Methane)] <- "T0"
col_annot$Oxygen[is.na(col_annot$Oxygen)]   <- "T0"
col_annot$Methane <- factor(col_annot$Methane, levels = c("No", "Yes", "T0"))
col_annot$Oxygen  <- factor(col_annot$Oxygen,  levels = c("No", "Yes", "T0"))

# ant_exp_map_MA.csv has a data-entry gap for this sample (Methane blank,
# Oxygen wrong); its own display name says +CH4/-O2, so fix it here for the
# annotation bar. Not applied upstream in 8.2 -- that sample was already
# excluded from the DESeq2 models by the Methane completeness filter there.
col_annot["E2R1_6_9cm_+CH4_-O2", "Methane"] <- "Yes"
col_annot["E2R1_6_9cm_+CH4_-O2", "Oxygen"]  <- "No"

ann_colors$Methane <- c(No = "#cccccc", Yes = "#000000", T0 = "#8856a7")
ann_colors$Oxygen  <- c(No = "#cccccc", Yes = "#1f78b4", T0 = "#8856a7")

#pheatmap::pheatmap(
#  vst_mat,
#  angle_col                = "90",
#  scale                    = "row",
#  cluster_rows             = TRUE,
#  cluster_cols             = TRUE,
#  clustering_distance_rows = "correlation",
#  clustering_distance_cols = "euclidean",
#  clustering_method        = "ward.D2",
#  annotation_col           = col_annot,
#  annotation_colors        = ann_colors,
#  color                    = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
#  fontsize_row             = 10,
#  fontsize_col             = 10,
#  res                      = 300,
#  filename                 = file.path(figures_dir, "pathway_heatmap_pheatmap.png"),
#  width                    = 14,
#  height                   = 11
#)

vst_mat_scaled <- t(scale(t(vst_mat)))  # equivalent to pheatmap's scale = "row"

ht <- Heatmap(
  vst_mat_scaled,
  name                        = "z-score",
  col                         = colorRamp2(c(-2, 0, 2), rev(brewer.pal(3, "RdBu"))),
  cluster_rows                = TRUE,
  cluster_columns             = TRUE,
  clustering_distance_rows    = "pearson",
  clustering_distance_columns = "euclidean",
  clustering_method_rows      = "ward.D2",
  clustering_method_columns   = "ward.D2",
  rect_gp                     = gpar(col = "black", lwd = 0.5),
  top_annotation              = HeatmapAnnotation(
    df  = col_annot,
    col = ann_colors,
    annotation_name_gp      = gpar(fontsize = 17),
    annotation_legend_param = list(
      labels_gp = gpar(fontsize = 17),
      title_gp  = gpar(fontsize = 17)
    )
  ),
  row_names_gp             = gpar(fontsize = 17),
  row_names_max_width      = max_text_width(rownames(vst_mat_scaled), gp = gpar(fontsize = 17)),
  column_names_gp          = gpar(fontsize = 17),
  column_names_rot         = 90,
  column_names_max_height  = unit(10, "cm"),
  heatmap_legend_param      = list(
    direction    = "horizontal",
    labels_gp    = gpar(fontsize = 17),
    title_gp     = gpar(fontsize = 17),
    legend_width = unit(6, "cm")
  )
)

png(file.path(figures_dir, "Figure8_pathway_heatmap.png"), width = 20, height = 16, units = "in", res = 300)
draw(
  ht,
  heatmap_legend_side    = "bottom",
  annotation_legend_side = "bottom",
  padding                = unit(c(10, 80, 2, 2), "mm")
)
dev.off()

pdf(file.path(figures_dir, "Figure8_pathway_heatmap.pdf"), width = 20, height = 16)
draw(
  ht,
  heatmap_legend_side    = "bottom",
  annotation_legend_side = "bottom",
  padding                = unit(c(10, 80, 2, 2), "mm")
)
dev.off()

write_csv(as.data.frame(vst_mat) %>% rownames_to_column("pathway"), file.path(diamond_dir, "pathway_vst_matrix.csv"))
write_csv(coldata %>% rownames_to_column("sample"), file.path(diamond_dir, "sample_metadata_used.csv"))
