# Genus-level PCoA (Bray-Curtis) ordination of the metagenomic (kaiju)
# community, faceted by Methane/Oxygen/Year and by depth/Site, plus a
# PERMANOVA (adonis2) test of those same variables.
# Saves PCoA_plot_genus_supplement.png/.pdf and
# Figure3_PCoA_plot_genus_depth_site.png/.pdf.
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed. Requires step_6_metagenomic_taxonomy/6.1_load_kaiju_data.R
# to have been run at least once, since this reuses its cached phyloseq
# rather than re-parsing the raw kaiju tsvs.

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
    "step_6_metagenomic_taxonomy/constants.R and the project's figures folder). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
taxonomy_dir <- file.path(project_root, "scripts", "step_6_metagenomic_taxonomy")
# -----------------------------------------------------------------------------

library(phyloseq)
library(ggplot2)
library(patchwork)
library(vegan)

source(file.path(taxonomy_dir, "constants.R"))  # kaiju_phyloseq_path, figures_dir

if (!file.exists(kaiju_phyloseq_path)) {
  stop(
    "Missing ", kaiju_phyloseq_path, ". Run scripts/step_6_metagenomic_taxonomy/",
    "6.1_load_kaiju_data.R first to build the cached kaiju phyloseq this ",
    "script reuses."
  )
}
kaiju_phyloseq <- readRDS(kaiju_phyloseq_path)

# Aggregate to genus and convert to relative abundance.
physeq_genus     <- tax_glom(kaiju_phyloseq, taxrank = "genus")
physeq_genus_rel <- transform_sample_counts(physeq_genus, function(x) x / sum(x))

# T0 (time-zero) samples have no Methane/Oxygen treatment applied yet;
# recode NA -> "T0" so they get their own ordination color/ellipse.
sample_data(physeq_genus_rel)$Year <- factor(sample_data(physeq_genus_rel)$Year)
sample_data(physeq_genus_rel)$Methane[is.na(sample_data(physeq_genus_rel)$Methane)] <- "T0"
sample_data(physeq_genus_rel)$Oxygen[is.na(sample_data(physeq_genus_rel)$Oxygen)]   <- "T0"

# This sample needed an explicit Methane fix here (see step_4's 4.3 for the
# same pattern). Oxygen was previously force-set to "No" here based on the
# sample's old, mislabeled name ("-O2"); metadata_ac_new.csv's own
# Oxygen="Yes" was correct all along -- the name has been fixed upstream in
# step_6's sample_id_map instead.
sample_data(physeq_genus_rel)["E2R1_6_9cm_+CH4_+O2", "Methane"] <- "Yes"

ord_pcoa <- ordinate(physeq_genus_rel, method = "PCoA", distance = "bray")

shared_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# =============================================================================
# Methane / Oxygen / Year panel
# =============================================================================
shared_coords <- coord_cartesian(xlim = c(-0.35, 0.45), ylim = c(-0.3, 0.3))

pA <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Methane") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Methane), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("A") + shared_coords

pB <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Oxygen") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Oxygen), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("B") + shared_coords

pC <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Year") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Year), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("C") + shared_coords

layout_treatment <- (pA | pB) / (pC + plot_spacer())
layout_treatment

ggsave(
  filename = file.path(project_root, "supplementary", "FigureS5_PCoA_plot_genus_supplement.png"),
  plot = layout_treatment, scale = 1, width = 10, height = 11, dpi = 600
)
ggsave(
  filename = file.path(project_root, "supplementary", "FigureS5_PCoA_plot_genus_supplement.pdf"),
  plot = layout_treatment, scale = 1, width = 10, height = 11, dpi = 600, device = cairo_pdf
)

# =============================================================================
# Depth / Site panel
# =============================================================================
shared_coords <- coord_cartesian(xlim = c(-0.3, 0.35), ylim = c(-0.25, 0.35))

pA <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "depth") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = depth), geom = "polygon", alpha = 0.15, level = 0.95) +
  labs(color = "Depth (cm)", fill = "Depth (cm)") +
  shared_theme + ggtitle("A") + shared_coords

pB <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Site") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Site), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("B") + shared_coords

layout_depth_site <- pA / pB
layout_depth_site

ggsave(
  filename = file.path(figures_dir, "Figure3_PCoA_plot_genus_depth_site.png"),
  plot = layout_depth_site, scale = 1, width = 10, height = 11, dpi = 600
)
ggsave(
  filename = file.path(figures_dir, "Figure3_PCoA_plot_genus_depth_site.pdf"),
  plot = layout_depth_site, scale = 1, width = 10, height = 11, dpi = 600, device = cairo_pdf
)

# =============================================================================
# PERMANOVA
# =============================================================================
dist_bray <- phyloseq::distance(physeq_genus_rel, method = "bray")
metadata  <- data.frame(sample_data(physeq_genus_rel))

metadata$Year    <- factor(metadata$Year)
metadata$Methane <- factor(metadata$Methane)
metadata$Oxygen  <- factor(metadata$Oxygen)
metadata$depth   <- factor(metadata$depth)
metadata$Site    <- factor(metadata$Site)

# Year dropped: depth is more biologically informative and the two are
# largely confounded (2023 samples use a different depth scheme than 2022).
adonis2(dist_bray ~ depth + Methane + Oxygen + Site, data = metadata, by = "margin", permutations = 999)
