# Genus-level PCoA (Bray-Curtis) ordination of the 16S community: one
# combined 5-panel figure (Methane, Oxygen, Year, depth, Site), plus a
# PERMANOVA (adonis2) test of those same variables. Mirrors
# step_5_pcoa/5.1_pcoa_ordination.R (kaiju version) so the two are directly
# comparable in the supplement.
# Saves PCoA_plot_genus_16s.png/.pdf.
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
    "constants.R and the project's supplementary folder). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
# -----------------------------------------------------------------------------

library(phyloseq)
library(ggplot2)
library(patchwork)
library(vegan)

source(file.path(script_dir, "constants.R"))  # phyloseq_16s_path, supplementary_dir

if (!file.exists(phyloseq_16s_path)) {
  stop(
    "Missing ", phyloseq_16s_path, ". Run scripts/step_4_16S_visualization/",
    "4.1_load_16s_data.R first to build the cached 16S phyloseq this ",
    "script reuses."
  )
}
physeq_16s <- readRDS(phyloseq_16s_path)

# Aggregate to genus and convert to relative abundance.
physeq_genus     <- tax_glom(physeq_16s, taxrank = "genus")
physeq_genus_rel <- transform_sample_counts(physeq_genus, function(x) x / sum(x))

# T0 (time-zero) samples have no Methane/Oxygen treatment applied yet;
# recode NA -> "T0" so they get their own ordination color/ellipse.
sample_data(physeq_genus_rel)$Year <- factor(sample_data(physeq_genus_rel)$Year)
sample_data(physeq_genus_rel)$Methane[is.na(sample_data(physeq_genus_rel)$Methane)] <- "T0"
sample_data(physeq_genus_rel)$Oxygen[is.na(sample_data(physeq_genus_rel)$Oxygen)]   <- "T0"

# ant_exp_map_MA.csv has a blank (not NA) Methane value for this sample. A
# blank string isn't caught by is.na() above, which is why this one point
# was showing up without a valid Methane group. (Oxygen was previously
# force-set to "No" here based on the sample's old, mislabeled name
# ("-O2"); ant_exp_map_MA.csv's own Oxygen="Yes" was correct all along --
# the name has been fixed upstream in step_6's sample_id_map instead.)
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
# Combined 5-panel figure: Methane / Oxygen / Year / depth / Site
# All five panels share the same ordination, so axes are left on their
# natural (auto) scale rather than manually cropped per panel.
# =============================================================================
pA <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Methane") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Methane), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("A")

pB <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Oxygen") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Oxygen), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("B")

pC <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Year") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Year), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("C")

pD <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "depth") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = depth), geom = "polygon", alpha = 0.15, level = 0.95) +
  labs(color = "Depth (cm)", fill = "Depth (cm)") +
  shared_theme + ggtitle("D")

pE <- plot_ordination(physeq_genus_rel, ord_pcoa, color = "Site") +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Site), geom = "polygon", alpha = 0.15, level = 0.95) +
  shared_theme + ggtitle("E")

layout_all <- (pA | pB | pC) / (pD | pE | plot_spacer())
layout_all

ggsave(
  filename = file.path(supplementary_dir, "FigureS1_PCoA_plot_genus_16s.png"),
  plot = layout_all, scale = 1, width = 15, height = 11, dpi = 600
)
ggsave(
  filename = file.path(supplementary_dir, "FigureS1_PCoA_plot_genus_16s.pdf"),
  plot = layout_all, scale = 1, width = 15, height = 11, dpi = 600, device = cairo_pdf
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
