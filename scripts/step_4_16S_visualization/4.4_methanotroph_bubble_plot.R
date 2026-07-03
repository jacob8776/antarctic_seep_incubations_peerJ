# Bubble plot of methanotroph/methanotroph-adjacent genera relative
# abundance across 16S samples, faceted by year, with location + metabolism
# annotation bars. Mirrors step_7_methanotroph_visualization/
# 7.1_methanotroph_bubble_plot.R (kaiju version) so the two are directly
# comparable in the supplement -- genus names are reconciled to canonical
# form first since 16S/SILVA naming differs from kaiju/NCBI for a couple of
# these taxa (see step_4/helpers.R).
# Saves methanotrophs_genus_16s.png/.pdf.
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
    "constants.R/helpers.R, step_8's reconcile_phyloseq.R, and the project's ",
    "supplementary folder). Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir      <- .get_script_dir()
project_root    <- normalizePath(file.path(script_dir, "..", ".."))
taxonomy_dir    <- file.path(project_root, "scripts", "step_6_metagenomic_taxonomy")
comparison_dir  <- file.path(project_root, "scripts", "step_8_comparing_16s_metagenomics")
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)

source(file.path(taxonomy_dir, "constants.R"))                  # sample_list (shared 34-sample ordering)
source(file.path(script_dir, "constants.R"))                    # phyloseq_16s_path, supplementary_dir
source(file.path(script_dir, "helpers.R"))                      # taxa_of_interest_16s, genus_order_16s, load_relative_abundance_long_16s
source(file.path(comparison_dir, "reconcile_phyloseq.R"))       # reconcile_phyloseq(), fix_lineage_from_genus()

loaded <- load_relative_abundance_long_16s(phyloseq_16s_path, taxa_of_interest_16s, sample_list, genus_order_16s)
relative_abundance_long <- loaded$data
sample_order <- loaded$sample_order

p <- ggplot(relative_abundance_long, aes(x = Sample, y = genus)) +
  geom_tile(aes(fill = Location), alpha = 0, show.legend = TRUE) +
  geom_point(aes(size = RelativeAbundance, color = RelativeAbundance)) +
  scale_size(range = c(1, 12), name = "Relative\nAbundance") +
  scale_color_viridis_c(option = "mako", direction = -1, name = "Relative\nAbundance") +
  scale_fill_manual(values = c("CCS" = "#1B9E77", "Jetty" = "#F0C808"), name = "Location") +
  guides(
    size = guide_legend(order = 1),
    color = guide_legend(order = 1),
    fill = guide_legend(order = 2, override.aes = list(alpha = 1))
  ) +
  facet_grid(~ Year, scales = "free_x", space = "free_x") +
  labs(x = "", y = "") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 13),
    axis.text.y = element_text(face = "italic", size = 15),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    legend.key = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

loc_data <- relative_abundance_long %>%
  distinct(Sample, Location, Year) %>%
  mutate(Sample = factor(Sample, levels = sample_order))

loc_bar <- ggplot(loc_data, aes(x = Sample, y = 1, fill = Location)) +
  geom_tile(width = 1) +
  scale_fill_manual(values = c("CCS" = "#1B9E77", "Jetty" = "#F0C808")) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_grid(~ Year, scales = "free_x", space = "free_x") +
  theme_void() +
  theme(legend.position = "none", strip.text = element_blank(), plot.margin = margin(0, 5, 0, 5))

metab_data <- relative_abundance_long %>%
  distinct(genus, Metabolism) %>%
  mutate(genus = factor(genus, levels = genus_order_16s))

metab_bar <- ggplot(metab_data, aes(x = 1, y = genus, fill = Metabolism)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = c("Aerobic" = "#2166AC", "Anaerobic" = "#B2182B"), name = "Metabolism") +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "cm"),
    plot.margin = margin(5, 0, 5, 0)
  )

final_plot <- (plot_spacer() + loc_bar + metab_bar + p) +
  plot_layout(widths = c(1, 20), heights = c(1, 20))

ggsave(
  filename = file.path(supplementary_dir, "FigureS4_methanotrophs_genus_16s.png"),
  plot = final_plot, width = 16, height = 9, dpi = 600
)
ggsave(
  filename = file.path(supplementary_dir, "FigureS4_methanotrophs_genus_16s.pdf"),
  plot = final_plot, width = 16, height = 9, dpi = 600, device = cairo_pdf
)
