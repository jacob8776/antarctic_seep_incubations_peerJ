# Top-15 genus stacked bar plots for Bacteria and Archaea, plus highlighted
# variants for sulfate reducers, Methanoperedens, and methanogens.
# Saves bac_rank_met.png and ar_rank_met.png.
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed. Requires 6.1_load_kaiju_data.R to have been run at least once.

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
    "constants.R/helpers.R and the project's kaiju_output/figures folders). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
# -----------------------------------------------------------------------------

library(dplyr)
library(phyloseq)
library(ggplot2)

source(file.path(script_dir, "constants.R"))
source(file.path(script_dir, "helpers.R"))

kaiju_phyloseq <- readRDS(kaiju_phyloseq_path)

kaiju_phyloseq_ar  <- subset_taxa(kaiju_phyloseq, superkingdom %in% c("Archaea"))
kaiju_phyloseq_bac <- subset_taxa(kaiju_phyloseq, superkingdom %in% c("Bacteria"))

top15_bac_rank <- get_top_taxa_by_rank(kaiju_phyloseq_bac, rank = "genus", top_n = 15)
top15_ar_rank  <- get_top_taxa_by_rank(kaiju_phyloseq_ar,  rank = "genus", top_n = 15)

# --- Summary for the stacked bar-plot paragraph -----------------------------
total_reads <- sum(otu_table(kaiju_phyloseq))
bac_reads   <- sum(otu_table(kaiju_phyloseq_bac))
ar_reads    <- sum(otu_table(kaiju_phyloseq_ar))

bac_genera <- sum(unique(as.character(tax_table(kaiju_phyloseq_bac)[, "genus"])) != "Unclassified")
ar_genera  <- sum(unique(as.character(tax_table(kaiju_phyloseq_ar)[,  "genus"])) != "Unclassified")

top15_bac_frac <- mean(sample_sums(top15_bac_rank))  # mean across samples
top15_ar_frac  <- mean(sample_sums(top15_ar_rank))

cat(sprintf(
  "Total classified reads: %s
  Bacteria: %s (%.1f%%) | %d genera
  Archaea:  %s (%.1f%%) | %d genera
Top-15 genera capture: Bacteria %.1f%%, Archaea %.1f%% of their domain (mean/sample)\n",
  format(total_reads, big.mark = ","),
  format(bac_reads, big.mark = ","), 100 * bac_reads / total_reads, bac_genera,
  format(ar_reads,  big.mark = ","), 100 * ar_reads  / total_reads, ar_genera,
  100 * top15_bac_frac, 100 * top15_ar_frac
))

# label_genus_with_phylum() is defined in helpers.R
top15_bac_rank <- label_genus_with_phylum(top15_bac_rank)
top15_ar_rank  <- label_genus_with_phylum(top15_ar_rank)

p_bac_rank <- plot_top_taxa_bar(top15_bac_rank, sample_order = sample_list, rank = "genus", title = "")
p_ar_rank  <- plot_top_taxa_bar(top15_ar_rank,  sample_order = sample_list, rank = "genus", title = "")

print(p_bac_rank)
print(p_ar_rank)

sulfate_reducers <- plot_top_taxa_bar(
  top15_bac_rank, sample_order = sample_list, rank = "genus", title = "",
  highlighted_taxa = c(
    "Desulfogranum (Thermodesulfobacteriota)", "Desulfomarina (Thermodesulfobacteriota)",
    "Desulfonema (Thermodesulfobacteriota)", "Desulfopila (Thermodesulfobacteriota)",
    "Desulforhopalus (Thermodesulfobacteriota)", "Desulfosarcina (Thermodesulfobacteriota)",
    "Desulfosediminicola (Thermodesulfobacteriota)"
  )
)

methanoperedens <- plot_top_taxa_bar(
  top15_ar_rank, sample_order = sample_list, rank = "genus", title = "",
  highlighted_taxa = c("Candidatus Methanoperedens (Euryarchaeota)")
)

methanogens <- plot_top_taxa_bar(
  top15_ar_rank, sample_order = sample_list, rank = "genus", title = "",
  highlighted_taxa = c(
    "Methanobacterium (Euryarchaeota)", "Methanocella (Euryarchaeota)",
    "Methanococcoides (Euryarchaeota)", "Methanoculleus (Euryarchaeota)",
    "Methanofollis (Euryarchaeota)", "Methanohalophilus (Euryarchaeota)",
    "Methanolobus (Euryarchaeota)", "Methanosarcina (Euryarchaeota)",
    "Methanothrix (Euryarchaeota)", "Methanobrevibacter (Euryarchaeota)"
  )
)

ggsave(file.path(figures_dir, "Figure4A_bac_rank_met.png"), plot = p_bac_rank, scale = 1, width = 10, height = 7, dpi = 600)
ggsave(file.path(figures_dir, "Figure4B_ar_rank_met.png"),  plot = p_ar_rank,  scale = 1, width = 10, height = 7, dpi = 600)
ggsave(file.path(figures_dir, "Figure4A_bac_rank_met.pdf"), plot = p_bac_rank, scale = 1, width = 10, height = 7, dpi = 600)
ggsave(file.path(figures_dir, "Figure4B_ar_rank_met.pdf"),  plot = p_ar_rank,  scale = 1, width = 10, height = 7, dpi = 600)
