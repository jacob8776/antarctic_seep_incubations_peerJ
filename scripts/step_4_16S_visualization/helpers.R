# Shared helper functions for the 16S methanotroph analysis scripts (4.3-4.5).
#
# These mirror step_6/step_7's kaiju versions but operate on reconciled
# (canonical) genus names, since 16S/SILVA uses different taxonomy strings
# than kaiju/NCBI for a couple of these taxa (most importantly the anaerobic
# ANME-2d / NC10 genera). Reconciliation itself (reconcile_phyloseq(),
# fix_lineage_from_genus(), taxonomy_map) comes from
# step_8_comparing_16s_metagenomics/reconcile_phyloseq.R.

library(dplyr)
library(tidyr)
library(phyloseq)

# Same as step_6's taxa_of_interest, translated to canonical names via
# reconcile_phyloseq.R's taxonomy_map. Only the two anaerobic genera change
# ("Candidatus X" -> "Ca. X"); everything else (including step_6's existing
# "methylogea"/"Methylovolum"/"Methyloholbius"/"Methylocaspa" typos, kept
# as-is for consistency with the kaiju-side list) is unchanged.
taxa_of_interest_16s <- c(
  "Methylococcus", "Methylomonas", "Methylobacter", "Methylomicrobium",
  "Methylosarcina", "Methylocaldum", "methylogea", "Methylosoma",
  "Methyloparacoccus", "Methyloglobulus", "Methyloprofundus", "Methylomarinum",
  "Methylovolum", "Methylomagnum", "Methylosphaera", "Methylothermus",
  "Methyloholbius", "Methylomarinovum", "Methylosinus", "Methylocystis",
  "Methylocella", "Methylocaspa", "Methyloferula",
  "Ca. Methanoperedens", "Ca. Methylomirabilis",
  "Methylacidiphilum", "Methylacidimicrobium"
)

genus_order_16s <- taxa_of_interest_16s

# Same as step_6's compute_relative_abundance(), but the Metabolism check
# uses the canonical anaerobic genus names ("Ca. X") instead of kaiju's
# "Candidatus X" -- physeq must already be reconciled (reconcile_phyloseq())
# before calling this.
compute_relative_abundance_16s <- function(physeq, taxa_of_interest) {
  taxonomy_df <- as.data.frame(tax_table(physeq))
  otu_ids <- rownames(taxonomy_df[taxonomy_df$genus %in% taxa_of_interest, ])
  otu_to_genus <- setNames(taxonomy_df$genus[match(otu_ids, rownames(taxonomy_df))], otu_ids)

  physeq_subset <- prune_taxa(otu_ids, physeq)
  if (nsamples(physeq_subset) == 0 || ntaxa(physeq_subset) == 0) {
    stop("compute_relative_abundance_16s(): subset has zero samples or taxa")
  }

  counts_df <- as.data.frame(otu_table(physeq_subset))
  counts_df$genus <- otu_to_genus[rownames(counts_df)]

  counts_df_summed <- counts_df %>%
    group_by(genus) %>%
    summarise(across(everything(), sum), .groups = "drop")

  total_counts <- colSums(as.data.frame(otu_table(physeq)))

  counts_df_summed %>%
    mutate(across(-genus, ~ . / total_counts[cur_column()])) %>%
    pivot_longer(cols = -genus, names_to = "Sample", values_to = "RelativeAbundance") %>%
    mutate(
      Location   = ifelse(grepl("^E3", Sample), "Jetty", "CCS"),
      Year       = ifelse(grepl("^E5|^E6", Sample), "2023", "2022"),
      Metabolism = ifelse(genus %in% c("Ca. Methanoperedens", "Ca. Methylomirabilis"),
                           "Anaerobic", "Aerobic")
    )
}

# Loads the cached 16S phyloseq, reconciles taxonomy to canonical names,
# computes relative abundance for taxa_of_interest_16s, and orders
# Sample/genus factors consistently. Returns list(data, sample_order).
load_relative_abundance_long_16s <- function(phyloseq_path, taxa_of_interest, sample_list, genus_order) {
  physeq <- readRDS(phyloseq_path)
  physeq <- reconcile_phyloseq(physeq)
  physeq <- fix_lineage_from_genus(physeq)

  rel_ab <- compute_relative_abundance_16s(physeq, taxa_of_interest) %>%
    mutate(
      Sample = factor(Sample, levels = sample_list),
      genus  = factor(genus, levels = genus_order)
    )

  sample_order <- rel_ab %>% distinct(Sample) %>% arrange(Sample) %>% pull(Sample)
  rel_ab$Sample <- factor(rel_ab$Sample, levels = sample_order)

  list(data = rel_ab, sample_order = sample_order)
}
