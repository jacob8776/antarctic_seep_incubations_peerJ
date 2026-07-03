# Shared helper functions for the kaiju taxonomy/visualization scripts
# (step_6_metagenomic_taxonomy and step_7_methanotroph_visualization).

library(dplyr)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Relative abundance (per sample, vs. total classified reads) of taxa_of_interest,
# summed to genus and reshaped long, with Site/Year/Location/Metabolism attached.
compute_relative_abundance <- function(physeq, taxa_of_interest) {
  taxonomy_df <- as.data.frame(tax_table(physeq))
  otu_ids <- rownames(taxonomy_df[taxonomy_df$genus %in% taxa_of_interest, ])
  otu_to_genus <- setNames(taxonomy_df$genus[match(otu_ids, rownames(taxonomy_df))], otu_ids)

  physeq_subset <- prune_taxa(otu_ids, physeq)
  if (nsamples(physeq_subset) == 0 || ntaxa(physeq_subset) == 0) {
    stop("compute_relative_abundance(): subset has zero samples or taxa")
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
      Metabolism = ifelse(genus %in% c("Candidatus Methanoperedens", "Candidatus Methylomirabilis"),
                           "Anaerobic", "Aerobic")
    )
}

# Loads the cached phyloseq object, computes relative abundance for
# taxa_of_interest, and orders Sample/genus factors consistently.
# Returns list(data, sample_order).
load_relative_abundance_long <- function(phyloseq_path, taxa_of_interest, sample_list, genus_order) {
  physeq <- readRDS(phyloseq_path)

  rel_ab <- compute_relative_abundance(physeq, taxa_of_interest) %>%
    mutate(
      Sample = factor(Sample, levels = sample_list),
      genus  = factor(genus, levels = genus_order)
    )

  sample_order <- rel_ab %>% distinct(Sample) %>% arrange(Sample) %>% pull(Sample)
  rel_ab$Sample <- factor(rel_ab$Sample, levels = sample_order)

  list(data = rel_ab, sample_order = sample_order)
}

# Top N taxa at a given taxonomic rank, expressed as relative abundance.
get_top_taxa_by_rank <- function(physeq, rank = "order", top_n = 15) {
  physeq_rel_abund <- transform_sample_counts(physeq, function(x) x / sum(x))
  physeq_rank <- tax_glom(physeq_rel_abund, taxrank = rank)
  top_taxa <- names(sort(taxa_sums(physeq_rank), decreasing = TRUE))[1:top_n]
  prune_taxa(top_taxa, physeq_rank)
}

extend_palette <- function(n) {
  colorRampPalette(brewer.pal(12, "Paired"))(n)
}

# Labels each genus with its phylum for the legend, e.g. "Woeseia (Proteobacteria)".
label_genus_with_phylum <- function(physeq) {
  tt <- as.data.frame(tax_table(physeq))
  tt$genus <- paste0(tt$genus, " (", tt$phylum, ")")
  tax_table(physeq) <- as.matrix(tt)
  physeq
}

# Stacked bar chart of `physeq` abundance at `rank`, with a location color
# bar on top. If highlighted_taxa is given, all other taxa are greyed out.
plot_top_taxa_bar <- function(physeq, rank, title, sample_order, highlighted_taxa = NULL, n_colors = 15) {
  physeq_melt <- psmelt(physeq)
  physeq_melt$Sample <- factor(physeq_melt$Sample, levels = sample_order)

  sample_sites <- physeq_melt %>%
    distinct(Sample, Exp) %>%
    mutate(site_label = ifelse(Exp == 1, "Jetty", "CCS"))

  taxa_levels <- levels(factor(physeq_melt[[rank]]))
  color_map <- setNames(extend_palette(n_colors)[seq_along(taxa_levels)], taxa_levels)

  if (!is.null(highlighted_taxa)) {
    color_map[!(taxa_levels %in% highlighted_taxa)] <- "grey80"
  }

  p_top <- ggplot(sample_sites, aes(x = Sample, y = 1, fill = site_label)) +
    geom_tile() +
    scale_fill_manual(values = c("CCS" = "#1B9E77", "Jetty" = "#F0C808"), name = "Location") +
    theme_void() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )

  p_main <- ggplot(physeq_melt, aes_string(x = "Sample", y = "Abundance", fill = rank)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = color_map) +
    labs(x = "Sample", y = "Relative Abundance") +
    guides(fill = guide_legend(title = rank)) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank()
    )

  (p_top / p_main) +
    plot_layout(heights = c(0.04, 1)) +
    plot_annotation(title = title, theme = theme(plot.title = element_text(size = 16, face = "bold")))
}
