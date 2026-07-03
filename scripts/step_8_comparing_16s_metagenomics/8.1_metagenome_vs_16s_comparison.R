# Compares 16S rRNA amplicon and shotgun metagenomic (kaiju) taxonomy for
# methanotroph/methanogen-relevant taxa: distinct-taxa counts by domain,
# per-taxon detection across both methods, and per-sample read/relative-
# abundance supplementary tables.
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed. Requires 6.1_load_kaiju_data.R (step_6) to have been run
# at least once, since the metagenomic side reuses its cached phyloseq
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
    "reconcile_phyloseq.R, step_6's constants.R, and the project's ",
    "updated_16s_run/kaiju_output/figures/supplementary folders). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
taxonomy_dir <- file.path(project_root, "scripts", "step_6_metagenomic_taxonomy")
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(phyloseq)

# sample_id_map / kaiju_phyloseq_path / figures_dir
source(file.path(taxonomy_dir, "constants.R"))
# reconcile_phyloseq() / fix_lineage_from_genus() / taxonomy_map / etc.
source(file.path(script_dir, "reconcile_phyloseq.R"))

run_16s_dir        <- file.path(project_root, "updated_16s_run")
asv_16s_path       <- file.path(run_16s_dir, "microbiome_analyst_16s_update.csv")
metadata_16s_path  <- file.path(run_16s_dir, "ant_exp_map_MA.csv")
supplementary_dir  <- file.path(project_root, "supplementary")
dir.create(supplementary_dir, showWarnings = FALSE)

if (!file.exists(kaiju_phyloseq_path)) {
  stop(
    "Missing ", kaiju_phyloseq_path, ". Run scripts/step_6_metagenomic_taxonomy/",
    "6.1_load_kaiju_data.R first to build the cached kaiju phyloseq this ",
    "script reuses for the metagenomic side of the comparison."
  )
}

# =============================================================================
# 16S PROCESSING
# =============================================================================
# The ASV table's first column holds each row's taxonomy string (not a
# sample name, despite the "#NAME" header) -- rows are ASVs, columns are
# per-sample counts.
data_16s <- read.csv(asv_16s_path, check.names = FALSE)
colnames(data_16s)[1] <- "taxon_string"

counts_16s   <- select(data_16s, -taxon_string)
taxonomy_16s <- select(data_16s, taxon_string)
colnames(taxonomy_16s) <- "taxon_name"

max_splits_16s <- max(str_count(taxonomy_16s$taxon_name, ";")) + 1
taxonomy_16s$taxon_name <- gsub(";$", "", taxonomy_16s$taxon_name)

rank_names_16s <- c("domain", "phylum", "class", "order", "family", "genus", "species")
if (max_splits_16s > 7) {
  rank_names_16s <- c(rank_names_16s, paste0("extra_", 1:(max_splits_16s - 7)))
}

clean_silva_prefixes <- function(taxonomy) {
  gsub("[a-z]__", "", taxonomy)
}

taxonomy_16s <- sapply(taxonomy_16s, clean_silva_prefixes) %>% as.data.frame()

taxonomy_16s_split <- taxonomy_16s %>%
  separate(taxon_name, into = rank_names_16s, sep = ";", fill = "right") %>%
  mutate_all(~ replace_na(., "Unclassified")) %>%
  mutate_all(~ gsub("__", "Unclassified", .)) %>%
  mutate(asv_id = 1:n())

counts_16s <- counts_16s %>% mutate(asv_id = 1:n())

otu_mat_16s <- counts_16s         %>% column_to_rownames("asv_id") %>% as.matrix()
tax_mat_16s <- taxonomy_16s_split %>% column_to_rownames("asv_id") %>% as.matrix()

# ant_exp_map_MA.csv is saved with a UTF-8 BOM (same issue fixed in
# step_6/5.1 and step_9) -- rename the first column positionally.
metadata_16s <- read.csv(metadata_16s_path, check.names = FALSE)
names(metadata_16s)[1] <- "NAME"
metadata_16s <- metadata_16s %>% filter(NAME != "Jan24_B", NAME != "Jan24_B2")
samples_df_16s <- metadata_16s %>% column_to_rownames("NAME")

physeq_16s <- phyloseq(
  otu_table(otu_mat_16s, taxa_are_rows = TRUE),
  tax_table(tax_mat_16s),
  sample_data(samples_df_16s)
)

# Apply the sample naming function (same sample_id_map as step_6/6/8), with
# the trailing "_S##" stripped since a few raw IDs lack it inconsistently.
stripped_map_16s <- setNames(sample_id_map, str_remove(names(sample_id_map), "_S\\d+$"))
raw_sample_16s    <- str_remove(sample_names(physeq_16s), "_S\\d+$")
new_sample_names_16s <- unname(stripped_map_16s[raw_sample_16s])
if (anyNA(new_sample_names_16s)) {
  stop(
    "sample_id_map is missing an entry for: ",
    paste(sample_names(physeq_16s)[is.na(new_sample_names_16s)], collapse = ", ")
  )
}
sample_names(physeq_16s) <- new_sample_names_16s

# =============================================================================
# METAGENOMIC PROCESSING (Kaiju short reads)
# =============================================================================
# Reuses the phyloseq already built (and sample-renamed) by
# step_6/6.1_load_kaiju_data.R rather than re-parsing the 34 raw tsvs here.
physeq_mg <- readRDS(kaiju_phyloseq_path)

# =============================================================================
# RECONCILE TAXONOMY
# =============================================================================
physeq_16s <- reconcile_phyloseq(physeq_16s)
physeq_mg  <- reconcile_phyloseq(physeq_mg)

# Fill in missing order/family from genus (especially NC10 in Kaiju)
physeq_16s <- fix_lineage_from_genus(physeq_16s)
physeq_mg  <- fix_lineage_from_genus(physeq_mg)

# =============================================================================
# Distinct taxa counts (Archaea vs Bacteria) -- order, family, genus
# =============================================================================
count_distinct_taxa <- function(physeq, rank, method_label,
                                mg_min_relabund = 0.00001,  # 0.001%
                                s16_min_reads = 2) {
  tax_df <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)
  otu_df <- as.data.frame(otu_table(physeq))
  top_rank_col <- intersect(c("domain", "superkingdom"), colnames(tax_df))[1]

  if (!rank %in% colnames(tax_df)) {
    return(tibble(domain = character(0), n_distinct = integer(0),
                  rank = rank, method = method_label))
  }

  is_mg <- grepl("Metagenomics|Kaiju", method_label, ignore.case = TRUE)

  # Per-sample relative abundances (denominator includes unclassified reads)
  lib_sizes <- colSums(otu_df)
  relabund_df <- sweep(otu_df, 2, lib_sizes, "/")

  tax_df$.taxon_at_rank <- tax_df[[rank]]
  tax_df$.domain        <- tax_df[[top_rank_col]]

  # Keep only Archaea/Bacteria rows with a classified taxon at this rank
  keep_rows <- tax_df$.domain %in% c("Archaea", "Bacteria") &
    tax_df$.taxon_at_rank != "Unclassified"

  tax_df       <- tax_df[keep_rows, , drop = FALSE]
  otu_df       <- otu_df[keep_rows, , drop = FALSE]
  relabund_df  <- relabund_df[keep_rows, , drop = FALSE]

  # For each unique (domain, taxon) pair, sum across matching rows and check
  # whether the taxon passes the threshold in any sample
  unique_taxa <- unique(tax_df[, c(".domain", ".taxon_at_rank")])

  passes <- mapply(function(dom, tx) {
    idx <- which(tax_df$.domain == dom & tax_df$.taxon_at_rank == tx)
    if (length(idx) == 1) {
      reads_per_sample    <- as.numeric(otu_df[idx, ])
      relabund_per_sample <- as.numeric(relabund_df[idx, ])
    } else {
      reads_per_sample    <- as.numeric(colSums(otu_df[idx, , drop = FALSE]))
      relabund_per_sample <- as.numeric(colSums(relabund_df[idx, , drop = FALSE]))
    }

    if (is_mg) {
      any(relabund_per_sample >= mg_min_relabund)
    } else {
      any(reads_per_sample >= s16_min_reads)
    }
  }, unique_taxa$.domain, unique_taxa$.taxon_at_rank)

  unique_taxa$passes <- passes

  unique_taxa %>%
    filter(passes) %>%
    group_by(domain = .domain) %>%
    summarise(n_distinct = n(), .groups = "drop") %>%
    mutate(rank = rank, method = method_label)
}

# =============================================================================
# Methanotroph detection -- order, family, genus
# All names use CANONICAL forms from the reconciliation table
# =============================================================================

# ---- Orders ----
methanotroph_orders <- c(
  # Aerobic -- Type I (Gammaproteobacteria)
  "Methylococcales",
  # Aerobic -- Type II (Alphaproteobacteria)
  "Rhizobiales/Hyphomicrobiales",
  # Aerobic -- Verrucomicrobia
  "Methylacidiphilales",
  # Anaerobic -- ANME-2/3
  "Methanosarcinales",
  # Anaerobic -- ANME-1
  "Ca. Methanophagales",
  # Anaerobic -- NC10
  "Methylomirabilales"
)

# ---- Families ----
methanotroph_families <- c(
  # Aerobic -- Type I
  "Methylococcaceae/Methylomonadaceae",
  "Methylothermaceae",
  # Aerobic -- Type II
  "Methylocystaceae",
  "Beijerinckiaceae",
  # Aerobic -- Verrucomicrobia
  "Methylacidiphilaceae",
  # Anaerobic -- ANME
  "Methanoperedenaceae",
  "Methanocomedenaceae",
  "Methanogasteraceae",
  "Methanosarcinaceae",
  "Methanophagaceae",
  # NC10
  "Methylomirabilaceae"
)

# ---- Genera ----
methanotroph_genera <- c(
  # --- Aerobic Type I (Methylococcaceae/Methylomonadaceae) ---
  "Methylococcus", "Methylomonas", "Methylobacter", "Methylomicrobium",
  "Methylosarcina", "Methylocaldum", "Methylogaea", "Methylosoma",
  "Methyloparacoccus", "Methyloglobulus", "Methyloprofundus", "Methylomarinum",
  "Methylovulum", "Methylomagnum", "Methylosphaera",
  # --- Aerobic Type I (Methylothermaceae) ---
  "Methylothermus", "Methylohalobius", "Methylomarinovum",
  # --- Aerobic Type II (Methylocystaceae) ---
  "Methylosinus", "Methylocystis",
  # --- Aerobic Type II (Beijerinckiaceae) ---
  "Methylocella", "Methylocapsa", "Methyloferula",
  # --- Aerobic Verrucomicrobia ---
  "Methylacidiphilum", "Methylacidimicrobium",
  # --- Anaerobic ANME-1 ---
  "Methanophaga", "Methanoalium",
  # --- Anaerobic ANME-2a ---
  "Methanocomedens",
  # --- Anaerobic ANME-2b ---
  "Methanomarinus",
  # --- Anaerobic ANME-2c ---
  "Methanogaster",
  # --- Anaerobic ANME-2d ---
  "Ca. Methanoperedens",
  # --- Anaerobic ANME-3 ---
  "Methanovorans",
  # --- NC10 ---
  "Ca. Methylomirabilis"
)

# =============================================================================
# Detection functions -- with filtering thresholds
# =============================================================================
# Metagenomics: taxon must have > 500 reads OR > 0.01% relative abundance
#               in at least 4 samples
# 16S:          taxon must have > 4 reads across at least 4 samples

methanotroph_detected <- function(physeq, target_taxa, rank, method_label,
                                  mg_min_relabund = 0.00001,  # 0.001%
                                  s16_min_reads = 2) {
  tax_df <- as.data.frame(tax_table(physeq))
  otu_df <- as.data.frame(otu_table(physeq))

  if (!rank %in% colnames(tax_df)) {
    return(tibble(taxon = target_taxa, detected = FALSE,
                  method = method_label, rank = rank))
  }

  # Library sizes include ALL reads (classified + unclassified)
  lib_sizes <- colSums(otu_df)
  relabund_df <- sweep(otu_df, 2, lib_sizes, "/")

  is_mg <- grepl("Metagenomics|Kaiju", method_label, ignore.case = TRUE)

  detected_vec <- sapply(target_taxa, function(t) {
    hits <- which(tolower(tax_df[[rank]]) == tolower(t))
    if (length(hits) == 0) return(FALSE)

    if (length(hits) == 1) {
      reads_per_sample    <- otu_df[hits, ]
      relabund_per_sample <- relabund_df[hits, ]
    } else {
      reads_per_sample    <- colSums(otu_df[hits, , drop = FALSE])
      relabund_per_sample <- colSums(relabund_df[hits, , drop = FALSE])
    }

    if (is_mg) {
      # Metagenomics: detected if >= 0.001% relative abundance in any sample
      any(relabund_per_sample >= mg_min_relabund)
    } else {
      # 16S: detected if >= 2 reads (non-singleton) in any sample
      any(reads_per_sample >= s16_min_reads)
    }
  })

  tibble(taxon = target_taxa, detected = detected_vec,
         method = method_label, rank = rank)
}

build_detection <- function(target_taxa, rank) {
  bind_rows(
    methanotroph_detected(physeq_16s, target_taxa, rank, "16S"),
    methanotroph_detected(physeq_mg,  target_taxa, rank, "Metagenomics (Kaiju)")
  )
}

det_order  <- build_detection(methanotroph_orders,   "order")
det_family <- build_detection(methanotroph_families, "family")
det_genus  <- build_detection(methanotroph_genera,   "genus")

ranks_to_compare <- c("order", "family", "genus")

counts_16s_summary <- bind_rows(lapply(ranks_to_compare, function(r) {
  count_distinct_taxa(physeq_16s, r, "16S")
}))
counts_mg_summary <- bind_rows(lapply(ranks_to_compare, function(r) {
  count_distinct_taxa(physeq_mg, r, "Metagenomics (Kaiju)")
}))

taxa_counts_combined <- bind_rows(counts_16s_summary, counts_mg_summary) %>%
  mutate(
    rank   = factor(rank,   levels = c("order", "family", "genus"),
                    labels = c("Order", "Family", "Genus")),
    method = factor(method, levels = c("16S", "Metagenomics (Kaiju)"),
                    labels = c("16S rRNA", "Metagenomics"))
  )

# =============================================================================
# Publication theme & plot functions
# =============================================================================

theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid        = element_blank(),
      panel.background  = element_blank(),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.text         = element_text(color = "black", size = base_size),
      axis.ticks        = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_blank(),
      plot.title        = element_text(size = base_size + 1, face = "bold",
                                       hjust = 0, margin = margin(b = 6)),
      plot.title.position = "plot",
      legend.position   = "top",
      legend.justification = "left",
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_size - 1),
      legend.key.size   = unit(10, "pt"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, -4, 0),
      plot.margin       = margin(4, 8, 4, 4)
    )
}

# =============================================================================
# FACETED PLOTS -- one bar plot, one heatmap
# =============================================================================

# ---- Faceted bar plot: distinct taxa by rank ----
p_bars <- ggplot(taxa_counts_combined,
                 aes(x = domain, y = n_distinct, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_distinct),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.5) +
  facet_wrap(~ rank, scales = "free_y") +
  scale_fill_manual(values = c("16S rRNA"     = "#2c7fb8",
                               "Metagenomics" = "#d95f0e")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(y = "Number of distinct taxa", x = NULL, fill = NULL) +
  theme_pub() +
  theme(
    axis.text.x    = element_text(face = "bold", size = 12),
    axis.text.y    = element_text(face = "plain", size = 10),
    axis.title.y   = element_text(size = 11),
    strip.text     = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "white", color = "black",
                                    linewidth = 0.3),
    legend.position = "top",
    legend.justification = "left",
    legend.text    = element_text(size = 11),
    legend.key.size = unit(14, "pt")
  )

ggsave(file.path(figures_dir, "Figure7_ABC_barplot_distinct_taxa.png"), p_bars,
       width = 12, height = 5, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "Figure7_ABC_barplot_distinct_taxa.pdf"), p_bars,
       width = 12, height = 5, dpi = 600, bg = "white")

# ---- Faceted heatmap: detection across order / family / genus ----

# Define aerobic vs anaerobic classification for each taxon
metabolism_map <- tribble(
  ~taxon, ~metabolism,
  # Orders
  "Methylococcales",                       "Aerobic",
  "Rhizobiales/Hyphomicrobiales",          "Aerobic",
  "Methylacidiphilales",                   "Aerobic",
  "Methanosarcinales",                     "Anaerobic",
  "Ca. Methanophagales",                   "Anaerobic",
  "Methylomirabilales",                    "Anaerobic",
  # Families
  "Methylococcaceae/Methylomonadaceae",    "Aerobic",
  "Methylothermaceae",                     "Aerobic",
  "Methylocystaceae",                      "Aerobic",
  "Beijerinckiaceae",                      "Aerobic",
  "Methylacidiphilaceae",                  "Aerobic",
  "Methanoperedenaceae",                   "Anaerobic",
  "Methanocomedenaceae",                   "Anaerobic",
  "Methanogasteraceae",                    "Anaerobic",
  "Methanosarcinaceae",                    "Anaerobic",
  "Methanophagaceae",                      "Anaerobic",
  "Methylomirabilaceae",                   "Anaerobic",
  # Genera -- aerobic Type I
  "Methylococcus", "Aerobic", "Methylomonas", "Aerobic",
  "Methylobacter", "Aerobic", "Methylomicrobium", "Aerobic",
  "Methylosarcina", "Aerobic", "Methylocaldum", "Aerobic",
  "Methylogaea", "Aerobic", "Methylosoma", "Aerobic",
  "Methyloparacoccus", "Aerobic", "Methyloglobulus", "Aerobic",
  "Methyloprofundus", "Aerobic", "Methylomarinum", "Aerobic",
  "Methylovulum", "Aerobic", "Methylomagnum", "Aerobic",
  "Methylosphaera", "Aerobic",
  # Genera -- aerobic Methylothermaceae
  "Methylothermus", "Aerobic", "Methylohalobius", "Aerobic",
  "Methylomarinovum", "Aerobic",
  # Genera -- aerobic Type II
  "Methylosinus", "Aerobic", "Methylocystis", "Aerobic",
  "Methylocella", "Aerobic", "Methylocapsa", "Aerobic",
  "Methyloferula", "Aerobic",
  # Genera -- aerobic Verrucomicrobia
  "Methylacidiphilum", "Aerobic", "Methylacidimicrobium", "Aerobic",
  # Genera -- anaerobic ANME
  "Methanophaga", "Anaerobic", "Methanoalium", "Anaerobic",
  "Methanocomedens", "Anaerobic", "Methanomarinus", "Anaerobic",
  "Methanogaster", "Anaerobic",
  "Ca. Methanoperedens", "Anaerobic",
  "Methanovorans", "Anaerobic",
  # Genera -- NC10
  "Ca. Methylomirabilis", "Anaerobic"
)

# Combine all detection data
det_all <- bind_rows(
  det_order  %>% mutate(rank = "Order"),
  det_family %>% mutate(rank = "Family"),
  det_genus  %>% mutate(rank = "Genus")
) %>%
  left_join(metabolism_map, by = "taxon") %>%
  mutate(
    rank   = factor(rank, levels = c("Order", "Family", "Genus")),
    method = factor(method,
                    levels = c("Metagenomics (Kaiju)", "16S"),
                    labels = c("Metagenomics", "16S rRNA"))
  )

# Build sort order: anaerobic first, then aerobic, alphabetical within each
taxon_levels <- det_all %>%
  select(taxon, rank, metabolism) %>%
  distinct() %>%
  mutate(metabolism = factor(metabolism, levels = c("Anaerobic", "Aerobic"))) %>%
  arrange(rank, metabolism, taxon) %>%
  pull(taxon) %>%
  unique()

det_all <- det_all %>%
  mutate(taxon = factor(taxon, levels = taxon_levels))

# Use numeric y-axis so we can control tile heights precisely:
#   y = 2 -> Metagenomics (height 1)
#   y = 1 -> 16S rRNA (height 1)
#   y = 0.15 -> annotation strip (height 0.3)
# Then relabel with scale_y_continuous breaks/labels

det_numeric <- det_all %>%
  mutate(y_pos = ifelse(method == "Metagenomics", 2, 1))

annot_data <- det_all %>%
  select(taxon, rank, metabolism) %>%
  distinct()

# Create a combined fill column so everything goes through one scale
det_plot <- det_numeric %>%
  mutate(fill_val = ifelse(detected, "Detected", "Not detected"))

annot_plot <- annot_data %>%
  mutate(fill_val = metabolism,
         y_pos = 0.15)

p_heatmap <- ggplot() +
  # Main detection tiles
  geom_tile(data = det_plot,
            aes(x = taxon, y = y_pos, fill = fill_val),
            color = "white", linewidth = 0.8,
            height = 0.9) +
  # Thin annotation strip
  geom_tile(data = annot_plot,
            aes(x = taxon, y = y_pos, fill = fill_val),
            color = "white", linewidth = 0.4,
            height = 0.25) +
  facet_wrap(~ rank, scales = "free_x",
             ncol = 1,
             strip.position = "left") +
  scale_fill_manual(
    values = c("Detected"     = "#1a5276",
               "Not detected" = "#f0f0f0",
               "Aerobic"      = "#e67e22",
               "Anaerobic"    = "#8e44ad"),
    breaks = c("Detected", "Not detected", "Aerobic", "Anaerobic"),
    name = NULL
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("16S rRNA", "Metagenomics"),
    expand = expansion(add = c(0.3, 0.3))
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_blank(),
    panel.spacing.y  = unit(1.2, "lines"),

    axis.text.x      = element_text(angle = 50, hjust = 1, vjust = 1,
                                    face = "italic", size = 9,
                                    color = "grey10"),
    axis.text.y      = element_text(face = "bold", size = 10,
                                    color = "grey10", margin = margin(r = 4)),
    axis.ticks       = element_blank(),

    strip.text.y.left = element_text(face = "bold", size = 12, angle = 0,
                                     margin = margin(r = 8)),
    strip.placement  = "outside",
    strip.background = element_blank(),

    legend.position      = "top",
    legend.justification = "left",
    legend.text          = element_text(size = 10, color = "grey10"),
    legend.key.size      = unit(14, "pt"),
    legend.key           = element_rect(color = "grey50", linewidth = 0.3),
    legend.margin        = margin(b = 8),

    plot.margin = margin(t = 6, r = 10, b = 6, l = 6)
  ) +
  guides(fill = guide_legend(nrow = 1))

ggsave(file.path(figures_dir, "Figure7_DEF_heatmap_detection.png"), p_heatmap,
       width = 12, height = 11, dpi = 600, bg = "white")
ggsave(file.path(figures_dir, "Figure7_DEF_heatmap_detection.pdf"), p_heatmap,
       width = 12, height = 11, dpi = 600, bg = "white")

message("Done! Two faceted plots saved to ", figures_dir, ":")
message("  barplot_distinct_taxa.png")
message("  heatmap_detection.png")

# =============================================================================
# Supplementary table: read counts per methanotroph genus per sample
# =============================================================================

# ---- Helper function to extract counts for target taxa at a given rank ----
extract_methanotroph_counts <- function(physeq, target_taxa, rank, method_label) {
  tax_df <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)
  otu_df <- as.data.frame(otu_table(physeq))

  if (!rank %in% colnames(tax_df)) {
    message("Rank '", rank, "' not found in ", method_label)
    return(NULL)
  }

  lib_sizes <- colSums(otu_df)

  results <- lapply(target_taxa, function(t) {
    hits <- which(tolower(tax_df[[rank]]) == tolower(t))
    if (length(hits) == 0) return(NULL)

    if (length(hits) == 1) {
      reads_per_sample <- as.numeric(otu_df[hits, ])
    } else {
      reads_per_sample <- as.numeric(colSums(otu_df[hits, , drop = FALSE]))
    }

    relabund_per_sample <- reads_per_sample / lib_sizes * 100

    tibble(
      taxon         = t,
      sample        = colnames(otu_df),
      reads         = reads_per_sample,
      rel_abund_pct = round(relabund_per_sample, 4)
    )
  })

  bind_rows(results) %>%
    mutate(method = method_label)
}

# ---- Extract counts for both methods at genus level ----
counts_16s_genera <- extract_methanotroph_counts(
  physeq_16s, methanotroph_genera, "genus", "16S")
counts_mg_genera <- extract_methanotroph_counts(
  physeq_mg, methanotroph_genera, "genus", "Metagenomics")

supp_table_long <- bind_rows(counts_16s_genera, counts_mg_genera)

# ---- Summary table: total reads, mean reads, n_samples detected, max relabund ----
supp_table_summary <- supp_table_long %>%
  group_by(method, taxon) %>%
  summarise(
    total_reads       = sum(reads),
    mean_reads        = round(mean(reads), 1),
    max_reads         = max(reads),
    n_samples_present = sum(reads > 0),
    mean_rel_abund    = round(mean(rel_abund_pct), 4),
    max_rel_abund     = round(max(rel_abund_pct), 4),
    .groups = "drop"
  ) %>%
  arrange(method, desc(total_reads))

print(supp_table_summary, n = 100)

write.csv(supp_table_summary,
          file.path(supplementary_dir, "TableS6_supplementary_methanotroph_genus_summary.csv"),
          row.names = FALSE)
write.csv(supp_table_long %>% pivot_wider(
  names_from = sample,
  values_from = c(reads, rel_abund_pct),
  values_fill = 0),
  file.path(supplementary_dir, "TableS7_supplementary_methanotroph_genus_per_sample.csv"),
  row.names = FALSE)

message("Supplementary tables saved to ", supplementary_dir, ":")
message("  supplementary_methanotroph_genus_summary.csv")
message("  supplementary_methanotroph_genus_per_sample.csv")
