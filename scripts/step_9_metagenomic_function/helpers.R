# Shared helper functions for the cycdb functional gene scripts (9.1-9.3).

library(dplyr)
library(readr)
library(stringr)
library(tibble)

# DIAMOND output filenames already carry the display sample name (they were
# renamed the same way as kaiju_output), e.g.
# "lane1-s001-...-1-E1R1_0_3cm_-CH4_+O2_processed_diamond.tsv" -> "E1R1_0_3cm_-CH4_+O2"
extract_sample_from_diamond_filename <- function(path) {
  raw_name <- tools::file_path_sans_ext(basename(path))
  # PCRE lookbehind must be fixed-length, so capture the sample name with a
  # group instead of a variable-length lookbehind.
  m <- regmatches(raw_name, regexec(
    "^lane1-s[0-9]{3}-indexN[0-9]+-S[0-9]+-[ACGT]+-[ACGT]+-[0-9]+-(.+)_processed_diamond$",
    raw_name, perl = TRUE
  ))[[1]]
  if (length(m) < 2) NA_character_ else m[2]
}

# Which pathway (list name in `lop`) a gene symbol belongs to, or NA.
check_string <- function(gene, lop) {
  for (list_name in names(lop)) {
    if (gene %in% lop[[list_name]]) return(list_name)
  }
  NA_character_
}

# Loads the cached pathway count matrix + experiment metadata, aligns them,
# and returns list(mat, coldata) ready for DESeq2. Shared by the deseq and
# heatmap scripts so the alignment/factor-coercion logic isn't duplicated.
load_pathway_matrix_and_coldata <- function(pathway_counts_path, exp_metadata_path, sample_id_map) {
  agg_counts_path_summary <- readRDS(pathway_counts_path)

  mat <- agg_counts_path_summary %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  mode(mat) <- "integer"

  # ant_exp_map_MA.csv is saved with a UTF-8 BOM (same issue as
  # metadata_ac_new.csv) -- rename the first column positionally.
  meta_raw <- read_csv(exp_metadata_path, show_col_types = FALSE, na = c("", "NA"))
  names(meta_raw)[1] <- "NAME"

  # sample_id_map is keyed on raw IDs with a trailing "_S##"; strip that to
  # match ant_exp_map_MA.csv's NAME column (some rows lack the suffix too).
  stripped_map <- setNames(sample_id_map, str_remove(names(sample_id_map), "_S\\d+$"))
  raw_sample <- str_remove(meta_raw$NAME, "_S\\d+$")
  meta_raw$sample <- ifelse(raw_sample %in% names(stripped_map), stripped_map[raw_sample], raw_sample)

  meta <- meta_raw %>% select(sample, everything(), -NAME)

  missing_in_meta <- setdiff(colnames(mat), meta$sample)
  missing_in_mat  <- setdiff(meta$sample, colnames(mat))
  if (length(missing_in_meta) > 0)
    warning("In count matrix but NOT metadata:\n  ", paste(missing_in_meta, collapse = "\n  "))
  if (length(missing_in_mat) > 0)
    warning("In metadata but NOT count matrix:\n  ", paste(missing_in_mat, collapse = "\n  "))

  shared <- intersect(colnames(mat), meta$sample)
  mat <- mat[, shared]

  coldata <- meta %>%
    filter(sample %in% shared) %>%
    arrange(match(sample, shared)) %>%
    column_to_rownames("sample") %>%
    mutate(
      Year      = factor(Year),
      Methane   = factor(Methane,   levels = c("No", "Yes")),
      Oxygen    = factor(Oxygen,    levels = c("No", "Yes")),
      Exp       = factor(Exp),
      depth     = factor(depth,     levels = c("0_3", "0_4", "3_6", "4_8", "6_9")),
      sediment  = factor(sediment,  levels = c("surface", "deep")),
      Site      = factor(Site),
      Oxidation = factor(Oxidation, levels = c("No", "Yes"))
    )

  stopifnot(all(rownames(coldata) == colnames(mat)))

  list(mat = mat, coldata = coldata)
}

# Writes a DESeq2 results() object to <out_dir>/deseq2_<label>.csv, sorted by padj.
write_deseq_results <- function(res, label, out_dir) {
  as.data.frame(res) %>%
    rownames_to_column("pathway") %>%
    arrange(padj) %>%
    write_csv(file.path(out_dir, paste0("deseq2_", label, ".csv")))
}

# Tidies a DESeq2 results() object into a plain data frame with significance flags.
tidy_deseq_results <- function(res, label) {
  as.data.frame(res) %>%
    rownames_to_column("pathway") %>%
    mutate(
      effect = label,
      sig    = !is.na(padj) & padj < 0.05,
      direction = case_when(
        sig & log2FoldChange > 0 ~ "Up",
        sig & log2FoldChange < 0 ~ "Down",
        TRUE ~ "n.s."
      )
    )
}

# Lollipop plot of significant (padj < padj_cutoff) pathways for one DESeq2 contrast.
lollipop_plot <- function(df, title, padj_cutoff = 0.05) {
  df_sig <- df %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    arrange(log2FoldChange) %>%
    mutate(pathway = factor(pathway, levels = pathway))

  if (nrow(df_sig) == 0) {
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(title, "(no significant pathways)")) + ggplot2::theme_void())
  }

  ggplot2::ggplot(df_sig, ggplot2::aes(x = log2FoldChange, y = pathway, color = log2FoldChange > 0)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = log2FoldChange, yend = pathway), linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40") +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
                                 labels = c(`TRUE` = "Higher", `FALSE` = "Lower"),
                                 name = NULL) +
    ggplot2::labs(title = title, x = "log2 fold change", y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
