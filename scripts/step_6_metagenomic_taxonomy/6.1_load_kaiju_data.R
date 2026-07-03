# Reads raw kaiju genus-level tsvs, builds the phyloseq object, applies
# updated sample names, and caches the result for scripts 5.2, 6.1, 6.2.
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed.

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
library(tidyr)
library(stringr)
library(tibble)
library(phyloseq)

source(file.path(script_dir, "constants.R"))

setwd(kaiju_dir)

tsv_column_names <- c("file", "percent", "reads", "taxon_id", "taxon_name")
files <- list.files(pattern = "*.tsv")

data_list <- lapply(files, function(f) {
  data <- read.csv(f, sep = "\t", col.names = tsv_column_names)
  data$taxon_name <- gsub("cellular organisms;", "", data$taxon_name)
  data
})
names(data_list) <- files

data_list <- lapply(data_list, function(df) filter(df, !is.na(taxon_id)))

# Recover the McM2X_... sample ID from the original kaiju.out path recorded
# in the "file" column.
extract_relevant_part <- function(df) {
  df$file <- str_extract(df$file, "McM.*(?=_kaiju\\.out)")
  df
}
data_list <- lapply(data_list, extract_relevant_part)

list_reads <- lapply(data_list, function(df) select(df, file, taxon_id, reads))
list_taxa  <- lapply(data_list, function(df) select(df, taxon_id, taxon_name))

combined_reads <- bind_rows(list_reads)

wide_reads <- combined_reads %>%
  pivot_wider(names_from = file, values_from = reads, values_fill = list(reads = 0))

combined_taxa <- list_taxa %>% bind_rows() %>% distinct()

# taxon_id must line up 1:1 and in the same order between the taxonomy and
# count tables before they can be bound into one phyloseq object.
stopifnot(identical(combined_taxa$taxon_id, wide_reads$taxon_id))

filter_taxonomy <- function(taxonomy_string) {
  terms <- strsplit(taxonomy_string, ";")[[1]]
  terms <- terms[!grepl("group$", terms) & !grepl("^unclassified", terms)]
  paste(terms, collapse = ";")
}

combined_taxa_filtered <- combined_taxa %>%
  mutate(taxon_name = sapply(taxon_name, filter_taxonomy)) %>%
  select(taxon_id, taxon_name) %>%
  mutate(taxon_name = gsub(";$", "", taxon_name))

taxonomy_columns <- c("superkingdom", "phylum", "class", "order", "family", "genus")

combined_taxa_split <- combined_taxa_filtered %>%
  separate(taxon_name, into = taxonomy_columns, sep = ";", fill = "right") %>%
  mutate_all(~ replace_na(., "Unclassified"))

otu_mat <- as.matrix(tibble::column_to_rownames(wide_reads, "taxon_id"))
tax_mat <- as.matrix(tibble::column_to_rownames(combined_taxa_split, "taxon_id"))

# metadata_ac_new.csv is saved with a UTF-8 BOM, which R's default parser
# folds into the first column's name (e.g. "X....NAME" instead of "NAME").
# Renaming positionally sidesteps however that mangling comes out.
metadata <- read.csv("metadata_ac_new.csv", check.names = FALSE)
names(metadata)[1] <- "NAME"
samples_df <- tibble::column_to_rownames(metadata, "NAME")

kaiju_phyloseq <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(samples_df)
)

# Apply the sample naming function: look up each raw sample ID and fail
# loudly rather than silently writing NA into a sample name.
new_sample_names <- unname(sample_id_map[sample_names(kaiju_phyloseq)])
if (anyNA(new_sample_names)) {
  stop(
    "sample_id_map is missing an entry for: ",
    paste(sample_names(kaiju_phyloseq)[is.na(new_sample_names)], collapse = ", ")
  )
}
sample_names(kaiju_phyloseq) <- new_sample_names

saveRDS(kaiju_phyloseq, kaiju_phyloseq_path)
kaiju_phyloseq
