# Reads the 16S ASV table + experiment metadata, builds the phyloseq object,
# applies the shared sample naming (step_6's sample_id_map), and caches
# full/Bacteria/Archaea phyloseq objects for scripts 4.2+.
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
    "constants.R, step_6's constants.R, and the project's updated_16s_run folder). ",
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
library(phyloseq)

source(file.path(taxonomy_dir, "constants.R"))  # sample_id_map
source(file.path(script_dir, "constants.R"))

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

# "domain" (not "superkingdom") for the top rank -- standard for SILVA/16S
# data and matches step_8's domain-vs-superkingdom auto-detection.
rank_names_16s <- c("domain", "phylum", "class", "order", "family", "genus", "species")
if (max_splits_16s > 7) {
  rank_names_16s <- c(rank_names_16s, paste0("extra_", 1:(max_splits_16s - 7)))
}

clean_silva_prefixes <- function(taxonomy) gsub("[a-z]__", "", taxonomy)
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
# step_6/step_8/step_9) -- rename the first column positionally.
metadata_16s <- read.csv(metadata_16s_path, check.names = FALSE)
names(metadata_16s)[1] <- "NAME"
metadata_16s <- metadata_16s %>% filter(NAME != "Jan24_B", NAME != "Jan24_B2")
samples_df_16s <- metadata_16s %>% column_to_rownames("NAME")

physeq_16s <- phyloseq(
  otu_table(otu_mat_16s, taxa_are_rows = TRUE),
  tax_table(tax_mat_16s),
  sample_data(samples_df_16s)
)

# Apply the sample naming function (same sample_id_map as steps 6-9), with
# the trailing "_S##" stripped since a few raw IDs lack it inconsistently.
stripped_map_16s <- setNames(sample_id_map, str_remove(names(sample_id_map), "_S\\d+$"))
raw_sample_16s <- str_remove(sample_names(physeq_16s), "_S\\d+$")
new_sample_names_16s <- unname(stripped_map_16s[raw_sample_16s])
if (anyNA(new_sample_names_16s)) {
  stop(
    "sample_id_map is missing an entry for: ",
    paste(sample_names(physeq_16s)[is.na(new_sample_names_16s)], collapse = ", ")
  )
}
sample_names(physeq_16s) <- new_sample_names_16s

physeq_16s_bac <- subset_taxa(physeq_16s, domain %in% c("Bacteria"))
physeq_16s_ar  <- subset_taxa(physeq_16s, domain %in% c("Archaea"))

saveRDS(physeq_16s,     phyloseq_16s_path)
saveRDS(physeq_16s_bac, phyloseq_16s_bac_path)
saveRDS(physeq_16s_ar,  phyloseq_16s_ar_path)

physeq_16s
