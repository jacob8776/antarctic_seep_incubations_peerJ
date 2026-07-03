# Shared constants for the 16S phyloseq scripts (4.1-4.2).
#
# Requires `project_root` to already be defined in the calling environment
# (each driver script sets it from its own location before sourcing this file).

run_16s_dir       <- file.path(project_root, "updated_16s_run")
asv_16s_path      <- file.path(run_16s_dir, "microbiome_analyst_16s_update.csv")
metadata_16s_path <- file.path(run_16s_dir, "ant_exp_map_MA.csv")

phyloseq_16s_path     <- file.path(run_16s_dir, "phyloseq_16s.rds")
phyloseq_16s_bac_path <- file.path(run_16s_dir, "phyloseq_16s_bac.rds")
phyloseq_16s_ar_path  <- file.path(run_16s_dir, "phyloseq_16s_ar.rds")

# 16S figures go to the supplement rather than the main figures/ folder.
supplementary_dir <- file.path(project_root, "supplementary")
dir.create(supplementary_dir, showWarnings = FALSE)
