# Scripts

> `.0`-numbered items (`2.0_processing_commands/`, `3.0_metagenome_processing_commands/`, `6.0_kaiju_commands/`, step_9's `diamond_cycdb_shortreads.sh`) are HPC commands, kept for provenance. Not runnable out of the box here.

## step_1_dic_script.R

Standalone. DIC (d13C) change per +CH4/-CH4 pairing, plotted by year with a site color strip. Also writes a supplementary table.
Run in RStudio (uses `rstudioapi` for its working directory), not `Rscript`.
Output: `figures/figure2_dic_rates.png/.pdf`, `supplementary/tableS5_dic_rates.csv`.

## step_2_16S_processing

* `fastp_run.sh`: adapter/quality trimming of raw 16S reads.
* `qiime_commands.txt`: QIIME2 commands for denoising, SILVA taxonomy, feature table.

## step_3_metagenomic_processing

* `pear_run.sh`: merges paired-end metagenomic reads.
* `prodigal_pear.sh`: ORF calling ahead of Kaiju (step_6) and DIAMOND (step_9).

## step_4_16S_visualization

* 4.1: builds the 16S phyloseq, applies sample naming, caches to `updated_16s_run/phyloseq_16s*.rds`.
* 4.2: top-15 genus bar plots (Bacteria/Archaea), same format as step_6.
* 4.3: genus PCoA (Bray-Curtis), 5-panel (Methane/Oxygen/Year/depth/Site) + PERMANOVA.
* 4.4: methanotroph bubble plot; genus names reconciled via step_8.
* `constants.R`/`helpers.R`: shared paths, taxa list, relative-abundance helpers.

## step_5_pcoa

* 5.1: metagenomic PCoA, faceted Methane/Oxygen/Year + depth/Site (two figures) + PERMANOVA. Reuses step_6's cached phyloseq.

## step_6_metagenomic_taxonomy

* `kaiju_short_reads.sh` etc.: Kaiju classification + summary tables.
* 6.1: builds the kaiju phyloseq, caches to `kaiju_output/kaiju_phyloseq.rds`. Run first -- used by 5.1, 6.2, 7.1, 7.2, 8.1, 9.1-9.3.
* 6.2: top-15 genus bar plots.
* 6.3: sample name + O2 crosswalk to `supplementary/`.
* `constants.R`/`helpers.R`: shared paths, sample name map, taxa list, plotting helpers (used across steps 4-7).

## step_7_methanotroph_visualization

* 7.1: methanotroph bubble plot, faceted by year, with location + metabolism bars.
* 7.2: Methanoperedens vs. Methylomirabilis -- trend/correlation plots, fold-change, residuals.

## step_8_comparing_16s_metagenomics

* `reconcile_phyloseq.R`: taxonomy synonym map + lineage corrections across SILVA/NCBI/GTDB.
* 8.1: reconciles 16S and metagenomic taxonomy, compares methanotroph detection and taxa counts. Bar plot, detection heatmap, supplementary tables.
* 8.2: exports the reconciliation tables to `supplementary/`.

## step_9_metagenomic_function

* `diamond_cycdb_shortreads.sh`: DIAMOND vs. cycDB.
* 9.1: parses DIAMOND output, maps genes to pathways, caches pathway x sample counts. Run before 9.2/9.3.
* 9.2: DESeq2 pathway differential abundance (full model + within-year). Results + lollipop plots to `deseq_output/`.
* 9.3: VST heatmap of pathway abundance, annotated by Year/Site/sediment/depth/Methane/Oxygen.
* `constants.R`/`helpers.R`: shared paths, pathway gene lists, DIAMOND-parsing helpers.
