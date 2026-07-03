# Shared constants for the kaiju taxonomy/visualization scripts
# (step_6_metagenomic_taxonomy and step_7_methanotroph_visualization).
#
# Requires `project_root` to already be defined in the calling environment
# (each driver script sets it from its own location before sourcing this file).

kaiju_dir           <- file.path(project_root, "kaiju_output", "genus_full_taxonomy")
kaiju_phyloseq_path <- file.path(project_root, "kaiju_output", "kaiju_phyloseq.rds")
figures_dir          <- file.path(project_root, "figures")

taxa_of_interest <- c(
  "Methylococcus", "Methylomonas", "Methylobacter", "Methylomicrobium",
  "Methylosarcina", "Methylocaldum", "methylogea", "Methylosoma",
  "Methyloparacoccus", "Methyloglobulus", "Methyloprofundus", "Methylomarinum",
  "Methylovolum", "Methylomagnum", "Methylosphaera", "Methylothermus",
  "Methyloholbius", "Methylomarinovum", "Methylosinus", "Methylocystis",
  "Methylocella", "Methylocaspa", "Methyloferula",
  "Candidatus Methanoperedens", "Candidatus Methylomirabilis",
  "Methylacidiphilum", "Methylacidimicrobium"
)

# Bubble/bar plots order genera the same way they're subset
genus_order <- taxa_of_interest

sample_list <- c(
  "E3_0_3cm_+CH4_+O2", "E3_0_3cm_-CH4_+O2",
  "E1R1_0_3cm_-CH4_+O2", "E1R1_3_6cm_-CH4_+O2", "E1R1_6_9cm_-CH4_+O2",
  "E1R2_0_3cm_-CH4_+O2", "E1R2_3_6cm_-CH4_+O2", "E1R2_6_9cm_-CH4_+O2",
  "E1R1_0_3cm_+CH4_+O2", "E1R1_3_6cm_+CH4_+O2", "E1R1_6_9cm_+CH4_+O2",
  "E1R2_0_3cm_+CH4_+O2", "E1R2_3_6cm_+CH4_+O2", "E1R2_6_9cm_+CH4_+O2",
  "E1T0_0_3cm", "E1T0_3_6cm", "E1T0_6-9cm",
  "E2R1_6_9cm_+CH4_-O2", "E2R2_6_9cm_+CH4_-O2", "E2R3_6_9cm_+CH4_-O2",
  "E4_0_3cm_+CH4_+O2", "E4_0_3cm_-CH4_+O2",
  "E5T0_0_4cm", "E5T0_4_8cm", "E5_0_4cm_+CH4_+O2", "E5_0_4cm_-CH4_+O2",
  "E5_4_8cm_+CH4_-O2", "E5_4_8cm_-CH4_-O2",
  "E6T0_0_4cm", "E6T0_4_8cm", "E6_0_4cm_+CH4_+O2", "E6_0_4cm_-CH4_+O2",
  "E6_4_8cm_+CH4_-O2", "E6_4_8cm_-CH4_-O2"
)

# Raw McM2X_... sample ID (as extracted from the kaiju "file" column) -> updated name.
# E2 "B" samples use "_6_9cm_" rather than "_B_" to match prior figures/manuscript text.
sample_id_map <- c(
  "McM22_E41_T16_0-3cm_S1"              = "E1R1_0_3cm_-CH4_+O2",
  "McM22_E41_T16_3-6cm_S2"              = "E1R1_3_6cm_-CH4_+O2",
  "McM22_E41_T16_6-9cm_S3"              = "E1R1_6_9cm_-CH4_+O2",
  "McM22_E42_T16_0-3cm_S4"              = "E1R2_0_3cm_-CH4_+O2",
  "McM22_E42_T16_3-6cm_S5"              = "E1R2_3_6cm_-CH4_+O2",
  "McM22_E42_T16_6-9cm_S6"              = "E1R2_6_9cm_-CH4_+O2",
  "McM22_E43_T16_0-3cm_S7"              = "E1R1_0_3cm_+CH4_+O2",
  "McM22_E43_T16_3-6cm_S8"              = "E1R1_3_6cm_+CH4_+O2",
  "McM22_E43_T16_6-9cm_S9"              = "E1R1_6_9cm_+CH4_+O2",
  "McM22_E44_T16_0-3cm_S10"             = "E1R2_0_3cm_+CH4_+O2",
  "McM22_E44_T16_3-6cm_S11"             = "E1R2_3_6cm_+CH4_+O2",
  "McM22_E44_T16_6-9cm_S12"             = "E1R2_6_9cm_+CH4_+O2",
  "McM22_E4_T0_0-3cm_S13"               = "E1T0_0_3cm",
  "McM22_E4_T0_3-6cm_S14"               = "E1T0_3_6cm",
  "McM22_E4_T0_6-9cm_S15"               = "E1T0_6-9cm",
  "McM22_E51_B_S16"                     = "E2R1_6_9cm_+CH4_-O2",
  "McM22_E52_B_S17"                     = "E2R2_6_9cm_+CH4_-O2",
  "McM22_E53_B_S18"                     = "E2R3_6_9cm_+CH4_-O2",
  "McM22_E1_A_S19"                      = "E3_0_3cm_+CH4_+O2",
  "McM22_E1_B_S20"                      = "E3_0_3cm_-CH4_+O2",
  "McM22_E31_T12_0-3_cm_S21"            = "E4_0_3cm_+CH4_+O2",
  "McM22_E32_T12_0-3_cm_S22"            = "E4_0_3cm_-CH4_+O2",
  "McM23_E12_T0_13CH4_S23"              = "E5T0_4_8cm",
  "McM23_E12_T0_13CH4_O2_S24"           = "E5T0_0_4cm",
  "McM23_E12_TF_CH4_S25"                = "E5_4_8cm_+CH4_-O2",
  "McM23_E12_TF_-CH4_S26"               = "E5_4_8cm_-CH4_-O2",
  "McM23_E12_TF_CH4_O2_S27"             = "E5_0_4cm_+CH4_+O2",
  "McM23_E12_TF_-CH4_O2_S28"            = "E5_0_4cm_-CH4_+O2",
  "McM23_E14_T0_13CH4_S29"              = "E6T0_4_8cm",
  "McM23_E14_T0_13CH4_O2_S30"           = "E6T0_0_4cm",
  "McM23_E14_TF_CH4_S31"                = "E6_4_8cm_+CH4_-O2",
  "McM23_E14_TF_-CH4_S32"               = "E6_4_8cm_-CH4_-O2",
  "McM23_E14_TF_CH4_O2_S33"             = "E6_0_4cm_+CH4_+O2",
  "McM23_E14_TF_-CH4_C6301_C6328O2_S34" = "E6_0_4cm_-CH4_+O2"
)
