# Exports sample name + measured O2 (mg/L) for the supplement.
# Sample names are sourced from the shared sample_id_map (this file's
# constants.R) rather than hardcoded again, so this always reflects the
# naming actually used throughout the project -- in particular the three
# E51/E52/E53 "B" samples use "_6_9cm_" (not "_B_") to match every other
# script and figure.
# Saves sample_name_map_o2.csv.
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
    "constants.R and the project's supplementary folder). ",
    "Run this script with Rscript, or open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
# -----------------------------------------------------------------------------

source(file.path(script_dir, "constants.R"))  # sample_id_map

supplementary_dir <- file.path(project_root, "supplementary")
dir.create(supplementary_dir, showWarnings = FALSE)

# Measured O2 (mg/L), keyed on the same raw sample IDs as sample_id_map.
# NA for T0 (time-zero) samples, which weren't measured.
o2_mgL <- c(
  "McM22_E41_T16_0-3cm_S1"              = 6.4,
  "McM22_E41_T16_3-6cm_S2"              = 6.8,
  "McM22_E41_T16_6-9cm_S3"              = 8.2,
  "McM22_E42_T16_0-3cm_S4"              = 7.3,
  "McM22_E42_T16_3-6cm_S5"              = 6.6,
  "McM22_E42_T16_6-9cm_S6"              = 7.0,
  "McM22_E43_T16_0-3cm_S7"              = 6.7,
  "McM22_E43_T16_3-6cm_S8"              = 7.8,
  "McM22_E43_T16_6-9cm_S9"              = 8.9,
  "McM22_E44_T16_0-3cm_S10"             = 6.6,
  "McM22_E44_T16_3-6cm_S11"             = 8.3,
  "McM22_E44_T16_6-9cm_S12"             = 7.3,
  "McM22_E4_T0_0-3cm_S13"               = NA,
  "McM22_E4_T0_3-6cm_S14"               = NA,
  "McM22_E4_T0_6-9cm_S15"               = NA,
  "McM22_E51_B_S16"                     = 0.07,
  "McM22_E52_B_S17"                     = 0.08,
  "McM22_E53_B_S18"                     = 8.02,
  "McM22_E1_A_S19"                      = 8.3,
  "McM22_E1_B_S20"                      = 9.2,
  "McM22_E31_T12_0-3_cm_S21"            = 7.4,
  "McM22_E32_T12_0-3_cm_S22"            = 7.32,
  "McM23_E12_T0_13CH4_S23"              = NA,
  "McM23_E12_T0_13CH4_O2_S24"           = NA,
  "McM23_E12_TF_CH4_S25"                = 3.69,
  "McM23_E12_TF_-CH4_S26"               = 1.07,
  "McM23_E12_TF_CH4_O2_S27"             = 5.12,
  "McM23_E12_TF_-CH4_O2_S28"            = 5.69,
  "McM23_E14_T0_13CH4_S29"              = NA,
  "McM23_E14_T0_13CH4_O2_S30"           = NA,
  "McM23_E14_TF_CH4_S31"                = 3.46,
  "McM23_E14_TF_-CH4_S32"               = 1.33,
  "McM23_E14_TF_CH4_O2_S33"             = 4.52,
  "McM23_E14_TF_-CH4_C6301_C6328O2_S34" = 6.56
)

sample_o2 <- data.frame(
  sample      = unname(sample_id_map),
  `O2 (mg/L)` = unname(o2_mgL[names(sample_id_map)]),
  check.names = FALSE
)

out_path <- file.path(supplementary_dir, "sample_name_map_o2.csv")
write.csv(sample_o2, out_path, row.names = FALSE, na = "")

message("Saved: ", out_path)
