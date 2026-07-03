# Methanoperedens vs. Methylomirabilis: trend line plot, correlation plot,
# and summary stats (fold-change, depth means, residuals).
# Saves methanoperedens_methylomirabilis.png.
# Self-locating: run with Rscript, or source/run in RStudio -- no manual
# setwd() needed. Requires 6.1_load_kaiju_data.R to have been run at least once.

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
    "step_6_metagenomic_taxonomy/constants.R+helpers.R and the project's ",
    "kaiju_output/figures folders). Run this script with Rscript, or ",
    "open/source it in RStudio."
  )
}

script_dir   <- .get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
taxonomy_dir <- file.path(project_root, "scripts", "step_6_metagenomic_taxonomy")
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)

source(file.path(taxonomy_dir, "constants.R"))
source(file.path(taxonomy_dir, "helpers.R"))

loaded <- load_relative_abundance_long(kaiju_phyloseq_path, taxa_of_interest, sample_list, genus_order)
relative_abundance_long <- loaded$data
sample_order <- loaded$sample_order

# Matched by substring so it's robust to the "Candidatus " prefix / minor naming variants.
mp_name <- grep("Methanoperedens",  levels(relative_abundance_long$genus), value = TRUE)[1]
mm_name <- grep("Methylomirabilis", levels(relative_abundance_long$genus), value = TRUE)[1]

two_genus <- relative_abundance_long %>%
  filter(genus %in% c(mp_name, mm_name)) %>%
  droplevels() %>%
  # Fill samples where a genus is absent with 0 so the lines stay continuous.
  complete(
    Sample = factor(sample_order, levels = sample_order),
    genus,
    fill = list(RelativeAbundance = 0)
  ) %>%
  group_by(Sample) %>%
  tidyr::fill(Year, Location, .direction = "downup") %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = sample_order))

p_lines <- ggplot(two_genus, aes(x = Sample, y = RelativeAbundance)) +
  geom_tile(aes(fill = Location), alpha = 0, show.legend = TRUE) +
  geom_line(aes(color = genus, group = genus), linewidth = 0.6) +
  geom_point(aes(color = genus, group = genus), size = 2.5, shape = 16) +
  scale_color_manual(values = setNames(c("#7570B3", "#D95F02"), c(mp_name, mm_name)), name = "Genus") +
  scale_fill_manual(values = c("CCS" = "#1B9E77", "Jetty" = "#F0C808"), name = "Location") +
  scale_y_continuous(labels = scales::label_number()) +
  guides(color = guide_legend(order = 1), fill = guide_legend(order = 2, override.aes = list(alpha = 1))) +
  facet_grid(~ Year, scales = "free_x", space = "free_x") +
  labs(x = "", y = "Relative abundance") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    axis.text.y      = element_text(size = 12),
    axis.title.y     = element_text(size = 14),
    legend.title     = element_text(size = 13, face = "bold"),
    legend.text      = element_text(face = "italic", size = 12),
    strip.background = element_rect(fill = "grey80", color = NA),
    strip.text       = element_text(face = "bold", size = 13)
  )

wide <- two_genus %>%
  group_by(Sample, Year, Location, genus) %>%
  summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop") %>%
  pivot_wider(names_from = genus, values_from = RelativeAbundance, values_fill = 0)

# Grab the two genus columns by name regardless of exact spelling.
mp_col <- grep("Methanoperedens",  names(wide), value = TRUE)[1]
mm_col <- grep("Methylomirabilis", names(wide), value = TRUE)[1]
wide$mp <- wide[[mp_col]]
wide$mm <- wide[[mm_col]]

p_corr <- ggplot(wide, aes(x = mp, y = mm)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", fill = "grey85") +
  geom_point(aes(color = Location, shape = Year), size = 3) +
  scale_color_manual(values = c("CCS" = "#1B9E77", "Jetty" = "#F0C808"), name = "Location") +
  scale_x_continuous(labels = scales::label_number()) +
  scale_y_continuous(labels = scales::label_number()) +
  labs(x = paste0(mp_name, " RA"), y = paste0(mm_name, " RA")) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(face = "italic", size = 14),
    axis.title.y = element_text(face = "italic", size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 12)
  )

loc_bar_two <- two_genus %>%
  distinct(Sample, Location, Year) %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%
  ggplot(aes(x = Sample, y = 1, fill = Location)) +
  geom_tile(width = 1) +
  scale_fill_manual(values = c("CCS" = "#1B9E77", "Jetty" = "#F0C808"), name = "Location") +
  scale_x_discrete(expand = c(0, 0)) +
  facet_grid(~ Year, scales = "free_x", space = "free_x") +
  theme_void() +
  theme(legend.position = "none", strip.text = element_blank(), plot.margin = margin(0, 5, 0, 5))

final_plot <- loc_bar_two / p_lines / p_corr +
  plot_layout(heights = c(0.06, 1, 1.2)) +
  plot_annotation(
    tag_levels = list(c("", "A", "B")),
    theme = theme(plot.tag = element_text(size = 16, face = "bold"))
  )

final_plot

ggsave(
  file.path(figures_dir, "Figure6_methanoperedens_methylomirabilis.png"),
  final_plot, width = 12, height = 10, dpi = 300
)
ggsave(
  file.path(figures_dir, "Figure6_methanoperedens_methylomirabilis.pdf"),
  final_plot, width = 12, height = 10, dpi = 300
)

ct_s <- suppressWarnings(cor.test(wide$mp, wide$mm, method = "spearman"))
cat(sprintf("Spearman rho = %.3f, p = %.3g", unname(ct_s$estimate), ct_s$p.value))

fold_change_2023 <- two_genus %>%
  filter(Year == "2023") %>%
  mutate(Depth = case_when(
    str_detect(Sample, "0_4cm") ~ "shallow",  # 0-4 cm
    str_detect(Sample, "4_8cm") ~ "deep",     # 4-8 cm
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Depth)) %>%
  group_by(genus, Depth) %>%
  summarise(mean_RA = mean(RelativeAbundance), .groups = "drop") %>%
  pivot_wider(names_from = Depth, values_from = mean_RA) %>%
  mutate(fold_change = deep / shallow)

print(fold_change_2023)

depth_means_2022 <- two_genus %>%
  filter(Year == "2022") %>%
  mutate(Depth = case_when(
    str_detect(Sample, "0_3cm") ~ "0-3 cm",
    str_detect(Sample, "3_6cm") ~ "3-6 cm",
    str_detect(Sample, "6_9cm|6-9cm") ~ "6-9 cm",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Depth)) %>%
  group_by(genus, Depth) %>%
  summarise(mean_RA = mean(RelativeAbundance), .groups = "drop") %>%
  pivot_wider(names_from = Depth, values_from = mean_RA)

print(depth_means_2022)

# Fit the same model the correlation plot draws.
fit <- lm(mm ~ mp, data = wide)

resid_tbl <- wide %>%
  mutate(
    predicted = predict(fit),
    residual  = mm - predicted,  # + = above line, - = below
    abs_resid = abs(residual)
  ) %>%
  select(Sample, Location, Year, mp, mm, predicted, residual, abs_resid) %>%
  arrange(desc(abs_resid))

print(resid_tbl, n = Inf)

resid_tbl %>%
  group_by(Location) %>%
  summarise(mean_abs_resid = mean(abs_resid), max_abs_resid = max(abs_resid), n = n(), .groups = "drop")

resid_tbl %>%
  summarise(mean_abs_resid = mean(abs_resid), max_abs_resid = max(abs_resid), n = n(), .groups = "drop")
