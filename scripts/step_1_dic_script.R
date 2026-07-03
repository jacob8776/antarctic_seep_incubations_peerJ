setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


data <- read.csv("../dic_rates.csv")


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(stringr)     # str_detect()
library(patchwork)   # plot_layout() / combining plots


data <- data %>%
  mutate(norm_d13C = d13C / (methane.1 * Sediment_g_mL))

# Difference in normalized d13C (Yes - No) per pairing.
# reframe() returns one row per pairing, so `pairing` is carried through.
diff_data <- data %>%
  group_by(pairing) %>%
  filter(any(methane == "Yes") & any(methane == "No" | methane == "no")) %>%
  reframe(diff_d13C = (norm_d13C[methane == "Yes"] - norm_d13C[methane == "No" | methane == "no"])/Days[methane == "Yes"])

# Figure labels, keyed on the pairing id (1-14). Nothing here depends on the
# #NAME or Name columns in the file, so renaming those can't break the plot.
pairing_labels <- c(
  `1`  = "E1R1_0_3cm_+O2",
  `2`  = "E1R1_3_6cm_+O2",
  `3`  = "E1R1_6_9cm_+O2",
  `4`  = "E1R2_0_3cm_+O2",
  `5`  = "E1R2_3_6cm_+O2",
  `6`  = "E1R2_6_9cm_+O2",
  `7`  = "E3_0_3cm_+O2",
  `8`  = "E4_0_3cm_+O2",
  `9`  = "E5_4_8cm_-O2*",
  `10` = "E5_0_4cm_+O2*",
  `11` = "E6_4_8cm_-O2*",
  `12` = "E6_0_4cm_+O2*",
  `13` = "E5E6_FSWC_-O2",
  `14` = "E5E6_FSWC_+O2"
)

diff_data$Name <- unname(pairing_labels[as.character(diff_data$pairing)])
if (any(is.na(diff_data$Name)))
  warning("No label for pairing(s): ",
          paste(diff_data$pairing[is.na(diff_data$Name)], collapse = ", "))

# Keep the full set (incl. FSWC controls) for the supplementary table at the end
all_diffs <- diff_data

# Quick look
ggplot(diff_data, aes(x = factor(Name), y = diff_d13C)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Experiment", y = expression(Delta*delta^13*C~DIC))

# Drop seawater controls
diff_data <- filter(diff_data, Name != "E5E6_FSWC_-O2", Name != "E5E6_FSWC_+O2")

# Bar order
diff_data$Name <- factor(diff_data$Name, levels = c("E3_0_3cm_+O2",
                                                    "E1R1_0_3cm_+O2",
                                                    "E1R1_3_6cm_+O2",
                                                    "E1R1_6_9cm_+O2",
                                                    "E1R2_0_3cm_+O2",
                                                    "E1R2_3_6cm_+O2",
                                                    "E1R2_6_9cm_+O2",
                                                    "E4_0_3cm_+O2",
                                                    "E5_0_4cm_+O2*",
                                                    "E5_4_8cm_-O2*",
                                                    "E6_0_4cm_+O2*",
                                                    "E6_4_8cm_-O2*"))


p <- ggplot(diff_data, aes(x = Name, y = diff_d13C)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(x = "", y = "sediment") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(p)

mean(diff_data$diff_d13C)
median(diff_data$diff_d13C)


plot_data <- diff_data %>%
  mutate(
    Year = case_when(str_detect(Name, "E5|E6") ~ "2023", TRUE ~ "2022"),
    Location = ifelse(str_detect(Name, "E3"), "Jetty", "CCS")
  )

y_label_complex <- expression(Delta * delta^13 * C ~ DIC ~x~ (g~or~mL ~ sediment~x~days)^-1)

p_main <- ggplot(plot_data, aes(y = Name, x = diff_d13C)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.7) +
  facet_grid(Year ~ ., scales = "free", space = "free_y") +
  scale_y_discrete(limits = rev) +
  labs(x = y_label_complex, y = "", title = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 12)
  )

p_strip <- ggplot(plot_data, aes(y = Name, x = 1, fill = Location)) +
  geom_tile() +
  facet_grid(Year ~ ., scales = "free", space = "free_y") +
  scale_fill_manual(values = c("Jetty" = "#FFD700", "CCS" = "#2E8B57")) +
  scale_y_discrete(limits = rev) +
  theme_void() +
  theme(strip.text = element_blank(), legend.position = "right")

final_plot <- p_main + p_strip + plot_layout(widths = c(10, 0.5))
print(final_plot)

ggsave("../figures/figure2_dic_rates.png", plot = final_plot, width = 5, height = 5, dpi = 300)
ggsave("../figures/figure2_dic_rates.pdf", plot = final_plot, width = 5, height = 5, dpi = 300, device = cairo_pdf)


# --- Supplementary table: DIC change per pairing (INCLUDING FSWC controls) -----
# Built from all_diffs (the full 14-pairing set, before controls were dropped).
# `d13C_DIC_change` is exactly the value plotted (Yes-No normalized d13C per day).
supp_order <- c("E3_0_3cm_+O2",
                "E1R1_0_3cm_+O2", "E1R1_3_6cm_+O2", "E1R1_6_9cm_+O2",
                "E1R2_0_3cm_+O2", "E1R2_3_6cm_+O2", "E1R2_6_9cm_+O2",
                "E4_0_3cm_+O2",
                "E5_0_4cm_+O2*", "E5_4_8cm_-O2*", "E6_0_4cm_+O2*", "E6_4_8cm_-O2*",
                "E5E6_FSWC_-O2", "E5E6_FSWC_+O2")   # controls last

supp_table <- all_diffs %>%
  transmute(
    Experiment      = as.character(Name),
    Year            = ifelse(str_detect(Name, "E5|E6"), "2023", "2022"),
    d13C_DIC_change = round(diff_d13C, 3)
  )
supp_table <- supp_table[match(supp_order, supp_table$Experiment), ]
rownames(supp_table) <- NULL

print(supp_table)

# --- Supplementary CSV ---------------------------------------------------------
write.csv(supp_table, "../supplementary/tableS5_dic_rates.csv", row.names = FALSE)
