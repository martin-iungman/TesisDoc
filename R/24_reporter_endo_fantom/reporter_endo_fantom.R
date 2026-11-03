# Exploratorio (no va a la tesis): actividad del reportero vs. actividad
# endogena promedio a traves de TODAS las muestras de FANTOM5 (tejidos +
# celulas primarias): reporter_endo_mean_fantom.jpg (mediana de la media
# por bin de actividad) y reporter_endo_mean_fantom_boxplot.jpg (boxplot
# por 20 bins de rango de actividad). Salida en
# figures/exploratorio_reporter_endo_fantom. Era la figura
# reporter_endo_fantom; su boxplot de la actividad endogena maxima paso a
# actividad_endogena_multi_tejido (panel D, R/23_actividad_endogena).
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv y data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# lee Analysis/Tables, ver R/00_prom_features/analysis_tables_exceptions.R).
# Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("Tissue specificity" section, fantom_mean_expresion).

library(tidyverse)
source("R/functions/plot_helpers.R")

out_dir <- "figures/exploratorio_reporter_endo_fantom"
dir.create(out_dir, showWarnings = FALSE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  mutate(across(c(mean_tpm, max_tpm), ~ replace_na(.x, 0)))
data <- add_mean_sw_bins(data)

# --- Mediana de actividad endogena media por bin de actividad -------------

p1 <- data %>%
  group_by(rep, mean_sw) %>%
  summarise(median_endo = median(mean_tpm), .groups = "drop") %>%
  ggplot(aes(mean_sw, median_endo)) +
  geom_point(col = "#AD343E") +
  geom_smooth(col = "#216869", method = "lm") +
  facet_wrap(~rep) +
  scale_y_log10() +
  labs(x = "Actividad media (bins por rango)", y = "Mediana de actividad media endógena (TPM)") +
  theme_bw(base_size = 20)
ggsave(file.path(out_dir, "reporter_endo_mean_fantom.jpg"), p1, width = 9, height = 6.75, units = "in")

# --- Boxplots por 20 bins de rango de actividad ----------------------------

tpm_rank_boxplot <- function(data, tpm_col, ylab) {
  data %>%
    group_by(rep) %>%
    mutate(bin = row_number(mean) %>% cut_number(20) %>% as.numeric() %>% as_factor()) %>%
    ggplot(aes(bin, .data[[tpm_col]])) +
    geom_boxplot(fill = "#AD343E", alpha = 0.5, outliers = FALSE) +
    facet_wrap(~rep) +
    scale_y_log10() +
    labs(x = "Actividad del reportero (bins por rango)", y = ylab) +
    theme_bw(base_size = 20)
}

p2 <- tpm_rank_boxplot(data, "mean_tpm", "Actividad endógena media entre tejidos (TPM)")
ggsave(file.path(out_dir, "reporter_endo_mean_fantom_boxplot.jpg"), p2, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
