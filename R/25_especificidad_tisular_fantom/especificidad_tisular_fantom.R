# especificidad_tisular_fantom (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Actividad endogena maxima (FANTOM5, todas las muestras) vs. actividad
# del reportero, separando promotores housekeeping (baja especificidad
# tisular) vs. tissue-specific (alta), segun tercil de sample_
# specificity_gini: tissue_specificty_boxplot.jpg (boxplots por 15 bins
# de actividad DENTRO de cada tercil) y tissue_specificity_residuals.jpg
# (residuos de un modelo lineal log10(tpm)~rank(actividad) por tercil,
# como medida de que tan bien la actividad del reportero predice la
# actividad endogena en cada grupo).
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv y data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("Tissue specificity" section).

library(tidyverse)
library(ggpubr)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "especificidad_tisular_fantom"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  mutate(max_tpm = replace_na(max_tpm, 0)) %>%
  filter(!is.na(sample_specificity_gini)) %>%
  group_by(rep) %>%
  mutate(specificity_tercile = cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) %>%
  ungroup()

# --- Boxplot: actividad endogena maxima por bin de actividad, dentro de cada tercil ---

panel_box <- data %>%
  filter(specificity_tercile != 2) %>%
  group_by(rep, specificity_tercile) %>%
  mutate(mean_bin = row_number(mean) %>% cut_number(15) %>% as.numeric()) %>%
  ungroup() %>%
  mutate(patron = ifelse(specificity_tercile == 1, "Housekeeping", "Tissue-specific")) %>%
  ggplot(aes(factor(mean_bin), max_tpm, fill = patron)) +
  geom_boxplot(outliers = FALSE, alpha = 0.6) +
  facet_wrap(~rep, scales = "free_x") +
  scale_y_log10() +
  labs(
    x = "Grupos ordenados de actividad media del reportero (bins de ~200)",
    y = "TPM (valor máximo a través de FANTOM5)", fill = "Patrón de expresión"
  ) +
  scale_fill_manual(values = c("Housekeeping" = "#1B8C8E", "Tissue-specific" = "#0D2C54")) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "tissue_specificty_boxplot.jpg"), panel_box, width = 9, height = 6.75, units = "in")

# --- Residuos de actividad endogena ~ rango de actividad, por tercil -----

lm_data <- data %>%
  filter(specificity_tercile != 2) %>%
  group_by(rep, specificity_tercile) %>%
  mutate(mean_rank = dense_rank(mean), log_tpm = log10(max_tpm + 0.1), grp = paste(rep, specificity_tercile, sep = "_")) %>%
  ungroup()

lm_l <- lm_data %>% split(.$grp) %>% map(~ lm(log_tpm ~ mean_rank, data = .x))

residuals_df <- lm_l %>%
  imap(~ tibble(residuals = residuals(.x), grp = .y) %>% separate(grp, into = c("rep", "specificity_tercile"), sep = "_")) %>%
  list_rbind() %>%
  mutate(patron = ifelse(specificity_tercile == "1", "Housekeeping", "Tissue-specific"))

panel_residuals <- residuals_df %>%
  ggplot(aes(residuals, fill = patron)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~rep) +
  scale_fill_manual(values = c("Housekeeping" = "#1B8C8E", "Tissue-specific" = "#0D2C54")) +
  labs(x = "Residuos (log10 TPM ~ rango de actividad)", y = "Densidad", fill = "Patrón de expresión") +
  ggpubr::theme_pubclean()
ggsave(file.path(out_dir, "tissue_specificity_residuals.jpg"), panel_residuals, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
