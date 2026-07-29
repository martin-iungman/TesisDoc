# reporter_endo_rank_dif (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Diferencia absoluta entre el rango de actividad del reportero y el
# rango de actividad endogena maxima (FANTOM5, todas las muestras),
# separando housekeeping vs. tissue-specific (tercil de sample_
# specificity_gini), contra un control de diferencias de rango
# permutadas al azar.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv y data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("Tissue specificity" section, fantom_max_expresion / rank_dif, lineas
# 799-848).

library(tidyverse)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

set.seed(1)

slug <- "reporter_endo_rank_dif"
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
  filter(specificity_tercile != 2) %>%
  mutate(
    rank_reporter = row_number(mean),
    rank_endo = row_number(max_tpm),
    patron = ifelse(specificity_tercile == 1, "Housekeeping", "Tissue-specific")
  ) %>%
  ungroup()

observed <- data %>%
  transmute(rep, patron, value = abs(rank_endo - rank_reporter))

control <- data %>%
  group_by(rep) %>%
  group_modify(~ tibble(value = abs(
    sample(seq_len(nrow(.x)), size = nrow(.x) * 1000, replace = TRUE) -
      sample(seq_len(nrow(.x)), size = nrow(.x) * 1000, replace = TRUE)
  ))) %>%
  ungroup() %>%
  mutate(patron = "Permutado")

plot_data <- bind_rows(observed, control) %>%
  mutate(patron = factor(patron, levels = c("Housekeeping", "Tissue-specific", "Permutado")))

p <- ggplot(plot_data, aes(value, col = patron, fill = patron, linetype = patron)) +
  geom_density(alpha = 0.5, linewidth = 0.7) +
  facet_wrap(~rep) +
  scale_color_manual(values = c("Housekeeping" = "#1B8C8E", "Tissue-specific" = "#0D2C54", "Permutado" = "darkred")) +
  scale_fill_manual(values = c("Housekeeping" = "#1B8C8E", "Tissue-specific" = "#0D2C54", "Permutado" = "transparent")) +
  scale_linetype_manual(values = c("Housekeeping" = "solid", "Tissue-specific" = "solid", "Permutado" = "dashed")) +
  labs(
    x = "Diferencia absoluta de rango (reportero vs. actividad endógena)",
    y = "Densidad", col = NULL, fill = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "reporter_endo_rank_dif.jpg"), p, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
