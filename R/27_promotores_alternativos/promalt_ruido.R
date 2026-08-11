# promalt_ruido (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Fig. 4B del paper. Ruido transcripcional (rango de varianza dentro de
# cada bin de actividad media) por categoria de promotor alternativo:
# "Sin promotores alternativos", "Promotor principal" o "Promotor
# secundario" - el secundario muestra mas ruido, particularmente en el
# rango medio de actividad.
#
# A diferencia del original, los bins de actividad (mean_sw) y el rango
# de varianza (var_rank_sw) se recalculan con los helpers compartidos
# add_mean_sw_bins()/add_var_rank_sw() (R/functions/plot_helpers.R, ya
# usados en toda la seccion R3) en vez de cut_width(mean_rank,100) -
# mismo criterio cualitativo (bins de ~100 promotores por rango de
# actividad, rango de varianza dentro de cada bin), simplificacion
# menor de como se define el borde del ultimo bin.
#
# Requiere: data/processed/activity_stats_highconf.tsv y data/processed/
# prom_alt_classification.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "Alternative promoters", lineas 1495-1513).

library(tidyverse)
library(ggpubr)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "promalt_ruido"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
promalt <- read_tsv("data/processed/prom_alt_classification.tsv", show_col_types = FALSE)

data_plot <- stats_highconf %>%
  left_join(promalt, by = "name", relationship = "many-to-many") %>%
  mutate(
    prom_alt2 = ifelse(prom_alt %in% c("independent", "correlated", "switch"), "Promotor secundario", prom_alt),
    prom_alt2 = ifelse(prom_alt2 == "unique", "Sin promotores alternativos", prom_alt2),
    prom_alt2 = ifelse(prom_alt2 == "Main promoter", "Promotor principal", prom_alt2)
  ) %>%
  filter(!is.na(prom_alt2), !prom_alt2 %in% c("non_detected", "unclassified")) %>%
  add_noise_rank()

data_plot %>%
  ggplot(aes(mean_sw, var_rank_sw, col = prom_alt2)) +
  geom_point(size = 0.2, alpha = 0.8) +
  facet_wrap(~rep) +
  labs(x = "Bin de actividad media (rango)", y = "Rango de varianza (dentro del bin)", col = "") +
  scale_color_manual(values = c(
    "Promotor principal" = "#16703A", "Promotor secundario" = "#5D9D1B",
    "Sin promotores alternativos" = "#B2D4FB"
  ), guide = guide_legend(override.aes = list(size = 3))) +
  geom_smooth(se = FALSE, aes(color = prom_alt2), show.legend = FALSE) +
  theme_pubr(base_size = 20)
ggsave(file.path(out_dir, "scatter_ruido_promalt.jpg"), width = 9, height = 6.75, units = "in")

message("Figura guardada en ", out_dir)
