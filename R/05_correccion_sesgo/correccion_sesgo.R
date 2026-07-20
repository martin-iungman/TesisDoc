# correccion_sesgo_muestreo (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Panel A: factores de correccion por esfuerzo de muestreo diferencial
# Panel B: scatter media-varianza y exclusion de promotores bimodales
#
# Requiere: data/processed/activity_stats_consistent.tsv (generado por
# R/01_activity_stats/build_activity_stats.R). Run from the TesisDoc repo root.

library(tidyverse)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/functions/perc_cells.R")

slug <- "correccion_sesgo_muestreo"
out_dir <- fig_dir(slug)

# --- Panel A: esfuerzo de muestreo -----------------------------------------

panel_a <- perc_cells %>%
  ggplot(aes(sample, perc)) +
  geom_col(fill = "#FFB400") +
  xlab("Población de células") +
  ylab("Factor de relativización") +
  theme_pubclean()

ggsave(file.path(out_dir, "panel_a_sample_effort.jpg"), panel_a, width = 9, height = 6.75, units = "in")

# --- Panel B: media vs varianza, exclusion de bimodales -------------------

stats_consistent <- read_tsv("data/processed/activity_stats_consistent.tsv", show_col_types = FALSE)

panel_b <- stats_consistent %>%
  arrange(bimodal) %>%
  ggplot(aes(mean, var, col = bimodal)) +
  geom_point(size = 0.5) +
  facet_wrap(~rep) +
  ylab("Varianza") +
  xlab("Media") +
  theme_pubr() +
  scale_color_manual(values = c("FALSE" = "#0D2C54", "TRUE" = "#1B8C8E")) +
  labs(col = "Bimodal") +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme(text = element_text(size = 20))

ggsave(file.path(out_dir, "panel_b_mean_var_bimodal.jpg"), panel_b, width = 9, height = 6.75, units = "in")

message("Paneles A y B guardados en ", out_dir)
