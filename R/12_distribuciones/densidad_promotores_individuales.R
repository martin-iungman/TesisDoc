# densidad_promotores_individuales (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Densidad de fluorescencia EGFP (citometria) para cada uno de los 8
# promotores validados individualmente, contra el control ("US"/Control)
# de fondo.
#
# Resync 2026-07-29: separado de distribuciones_expresion (que paso a Fig.
# R1.6) en su propia figura - antes era un panel "extra" sin numero propio.
#
# Requiere: ~519MB de datos de citometria leidos directo de
# transcriptional_library (ver R/00_prom_features/heavy_data_paths.R).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "densidad_promotores_individuales"
out_dir <- fig_dir(slug)

df <- load_citometry_stable_validation(path_citometry_stable_validation)

p <- df %>%
  filter(name != "Control") %>%
  ggplot() +
  geom_density(mapping = aes(Comp_FL2_A, col = name)) +
  scale_x_log10() +
  facet_wrap(~name) +
  geom_density(data = df %>% filter(name == "Control") %>% select(-name), mapping = aes(Comp_FL2_A), col = "grey") +
  xlab("Señal de EGFP\n(unidades de fluorescencia relativa)") +
  ylab("Densidad") +
  theme_pubclean() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "densidad_individual_promotores.jpg"), p, width = 9, height = 6.75, units = "in")

message("Figura guardada en ", out_dir)
