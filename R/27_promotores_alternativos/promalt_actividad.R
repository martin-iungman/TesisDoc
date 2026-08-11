# promalt_actividad (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Fig. 4A del paper (mitad derecha - el esquema a mano de la izquierda,
# Only/Main/Secondary, queda pendiente). Actividad media del reportero
# por categoria de promotor alternativo: "Sin promotores alternativos"
# (unique - un solo promotor concentra >99.9% de la actividad del gen),
# "Promotor principal" (Main) o "Promotor secundario" (agrupa switch/
# correlated/independent), con comparaciones de Wilcoxon entre las tres
# categorias.
#
# Requiere: data/processed/activity_stats_highconf.tsv (ver
# R/01_activity_stats/build_activity_stats.R) y data/processed/
# prom_alt_classification.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "Alternative promoters", lineas 1410-1437).

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "promalt_actividad"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
promalt <- read_tsv("data/processed/prom_alt_classification.tsv", show_col_types = FALSE)

# 62 promotores aparecen 2 veces en prom_alt_classification.tsv: son a
# la vez "unclassified" (clasificacion primaria ambigua, cumplen solo
# uno de los dos criterios de "Main") y "independent"/"correlated"/
# "switch" (su clasificacion secundaria contra el Main real del gen) -
# no son categorias excluyentes en el diseno original. El filtro de
# abajo descarta la fila "unclassified" en todos los casos, asi que no
# hay duplicados en el resultado final.
data_plot <- stats_highconf %>%
  left_join(promalt, by = "name", relationship = "many-to-many") %>%
  filter(!is.na(prom_alt), !prom_alt %in% c("unclassified", "non_detected")) %>%
  mutate(
    prom_alt2 = ifelse(prom_alt %in% c("independent", "correlated", "switch"), "Promotor secundario", prom_alt),
    prom_alt2 = ifelse(prom_alt2 == "unique", "Sin promotores alternativos", prom_alt2),
    prom_alt2 = ifelse(prom_alt2 == "Main promoter", "Promotor principal", prom_alt2),
    prom_alt2 = factor(prom_alt2, levels = c("Sin promotores alternativos", "Promotor principal", "Promotor secundario"))
  )

p <- ggviolin(data_plot,
  x = "prom_alt2", y = "mean", fill = "prom_alt2", draw_quantiles = 0.5,
  add = "median_q1q3", alpha = 0.85,
  palette = c(
    "Promotor principal" = "#35B166", "Promotor secundario" = "#5D9D1B",
    "Sin promotores alternativos" = "#B2D4FB"
  )
)
p <- facet(p, facet.by = "rep", ncol = 2)
p +
  stat_compare_means(
    comparisons = list(
      c("Promotor principal", "Promotor secundario"),
      c("Promotor principal", "Sin promotores alternativos"),
      c("Sin promotores alternativos", "Promotor secundario")
    ),
    label = "p.signif", method = "wilcox.test", na.rm = TRUE
  ) +
  theme_pubclean() +
  labs(fill = "Promotores alternativos", x = NULL, y = "Actividad media transcripcional") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size = 15)) +
  ylim(0, 7.5)
ggsave(file.path(out_dir, "violin_actividad_promalt.jpg"), width = 12, height = 6.75, units = "in")

message("Figura guardada en ", out_dir)
