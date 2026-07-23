# shape_especificidad_tisular (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Panel A: densidad del shape del promotor (ancho interquantile, CAGEr).
# Panel B: densidad de la especificidad tisular (indice de Gini, FANTOM5).
# En ambos, lineas punteadas marcando los cortes de terciles (mismo
# criterio usado para clasificar "angosto/intermedio/ancho" y "baja/
# media/alta especificidad" en R/09_coocurrencia y R/15_summary_features:
# cut_number(x, 3)).
#
# Requiere data/processed/prom_df.tsv (ver R/00_prom_features/
# build_prom_features.R - interquantile_width y sample_specificity_gini
# vienen de las excepciones documentadas en
# R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.

library(tidyverse)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "shape_especificidad_tisular"
out_dir <- fig_dir(slug)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

density_terciles <- function(data, var, xlab, fill) {
  x <- data[[var]]
  cortes <- quantile(x, c(1 / 3, 2 / 3), na.rm = TRUE)
  ggplot(data, aes(.data[[var]])) +
    geom_density(fill = fill, alpha = 0.6, col = fill, linewidth = 0.8) +
    geom_vline(xintercept = cortes, linetype = "dashed", linewidth = 0.5, col = "#333333") +
    labs(x = xlab, y = "Densidad") +
    theme_pubclean(base_size = 16)
}

# --- Panel A: shape del promotor (ancho interquantile) ---------------------

panel_a <- density_terciles(prom_df, "interquantile_width", "Ancho interquantile (pb)", thesis_clr)

# --- Panel B: especificidad tisular (indice de Gini) -----------------------

panel_b <- density_terciles(prom_df, "sample_specificity_gini", "Especificidad tisular (índice de Gini)", "#AD343E")

# --- Save individual panels + combined figure -----------------------------

ggsave(file.path(out_dir, "panel_a_shape.jpg"), panel_a, width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "panel_b_especificidad_tisular.jpg"), panel_b, width = 9, height = 6.75, units = "in")

combined <- panel_a + panel_b + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, paste0(slug, ".jpg")), combined, width = 14, height = 6, units = "in")

message("Figuras guardadas en ", out_dir)
