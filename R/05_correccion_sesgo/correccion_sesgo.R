# correccion_sesgo_muestreo (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Panel A: fraccion de celulas recolectadas en cada gate (esfuerzo de muestreo)
# Panel B: scatter media-varianza y exclusion de promotores bimodales
#
# Requiere: data/processed/activity_stats_consistent.tsv (generado por
# R/01_activity_stats/build_activity_stats.R). Run from the TesisDoc repo root.

library(tidyverse)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/functions/perc_cells.R")

slug <- "correccion_sesgo_muestreo"
out_dir <- fig_dir(slug)
base_size <- 20

save_panel <- function(p, name) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = 9, height = 6.75, units = "in")
  }
}

# --- Panel A: esfuerzo de muestreo -----------------------------------------

panel_a <- perc_cells %>%
  ggplot(aes(sample, perc)) +
  geom_col(fill = thesis_clr) +
  xlab("Fracción de células") +
  ylab("Proporción de las células\ntotales en la fracción") +
  theme_pubclean(base_size = base_size) +
  theme(axis.text = element_text(colour = "black"))

save_panel(panel_a, "panel_a_sample_effort")

# --- Panel B: media vs varianza, exclusion de bimodales -------------------

stats_consistent <- read_tsv("data/processed/activity_stats_consistent.tsv", show_col_types = FALSE)

panel_b <- stats_consistent %>%
  arrange(bimodal) %>%
  mutate(bimodal = factor(bimodal, levels = c(FALSE, TRUE), labels = c("No bimodal", "Bimodal"))) %>%
  ggplot(aes(mean, var, col = bimodal)) +
  geom_point(size = 0.5) +
  facet_wrap(~rep) +
  scale_x_continuous(breaks = 1:6, limits = c(1, 6)) +
  ylab("Varianza") +
  xlab("Media") +
  scale_color_manual(values = c("No bimodal" = "#0D2C54", "Bimodal" = "#1B8C8E")) +
  labs(col = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  theme_pubr(base_size = base_size)

save_panel(panel_b, "panel_b_mean_var_bimodal")

# --- Figura compuesta ---------------------------------------------------------

compuesta <- (panel_a | panel_b) +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("correccion_sesgo_muestreo.", ext)), compuesta, width = 16, height = 6.75, units = "in")
}

message("Figura y paneles guardados en ", out_dir)
