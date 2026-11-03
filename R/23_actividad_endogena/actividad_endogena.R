# actividad_endogena_multi_tejido (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Actividad del promotor endogeno (CAGE TPM, FANTOM5) vs. actividad del
# promotor reportero, por bin de actividad (mean_sw) y replica: HEK293 (A),
# HeLa (B) y musculo esqueletico (C); y actividad endogena maxima entre
# todas las muestras de FANTOM5 (tejidos + celulas primarias) por 20 bins
# de rango de actividad del reportero (D, boxplot). Figura completa
# actividad_endogena_multi_tejido.jpg/.pdf (2x2).
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv (hek_tpm, ver R/00_prom_features/build_prom_features.R),
# data/processed/endo_cage_activity.tsv (hela_tpm, muscle_tpm, ver
# R/00_prom_features/build_endo_cage_activity.R) y data/processed/
# fantom_endo_activity_summary.tsv (max_tpm, ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION, lee
# Analysis/Tables, ver R/00_prom_features/analysis_tables_exceptions.R).
# Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/cell_line_cage.qmd
# (HeLa y musculo esqueletico), hek_cage.qmd (HEK293) y final_github.R
# ("Tissue specificity", fantom_mean_expresion; el boxplot venia de
# R/24_reporter_endo_fantom, que ahora solo deja exploratorios).

library(tidyverse)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "actividad_endogena_multi_tejido"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
endo_cage <- read_tsv("data/processed/endo_cage_activity.tsv", show_col_types = FALSE)
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(endo_cage, by = "seq_id") %>%
  left_join(fantom_summary, by = "seq_id") %>%
  mutate(across(c(hek_tpm, hela_tpm, muscle_tpm, max_tpm), ~ replace_na(.x, 0)))
data <- add_mean_sw_bins(data)

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}

p_hek <- endo_activity_scatter(data, "hek_tpm", "Línea celular HEK293")
p_hela <- endo_activity_scatter(data, "hela_tpm", "Línea celular HeLa")
p_muscle <- endo_activity_scatter(data, "muscle_tpm", "Músculo esquelético")

save_jpg_pdf(p_hek, "HEK293_mean", 9, 6.75)
save_jpg_pdf(p_hela, "Hela_mean", 9, 6.75)
save_jpg_pdf(p_muscle, "skeletal_muscle_mean", 9, 6.75)

# Boxplot de la actividad endogena maxima por 20 bins de rango de actividad
# del reportero (sin outliers); en el eje x solo algunos bins para que las
# etiquetas no se pisen.
p_max <- data %>%
  group_by(rep) %>%
  mutate(bin = row_number(mean) %>% cut_number(20) %>% as.numeric() %>% as_factor()) %>%
  ggplot(aes(bin, max_tpm)) +
  geom_boxplot(fill = "#AD343E", alpha = 0.5, outliers = FALSE) +
  facet_wrap(~rep) +
  scale_y_log10() +
  scale_x_discrete(breaks = c("1", "5", "10", "15", "20")) +
  labs(
    x = "Bins de actividad del promotor reportero (por rango)",
    y = "TPM (valor máximo entre\nmuestras CAGE FANTOM5)",
    title = "Máximo entre muestras FANTOM5"
  ) +
  theme_bw(base_size = 16) +
  theme(panel.spacing = unit(1.5, "lines"))
save_jpg_pdf(p_max, "reporter_endo_max_fantom_boxplot", 9, 6.75)

figura <- p_hek + p_hela + p_muscle + p_max +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
save_jpg_pdf(figura, slug, 18, 13.5)

message("Figuras guardadas en ", out_dir)
