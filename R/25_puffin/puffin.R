# puffin_selectivity_predscore (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Validacion contra el modelo de secuencia PUFFIN (Dudnyk et al. 2024):
# 4 figuras, actividad del reportero (bins por rango) vs. dos scores de
# PUFFIN -
#   - selectivity_dudnik: selectividad tisular predicha por Puffin-D
#     sobre el gen endogeno correspondiente (dato externo publicado,
#     tabla suplementaria del paper).
#   - teoTSS_score: score de PUFFIN en la posicion del TSS (pos=16 de
#     252pb) al correr el modelo sobre nuestras propias secuencias de
#     la library (EXCEPCION, ver mas abajo).
# scatter (bins de 100) + boxplot (20 bins de rango) para cada score.
#
# Requiere: data/processed/activity_stats_highconf.tsv,
# data/processed/prom_df.tsv, data/external/dudnyk2024_sup.xlsx,
# data/external/EPD/human38_epdnew.bed. Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "#PUFFIN") y seq_features.qmd (join selectivity_dudnik contra
# EPD). EXCEPCION: teoTSS_score se lee de Analysis/Tables/
# Dudnyk_puffin_prediction_summ.tsv (path_puffin_pred en
# R/00_prom_features/analysis_tables_exceptions.R) - el pipeline que
# genera esa tabla desde cero (correr el modelo puffin.py sobre nuestra
# library) esta portado pero NO ejecutado en
# R/25_puffin/puffin_prediction_model.R, ver ese archivo para el motivo.

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

# GenomicRanges/S4Vectors define S4 generics that shadow dplyr's
# select/filter/rename; pin them to the dplyr versions (same guard as
# R/00_prom_features/build_prom_features.R).
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

source("R/00_prom_features/analysis_tables_exceptions.R")
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "puffin_selectivity_predscore"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

# --- selectivity_dudnik: join externo contra la tabla suplementaria de
# Dudnyk et al. 2024 (Puffin-D), por posicion de TSS contra EPD ------------

epd <- import.bed("data/external/EPD/human38_epdnew.bed")
selectivity_puffin <- readxl::read_xlsx("data/external/dudnyk2024_sup.xlsx", sheet = 5, skip = 1) %>%
  rename(
    start = TSS,
    selectivity_dudnik = `Selectivity (Puffin-D)`,
    motif_selectivity_dudnik = `Motif Selectivity`,
    dispersion_tissue_dudnik = `Dispersion index (log10)`,
    mean_expr_log10_dudnik = `Mean expression (log10)`
  ) %>%
  select(-geneID, -geneName) %>%
  distinct()
selectivity_puffin$end <- selectivity_puffin$start
gr_selectivity_puffin <- GRanges(selectivity_puffin)
selectivity_puffin_lib <- plyranges::join_overlap_intersect_directed(
  promoters(epd, upstream = 10, downstream = 10), gr_selectivity_puffin
) %>%
  values() %>%
  as_tibble() %>%
  add_count(name) %>%
  filter(n == 1) %>%
  select(-c(n, score))
prom_df <- left_join(prom_df, selectivity_puffin_lib, by = "name")

# --- teoTSS_score: prediccion de PUFFIN sobre nuestra library (EXCEPCION) --

df_pred <- read_tsv(path_puffin_pred, show_col_types = FALSE) %>% select(seq_id, teoTSS_score)
prom_df <- left_join(prom_df, df_pred, by = "seq_id")

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))
data_sw <- add_mean_sw_bins(data)

# --- Scatter (bins de 100 por rango), scores contra actividad -------------

puffin_scatter <- function(data_sw, score_col, ylab, titulo) {
  data_sw %>%
    filter(!is.na(.data[[score_col]])) %>%
    group_by(rep, mean_sw) %>%
    summarise(puffin = median(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(mean_sw, puffin)) +
    geom_point(col = "#14AFB2") +
    geom_smooth(col = "#216869", method = "lm") +
    facet_wrap(~rep) +
    labs(x = "Actividad media (bins de 100, ordenados por rango)", y = ylab, title = titulo) +
    theme_bw(base_size = 20)
}

ggsave(
  file.path(out_dir, "PUFFIN_selectivity_scatter.jpg"),
  puffin_scatter(data_sw, "selectivity_dudnik", "Índice de selectividad (Puffin-D)", "Selectividad PUFFIN"),
  width = 9, height = 6.75, units = "in"
)
ggsave(
  file.path(out_dir, "PUFFIN_predscore_scatter.jpg"),
  puffin_scatter(data_sw, "teoTSS_score", "Score de predicción (mediana)", "Predicción PUFFIN sobre el TSS"),
  width = 9, height = 6.75, units = "in"
)

# --- Boxplot (20 bins de rango) --------------------------------------------

puffin_boxplot <- function(data, score_col, ylab, titulo) {
  data %>%
    filter(!is.na(.data[[score_col]])) %>%
    group_by(rep) %>%
    mutate(mean_sw = row_number(mean) %>% cut_number(20) %>% as.numeric() %>% as_factor()) %>%
    ggplot(aes(mean_sw, .data[[score_col]], group = mean_sw)) +
    geom_boxplot(alpha = 0.6, fill = "#216869", outliers = FALSE) +
    facet_wrap(~rep, scale = "free_x") +
    scale_x_discrete(breaks = c(1, 5, 10, 15, 20)) +
    labs(x = "Actividad media (bins de 20, ordenados por rango)", y = ylab, title = titulo) +
    theme_bw(base_size = 20)
}

ggsave(
  file.path(out_dir, "PUFFIN_selectivity_boxplot.jpg"),
  puffin_boxplot(data, "selectivity_dudnik", "Índice de selectividad (Puffin-D)", "Selectividad PUFFIN"),
  width = 9, height = 6.75, units = "in"
)
ggsave(
  file.path(out_dir, "PUFFIN_predscore_boxplot.jpg"),
  puffin_boxplot(data, "teoTSS_score", "Score de predicción", "Predicción PUFFIN sobre el TSS"),
  width = 9, height = 6.75, units = "in"
)

# --- Selectividad PUFFIN por patron de expresion (housekeeping vs.
# tissue-specific, tercil de sample_specificity_gini) - mismo criterio de
# tercil que R4.3/R4.4 (R/25_especificidad_tisular_fantom,
# R/26_reporter_endo_rank_dif) --------------------------------------------

data_selectivity_pattern <- data %>%
  filter(!is.na(selectivity_dudnik), !is.na(sample_specificity_gini)) %>%
  group_by(rep) %>%
  mutate(specificity_tercile = cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) %>%
  ungroup() %>%
  filter(specificity_tercile != 2) %>%
  mutate(patron = ifelse(specificity_tercile == 1, "Housekeeping", "Tissue-specific")) %>%
  group_by(rep, patron) %>%
  mutate(activity_bin = row_number(mean) %>% cut_number(30) %>% as.numeric()) %>%
  group_by(rep, patron, activity_bin) %>%
  summarise(selectivity_dudnik = median(selectivity_dudnik), .groups = "drop")

p_selectivity_pattern <- data_selectivity_pattern %>%
  ggplot(aes(activity_bin, selectivity_dudnik, col = patron, fill = patron)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm") +
  facet_wrap(~rep) +
  scale_color_manual(values = c("Housekeeping" = "#1B8C8E", "Tissue-specific" = "#0D2C54")) +
  scale_fill_manual(values = c("Housekeeping" = "#1B8C8E", "Tissue-specific" = "#0D2C54")) +
  labs(
    x = "Bins de actividad media del promotor reportero",
    y = "Índice de selectividad PUFFIN",
    title = "Selectividad de promotor predicha vs. actividad",
    col = "Patrón de expresión", fill = "Patrón de expresión"
  ) +
  theme_bw(base_size = 20)
ggsave(file.path(out_dir, "PUFFIN_selectivity_by_pattern.jpg"), p_selectivity_pattern, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
