# coocurrencia_motivos (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Co-ocurrencia entre features de promotores: frecuencia de cada feature
# sobre su propio universo no-NA (el universo difiere entre features,
# plot_feature_frequency()), heatmap de coeficiente phi (correlacion de
# Pearson sobre variables binarias - contempla co-presencia y co-ausencia,
# apropiado para prevalencias moderadas) y heatmap de log2(lift)
# (observado/esperado bajo independencia, P(A∩B)/(P(A)*P(B)) - a
# diferencia de phi no tiene techo cuando las prevalencias estan lejos de
# 50/50, revela relaciones fuertes entre features raras que phi subestima;
# pseudocount de Haldane-Anscombe para evitar lift=0 en pares mutuamente
# excluyentes), cada uno con todas las features y solo con las
# estadisticamente significativas (misma direccion de efecto en ambas
# replicas y ambos limites del IC, ver R7/activity_summary.tsv) - 6
# figuras en total. Jaccard y overlap (Szymkiewicz-Simpson) se evaluaron
# tambien (ver plot_cooccurrence_jaccard()/plot_cooccurrence_overlap() en
# R/functions/plot_helpers.R) pero no se incluyen en esta figura final.
#
# build_cooccurrence_features()/plot_cooccurrence_phi()/
# plot_cooccurrence_lift()/plot_feature_frequency() en
# R/functions/plot_helpers.R, compartidos con R5.7 (promalt_coocurrencia).
# Incluye "Promotor unidireccional" (orientacion de PRO-cap, ver Fig. M14
# / R/00_prom_features/build_promoter_orientation.R) como feature mas de
# build_cooccurrence_features() - resync 2026-09-30: dejo de agregarse
# localmente para estar disponible en todos los analisis de features
# (R7/activity_summary, R3.2/ruido_summary, R5.7/promalt), no solo aca.
#
# Requiere: data/processed/prom_df.tsv (ver R/00_prom_features) y
# data/processed/activity_summary.tsv (ver
# R/15_summary_features/summary_features.R). Run from the TesisDoc repo root.

library(tidyverse)
library(fastDummies)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")

slug <- "coocurrencia_motivos"
out_dir <- fig_dir(slug)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
act_summary <- read_tsv("data/processed/activity_summary.tsv", show_col_types = FALSE)
act_signif_features <- act_summary %>%
  mutate(sign = sign(estimate)) %>%
  filter(pval_corr < 0.05) %>%
  count(feature, sign) %>%
  filter(n == 6) %>%
  pull(feature)

tidy_data <- build_cooccurrence_features(prom_df)

save_pair <- function(data_for_matrix, suffix) {
  feat_mat <- data_for_matrix %>%
    select(-seq_id) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  ggsave(file.path(out_dir, paste0("cooccurrence_phi", suffix, ".jpg")), plot_cooccurrence_phi(feat_mat), width = 9, height = 6.75, units = "in")
  ggsave(file.path(out_dir, paste0("cooccurrence_lift", suffix, ".jpg")), plot_cooccurrence_lift(feat_mat), width = 9, height = 6.75, units = "in")
  freq_height <- pmax(4, ncol(feat_mat) * 0.28 + 1.5)
  ggsave(file.path(out_dir, paste0("feature_frequency", suffix, ".jpg")), plot_feature_frequency(feat_mat), width = 8, height = freq_height, units = "in")
}

# --- 6 figuras: todas las features / solo las significativas -------------

save_pair(tidy_data, "")
save_pair(tidy_data %>% select(seq_id, any_of(act_signif_features)), "_signif")

message("Figuras guardadas en ", out_dir)
