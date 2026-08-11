# coocurrencia_motivos (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Co-ocurrencia entre features de promotores: heatmaps de coeficiente phi
# (correlacion de Pearson sobre variables binarias) y similitud de
# Jaccard, cada uno con todas las features y solo con las
# estadisticamente significativas (misma direccion de efecto en ambas
# replicas y ambos limites del IC, ver R7/activity_summary.tsv) - 4
# figuras en total.
#
# build_cooccurrence_features()/plot_cooccurrence_phi()/
# plot_cooccurrence_jaccard() en R/functions/plot_helpers.R, compartidos
# con R5.7 (promalt_coocurrencia).
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
  ggsave(file.path(out_dir, paste0("cooccurrence_jaccard", suffix, ".jpg")), plot_cooccurrence_jaccard(feat_mat), width = 9, height = 6.75, units = "in")
}

# --- 4 figuras: todas las features / solo las significativas -------------

save_pair(tidy_data, "")
save_pair(tidy_data %>% select(seq_id, any_of(act_signif_features)), "_signif")

message("Figuras guardadas en ", out_dir)
