# actividad_endogena_multi_tejido (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Actividad del promotor endogeno (CAGE TPM, FANTOM5) vs. actividad del
# promotor reportero, por bin de actividad (mean_sw) y replica: 3
# paneles, HEK293, HeLa y musculo esqueletico.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv (hek_tpm, ver R/00_prom_features/build_prom_features.R) y
# data/processed/endo_cage_activity.tsv (hela_tpm, muscle_tpm, ver
# R/00_prom_features/build_endo_cage_activity.R). Run from the TesisDoc
# repo root.
#
# Ported from transcriptional_library/Analysis/scripts/cell_line_cage.qmd
# (HeLa y musculo esqueletico) y hek_cage.qmd (HEK293).

library(tidyverse)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "actividad_endogena_multi_tejido"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
endo_cage <- read_tsv("data/processed/endo_cage_activity.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(endo_cage, by = "seq_id") %>%
  mutate(across(c(hek_tpm, hela_tpm, muscle_tpm), ~ replace_na(.x, 0)))
data <- add_mean_sw_bins(data)

p_hek <- endo_activity_scatter(data, "hek_tpm", "HEK293 cell line")
p_hela <- endo_activity_scatter(data, "hela_tpm", "HeLa cell line")
p_muscle <- endo_activity_scatter(data, "muscle_tpm", "Skeletal muscle")

ggsave(file.path(out_dir, "HEK293_mean.jpg"), p_hek, width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "Hela_mean.jpg"), p_hela, width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "skeletal_muscle_mean.jpg"), p_muscle, width = 9, height = 6.75, units = "in")

combined <- p_hek / p_hela / p_muscle + plot_annotation(title = "Endogenous promoter usage vs. activity")
ggsave(file.path(out_dir, paste0(slug, ".jpg")), combined, width = 9, height = 18, units = "in")

message("Figuras guardadas en ", out_dir)
