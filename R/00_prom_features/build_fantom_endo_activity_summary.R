# Shared intermediate: actividad CAGE endogena (FANTOM5) resumida por
# promotor, a traves de TODAS las muestras de tejido + celula primaria
# (no un solo tipo celular como build_endo_cage_activity.R) - mean_tpm
# (promedio) y max_tpm (maximo), usados por R4.2/R4.3/R4.4.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("Tissue specificity" section, fantom_mean_expresion/fantom_max_expresion).
#
# EXCEPTION: lee transcriptional_library/Analysis/Tables/{tissue,
# primary_cell}_CAGE_activity.tsv - la misma excepcion ya autorizada
# para sample_CAGE_activity.tsv (especificidad tisular en R7), ver
# R/00_prom_features/analysis_tables_exceptions.R. PENDING: reconstruir
# desde datos crudos.
#
# Run from the TesisDoc repo root. Output:
# data/processed/fantom_endo_activity_summary.tsv (seq_id, mean_tpm, max_tpm)

library(tidyverse)
source("R/00_prom_features/analysis_tables_exceptions.R")

tissue <- read_tsv(path_tissue_cage_activity, show_col_types = FALSE) %>% mutate(group = "tissue")
cell <- read_tsv(path_primary_cell_cage_activity, show_col_types = FALSE) %>% mutate(group = "primary_cell")

fantom <- bind_rows(tissue, cell) %>% mutate(tpm = 1e6 * counts / libsize)

fantom_summary <- fantom %>%
  group_by(name) %>%
  summarise(mean_tpm = mean(tpm), max_tpm = max(tpm)) %>%
  rename(seq_id = name)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(fantom_summary, "data/processed/fantom_endo_activity_summary.tsv")
message("Wrote data/processed/fantom_endo_activity_summary.tsv with ", nrow(fantom_summary), " rows")
