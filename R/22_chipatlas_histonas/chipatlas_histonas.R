# chipatlas_histonas (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Efecto de marcas de histona (ChIP-Atlas, union en cualquier
# experimento) sobre el ruido transcripcional: chipatlas_noise.jpg
# (AUC por marca), noise_scatter_Hist_H3K4me3.jpg (ruido vs. actividad
# por presencia de H3K4me3) y noise_scatter_Hist_H3K4me3_nonCGI.jpg
# (lo mismo pero solo dentro de promotores sin isla CpG, con bins de
# actividad recalculados desde cero en ese subconjunto).
#
# Requiere: data/external/allPeaks_chipatlas_counted.tsv (ver
# R/00_prom_features/build_chipatlas_hits.R), data/processed/
# activity_stats_highconf.tsv y data/processed/prom_df.tsv. Run from
# the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("ChIP-Atlas" section) y Histone_chipatlas.qmd. Se corrigieron 2 bugs
# presentes en ambos scripts originales: (1) el archivo crudo de picos
# tiene un prefijo numerico sin sentido antes de cada seq_id que
# arruinaba el count() de n_samples si no se separaba (ver
# R/00_prom_features/build_chipatlas_hits.R); (2)
# `explist[duplicated(explist$X4),]` hace lo opuesto de lo que parece
# buscarse - descarta los mapeos limpios feature->grupo (los que
# aparecen una sola vez tras el distinct() previo) y solo deja los
# ambiguos, vaciando la tabla casi por completo (H3K4me3 solo mapea a
# "Histone", una unica vez, así que quedaba afuera). Reemplazado por
# distinct(X4, .keep_all=TRUE).

library(tidyverse)
library(ggpubr)
library(pROC)
library(ggnewscale)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "chipatlas_histonas"
out_dir <- fig_dir(slug)

# --- Data prep: ruido x union a cada marca de histona ---------------------

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))
data <- add_noise_rank(data)

chipatlas <- read_tsv("data/external/allPeaks_chipatlas_counted.tsv", col_names = c("seq_id", "feature", "n_samples"), show_col_types = FALSE)
explist <- read_tsv(path_chipatlas_explist, col_names = FALSE, show_col_types = FALSE) %>%
  select(X3, X4) %>%
  distinct(X4, .keep_all = TRUE)
chipatlas <- chipatlas %>%
  left_join(explist, by = c("feature" = "X4")) %>%
  rename(group = X3)

histones <- chipatlas %>% filter(group == "Histone")

binary_hist <- histones %>%
  distinct(seq_id, feature) %>%
  mutate(value = TRUE, feature = paste0("Hist_", feature)) %>%
  pivot_wider(names_from = feature, values_from = value, values_fill = FALSE)

hist_data <- data %>%
  left_join(binary_hist, by = "seq_id") %>%
  mutate(across(starts_with("Hist_"), ~ replace_na(.x, FALSE)))

# Marcas con >100 promotores marcados en AMBAS replicas.
nTF <- hist_data %>%
  group_by(rep) %>%
  summarise(across(starts_with("Hist_"), sum, na.rm = TRUE)) %>%
  pivot_longer(starts_with("Hist_")) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  distinct(name)

vbles <- nTF$name
repname <- unique(hist_data$rep)

# --- AUC (IC DeLong) de cada marca prediciendo ruido alto -----------------

auc_df <- map(repname, function(r) {
  df <- hist_data %>% filter(rep == r)
  map(vbles, function(f) {
    roc_obj <- pROC::roc(df[[f]], df[["var_rank_sw"]], ci = TRUE, direction = ">", quiet = TRUE)
    tibble(feature = f, rep = r, AUC = as.numeric(roc_obj$ci[2]), ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]))
  }) %>% list_rbind()
}) %>% list_rbind() %>%
  mutate(feature = str_remove(feature, "^Hist_"))

write_tsv(auc_df, file.path(out_dir, "chipatlas_noise_auc.tsv"))

panel <- plot_noise_auc_summary(auc_df, show_labels = length(unique(auc_df$feature)) <= 40)
ggsave(file.path(out_dir, "chipatlas_noise.jpg"), panel, width = 9, height = 6.75, units = "in")

# --- H3K4me3: scatter de ruido vs. actividad ------------------------------

hist_data_h3k4me3 <- hist_data %>% rename(H3K4me3 = Hist_H3K4me3)

ggsave(
  file.path(out_dir, "noise_scatter_Hist_H3K4me3.jpg"),
  noise_scatter(hist_data_h3k4me3, "H3K4me3", "Ruido vs. actividad, por H3K4me3", "H3K4me3"),
  width = 9, height = 6.75, units = "in"
)

# --- H3K4me3, solo promotores sin isla CpG --------------------------------
# mean_sw/var_rank_sw recalculados desde cero dentro del subconjunto
# no-CGI (no reutiliza los del dataset completo).

hist_data_h3k4me3_nonCGI <- data %>%
  filter(!CGI) %>%
  select(-mean_sw, -var_rank_sw) %>%
  left_join(binary_hist %>% select(seq_id, any_of("Hist_H3K4me3")), by = "seq_id") %>%
  mutate(H3K4me3 = replace_na(Hist_H3K4me3, FALSE)) %>%
  add_noise_rank()

ggsave(
  file.path(out_dir, "noise_scatter_Hist_H3K4me3_nonCGI.jpg"),
  noise_scatter(hist_data_h3k4me3_nonCGI, "H3K4me3", "Ruido de H3K4me3 en promotores sin isla CpG", "H3K4me3"),
  width = 9, height = 6.75, units = "in"
)

message("Figuras guardadas en ", out_dir)
