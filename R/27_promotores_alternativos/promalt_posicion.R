# promalt_posicion (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Fig. S7B del paper. Diferencias marginales en actividad y ruido segun
# si el promotor secundario esta rio arriba (upstream) o rio abajo
# (downstream) del principal, dentro del mismo gen.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_alt_pairs.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R) y data/raw/library.bed. Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "Alternative promoters", lineas 1634-1674).

library(tidyverse)
library(ggpubr)
library(rtracklayer)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

slug <- "promalt_posicion"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
pairs_df <- read_tsv("data/processed/prom_alt_pairs.tsv", show_col_types = FALSE)
lib <- import.bed("data/raw/library.bed") %>% as_tibble()
lib$name <- str_remove(lib$name, "^FP.{6}_")

data_ranked <- add_noise_rank(stats_highconf)

pair_proms <- pairs_df %>%
  distinct(name, main_name) %>%
  inner_join(data_ranked %>% select(name, rep, mean, var_rank_sw), by = "name") %>%
  inner_join(
    data_ranked %>% select(name, rep, mean, var_rank_sw) %>%
      rename(main_name = name, main_mean = mean, main_var_rank = var_rank_sw),
    by = c("main_name", "rep")
  )

ordered_pair_proms <- pair_proms %>%
  left_join(lib %>% select(name, strand, start), by = "name") %>%
  left_join(lib %>% select(name, start), by = c("main_name" = "name"), suffix = c("", "_main")) %>%
  mutate(order_from_main = case_when(
    strand == "+" & start < start_main ~ "Río arriba",
    strand == "+" & start >= start_main ~ "Río abajo",
    strand != "+" & start > start_main ~ "Río arriba",
    TRUE ~ "Río abajo"
  ))

p_mean <- ggviolin(ordered_pair_proms,
  x = "order_from_main", y = "mean", fill = "order_from_main", draw_quantiles = 0.5,
  add = "median_q1q3", alpha = 0.7, palette = c("#3D518C", "#14AFB2")
)
p_mean <- facet(p_mean, facet.by = "rep", ncol = 2)
p_mean +
  stat_compare_means(comparisons = list(c("Río arriba", "Río abajo")), label = "p.signif", method = "wilcox.test", na.rm = TRUE) +
  theme_pubclean() +
  labs(y = "Actividad media transcripcional", x = "Posición relativa al promotor principal", fill = NULL) +
  theme(text = element_text(size = 15)) +
  ylim(0, 7.5)
ggsave(file.path(out_dir, "promalt_posicion_actividad.jpg"), width = 12, height = 6.75, units = "in")

p_noise <- ggviolin(ordered_pair_proms,
  x = "order_from_main", y = "var_rank_sw", fill = "order_from_main", draw_quantiles = 0.5,
  add = "median_q1q3", alpha = 0.7, palette = c("#3D518C", "#14AFB2")
)
p_noise <- facet(p_noise, facet.by = "rep", ncol = 2)
p_noise +
  stat_compare_means(comparisons = list(c("Río arriba", "Río abajo")), label = "p.signif", method = "wilcox.test", na.rm = TRUE) +
  theme_pubclean() +
  labs(y = "Ruido transcripcional (rango de varianza)", x = "Posición relativa al promotor principal", fill = NULL) +
  theme(text = element_text(size = 15))
ggsave(file.path(out_dir, "promalt_posicion_ruido.jpg"), width = 12, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
