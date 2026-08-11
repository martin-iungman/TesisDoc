# promalt_pairwise_dif (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Fig. S7A del paper. Los promotores alternativos de un mismo gen no son
# mas similares entre si que pares al azar: densidad de la diferencia
# absoluta de actividad (arriba) y de ruido (abajo) entre pares
# Secundario-Principal observados, contra una distribucion de diferencias
# permutadas (recombinando al azar los pares, n=10 permutaciones).
#
# La leyenda del paper dice "differences between the RANKS of activity",
# pero el codigo original grafica la diferencia de la media CRUDA
# (true_dif_mean, no true_dif_mean_rank - esta ultima se calcula pero
# nunca se usa). A pedido del autor se generan ambas versiones para
# decidir cual usar: pairwise_dif_promalt_actividad_cruda.jpg (lo que
# realmente hace el codigo) y pairwise_dif_promalt_actividad_rank.jpg
# (lo que dice la leyenda). Para ruido no hay ambiguedad - el original
# solo usa la version de rank (true_dif_var_rank).
#
# Requiere: data/processed/activity_stats_highconf.tsv y data/processed/
# prom_alt_pairs.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "Alternative promoters", lineas 1560-1632).

library(tidyverse)
library(ggpubr)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

set.seed(1)

slug <- "promalt_pairwise_dif"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
pairs_df <- read_tsv("data/processed/prom_alt_pairs.tsv", show_col_types = FALSE)

data_ranked <- add_noise_rank(stats_highconf)

pair_proms <- pairs_df %>%
  distinct(name, main_name) %>%
  inner_join(data_ranked %>% select(name, rep, mean_rank_sw, mean, var_rank_sw), by = "name") %>%
  inner_join(
    data_ranked %>% select(name, rep, mean_rank_sw, mean, var_rank_sw) %>%
      rename(main_name = name, main_mean_rank = mean_rank_sw, main_mean = mean, main_var_rank = var_rank_sw),
    by = c("main_name", "rep")
  ) %>%
  mutate(
    true_dif_mean_rank = abs(mean_rank_sw - main_mean_rank),
    true_dif_mean = abs(mean - main_mean),
    true_dif_var_rank = abs(var_rank_sw - main_var_rank)
  ) %>%
  group_split(rep)

# n=10 permutaciones (reordenamientos completos, sin reemplazo) del valor
# del secundario contra el principal fijo - da una muestra "Permutado" 10x
# mas grande que la observada, para una densidad mas suave.
permute_diff <- function(pair_proms, main_col, alt_col, n = 10) {
  map(pair_proms, function(df) {
    map(seq_len(n), ~ abs(df[[main_col]] - sample(df[[alt_col]]))) %>%
      set_names(paste0("perm", seq_len(n))) %>%
      as_tibble() %>%
      pivot_longer(everything(), values_to = "value") %>%
      mutate(dif = "Permutado", rep = unique(df$rep))
  })
}

pairwise_density <- function(true_dif_col, main_col, alt_col, titulo, xlab) {
  permuted <- permute_diff(pair_proms, main_col, alt_col)
  observed <- map(pair_proms, ~ tibble(value = .x[[true_dif_col]], dif = "Observado", rep = unique(.x$rep)))
  bind_rows(observed, permuted) %>%
    ggplot(aes(value, fill = dif, linetype = rep)) +
    geom_density(alpha = 0.5) +
    labs(x = xlab, fill = "", linetype = "Réplica") +
    scale_fill_manual(values = c("Observado" = "#3D518C", "Permutado" = "#14AFB2")) +
    theme_pubclean() +
    ggtitle(titulo)
}

p_mean_crudo <- pairwise_density("true_dif_mean", "main_mean", "mean",
  "Similitud de actividad entre promotores alternativos", "Diferencia de actividad media (par Secundario-Principal)"
)
ggsave(file.path(out_dir, "pairwise_dif_promalt_actividad_cruda.jpg"), p_mean_crudo, width = 10, height = 10, units = "in")

p_mean_rank <- pairwise_density("true_dif_mean_rank", "main_mean_rank", "mean_rank_sw",
  "Similitud de actividad entre promotores alternativos (rank)", "Diferencia de rango de actividad (par Secundario-Principal)"
)
ggsave(file.path(out_dir, "pairwise_dif_promalt_actividad_rank.jpg"), p_mean_rank, width = 10, height = 10, units = "in")

p_ruido <- pairwise_density("true_dif_var_rank", "main_var_rank", "var_rank_sw",
  "Similitud de ruido entre promotores alternativos", "Diferencia de rango de varianza (par Secundario-Principal)"
)
ggsave(file.path(out_dir, "pairwise_dif_promalt_ruido.jpg"), p_ruido, width = 10, height = 10, units = "in")

message("Figuras guardadas en ", out_dir)
