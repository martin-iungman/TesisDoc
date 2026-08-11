# promalt_clasificacion_secundaria (ver docs/mapping_figuras.csv para
# el numero de figura vigente)
# Fig. S7C del paper. Distribucion de la correlacion de Pearson entre
# la actividad endogena del par Secundario-Principal, para cada clase
# de par (switch/correlated/independent), separado segun si el gen
# tiene o no casos de "switch" (el secundario supera al principal en
# alguna muestra).
#
# Requiere: data/processed/prom_alt_pairs.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "Alternative promoters", lineas 1394-1408).

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "promalt_clasificacion_secundaria"
out_dir <- fig_dir(slug)

pairs_df <- read_tsv("data/processed/prom_alt_pairs.tsv", show_col_types = FALSE)

short_df <- pairs_df %>%
  distinct(name, main_name, cor_pearson, cor_spearman, prom_alt, N_switch)

p <- short_df %>%
  mutate(
    with_switch_cases = ifelse(N_switch > 0, "Con casos de alternancia", "Sin casos de alternancia"),
    prom_alt = recode(prom_alt, independent = "Independiente", correlated = "Correlacionado", switch = "Alternancia") %>%
      factor(levels = c("Independiente", "Correlacionado", "Alternancia"))
  ) %>%
  ggplot(aes(cor_pearson, fill = prom_alt)) +
  facet_wrap(~with_switch_cases, nrow = 2) +
  geom_density(alpha = 0.75) +
  theme_pubclean() +
  labs(x = "Correlación de Pearson entre la actividad\ndel par Secundario-Principal", y = "Densidad", fill = "Clase de\npromotor secundario") +
  scale_fill_manual(values = c("#48d6d9", "#14AFB2", "#216869"))
ggsave(file.path(out_dir, "secondary_classification.jpg"), p, width = 9, height = 6.75, units = "in")

message("Figura guardada en ", out_dir)
