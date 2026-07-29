# histograma_gates_kiaa0753 (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Histograma de reconstruccion de distribucion (aka Fig 1B) para
# KIAA0753_1: cuentas normalizadas por gate de sorting (FACS), ambas
# replicas diferenciadas por alpha.
#
# Resync 2026-07-29: separado de distribuciones_expresion (que paso a Fig.
# R1.6) en su propia figura - docs/Fig R1.pptx le dio slide propio.
#
# Portado de transcriptional_library/Analysis/scripts/stable_validation.qmd
# (ultimo chunk: gates_df %>% ... %>% ggplot(aes(sample, counts_rel,
# fill=name)) + geom_col() + facet_wrap(~name+rep)). Construido enteramente
# desde data/processed/data_long.tsv (counts_norm por gate/replica/seq_id,
# ver R/01_activity_stats/build_activity_stats.R) - no requiere ninguna
# excepcion, a diferencia del gates_counts.tsv del original.
#
# Requiere: data/processed/data_long.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "histograma_gates_kiaa0753"
out_dir <- fig_dir(slug)

gates_hist <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE) %>%
  filter(name == "KIAA0753_1") %>%
  group_by(rep, name) %>%
  mutate(counts_rel = counts_norm / sum(counts_norm)) %>%
  ungroup()

p <- gates_hist %>%
  ggplot(aes(factor(sample), counts_rel, alpha = rep, group = rep)) +
  geom_col(fill = "#AD343E", position = "dodge") +
  scale_alpha_manual(values = c("Rep 1" = 0.5, "Rep 2" = 1)) +
  labs(x = "Gate de fluorescencia EGFP", y = "Cuentas relativas", alpha = "Réplica", title = "KIAA0753_1") +
  theme_pubclean()
ggsave(file.path(out_dir, "histograma_gates_KIAA0753.jpg"), p, width = 9, height = 6.75, units = "in")

message("Figura guardada en ", out_dir)
