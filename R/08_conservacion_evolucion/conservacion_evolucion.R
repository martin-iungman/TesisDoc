# conservacion_evolucion_library (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Panel A: PhyloP score a lo largo del metapromotor
# Panel B: clasificacion de conservacion funcional humano-raton (Young et al. 2015)
#
# Requiere: data/processed/prom_df.tsv y data/processed/phylo100_1pb_lib.tsv
# (generados por R/00_prom_features/build_prom_features.R). Run from the
# TesisDoc repo root.

library(tidyverse)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "conservacion_evolucion_library"

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
phylo1pb_prom_df <- read_tsv("data/processed/phylo100_1pb_lib.tsv", show_col_types = FALSE)

# --- Panel A: PhyloP score por posicion relativa al TSS -------------------

panel_a <- phylo1pb_prom_df %>%
  group_by(relpos) %>%
  summarise(med_score = median(score, na.rm = TRUE)) %>%
  ggplot(aes(16 - relpos, med_score)) +
  geom_line(linewidth = 0.4, col = thesis_clr) +
  geom_vline(xintercept = c(-50, -150), linetype = "dashed", linewidth = 0.3) +
  xlab("Posición relativa al TSS") +
  ylab("Mediana del PhyloP score") +
  theme_pubclean(base_size = 16)

# --- Panel B: conservacion funcional humano-raton (Young 2015) -----------

turnover_labels <- tibble(
  cambio = c("Ambigua", "Cambio\nde tejido", "Inserción\nen humanos",
             "Consistente", "Deleción\nen ratón", "Actividad\ndisminuida en ratón"),
  turnover = c("ambiguous", "expression-turnover", "human-inserted",
               "matched", "mouse-deleted", "mouse-diminished")
)

turnover_counts <- prom_df %>%
  filter(type == "promoter", !is.na(turnover)) %>%
  left_join(turnover_labels, by = "turnover") %>%
  count(cambio) %>%
  mutate(total_n = sum(n)) %>%
  arrange(n) %>%
  mutate(cambio = fct_inorder(cambio))

panel_b <- turnover_counts %>%
  count_bar_labeled(cambio, n, total_n, label_nudge_frac = 0.04, coord_flip = TRUE) +
  labs(y = "Frecuencia", x = "Conservación funcional ratón-humano")

# --- Save individual panels + combined figure -----------------------------

out_dir <- fig_dir(slug)
ggsave(file.path(out_dir, "panel_a_phylop.jpg"), panel_a, width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "panel_b_turnover.jpg"), panel_b, width = 9, height = 6.75, units = "in")

combined <- panel_a / panel_b + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, paste0(slug, ".jpg")), combined, width = 9, height = 12, units = "in")
