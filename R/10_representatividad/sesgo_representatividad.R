# sesgo_representatividad (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Panel A: venn de comparacion con el reportero de splicing pre-sorting
#   (proyecto 2022) + test de Fisher de asociacion.
# Panel B/C: sesgo de contenido G+C en las lecturas secuenciadas
#   (2 representaciones alternativas: lineal y de barras) y violin de
#   G+C para promotores secuenciados vs no secuenciados.
#
# Requiere: data/processed/prom_df.tsv, data/processed/activity_stats_full.tsv,
# data/processed/data_long.tsv (ver R/00_prom_features y R/01_activity_stats)
# y data/external/2022_presort_counts.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "sesgo_representatividad"
out_dir <- fig_dir(slug)
clrs <- ghibli::ghibli_palette(name = "PonyoMedium")

all_seq <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE) %>% filter(!seq_id %in% dupl_id)
data_long <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
full_stats <- read_tsv("data/processed/activity_stats_full.tsv", show_col_types = FALSE)
spike_id <- c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")

# --- Panel A: venn contra reportero de splicing 2022 (pre-sorting) -------

presort_2022 <- read_tsv("data/external/2022_presort_counts.tsv", show_col_types = FALSE)
venn_sets <- list(
  "Reportero de splicing\n(pre-sorting)" = presort_2022 %>%
    filter(original_count > 0, seqname %in% all_seq$seq_id) %>%
    pull(seqname),
  "Reportero\ntranscripcional" = data_long %>%
    filter(counts > 0, !seq_id %in% spike_id) %>%
    pull(seq_id) %>%
    unique(),
  "Library" = all_seq$seq_id
)

VennDiagram::venn.diagram(
  venn_sets,
  file.path(out_dir, "panel_a_venn_presort2022.png"),
  imagetype = "png",
  disable.logging = TRUE,
  alpha = c(.5, .5, .5),
  fill = c("#AADAD4", "#358AAA", "#dfedf6"),
  cat.cex = c(0, 0, 0)
)

fisher_result <- tibble(
  all = venn_sets$Library,
  splicing = venn_sets$Library %in% venn_sets$`Reportero de splicing\n(pre-sorting)`,
  transcr = venn_sets$Library %in% venn_sets$`Reportero\ntranscripcional`
) %>%
  select(-all) %>%
  table() %>%
  fisher.test()
message("Fisher test (splicing x transcripcional): p = ", format.pval(fisher_result$p.value), ", OR = ", round(fisher_result$estimate, 3))

# --- Panel B: sesgo de G+C en lecturas secuenciadas -----------------------

tmp <- bind_rows(prom_df %>% mutate(rep = "Rep 1"), prom_df %>% mutate(rep = "Rep 2")) %>%
  left_join(full_stats, by = c("seq_id", "rep")) %>%
  filter(!seq_id %in% spike_id) %>%
  mutate(sum_norm_counts = replace_na(sum_norm_counts, 0))
tmp$bin_gc <- as.numeric(as.character(cut_interval(tmp$g_c,
  length = 0.02,
  labels = seq(min(tmp$g_c) - .01, max(tmp$g_c) + .01, by = 0.02)
)))
tmp <- tmp %>%
  group_by(bin_gc, rep) %>%
  summarise(sum = sum(sum_norm_counts), n = n(), .groups = "drop")
tmp <- tmp %>%
  group_by(rep) %>%
  mutate(sum_rel = sum / sum(sum), n_rel = n / sum(n)) %>%
  select(sum_rel, bin_gc, n_rel, rep)

panel_b_lineal <- tmp %>%
  pivot_longer(c(sum_rel, n_rel), values_to = "counts", names_to = "data") %>%
  ggplot(aes(bin_gc, counts, col = data)) +
  geom_line() +
  facet_wrap(~rep) +
  xlab("Contenido de [G+C]") +
  ylab("Proporción de lecturas secuenciadas") +
  scale_color_manual(values = clrs[c(3, 5)], labels = c("Distribución uniforme", "Observada")) +
  theme_pubr() +
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(file.path(out_dir, "panel_b_gc_bias_lineal.jpg"), panel_b_lineal, width = 9, height = 6.75, units = "in")

panel_b_bars <- tmp %>%
  ggplot(aes(x = bin_gc)) +
  geom_col(aes(y = n_rel), fill = "#AADAD4") +
  geom_col(aes(y = sum_rel), alpha = 0.5, fill = "#358AAA") +
  facet_wrap(~rep) +
  xlab("Contenido de [G+C]") +
  ylab("Proporción") +
  theme_pubr() +
  theme(text = element_text(size = 25))
ggsave(file.path(out_dir, "panel_b_gc_bias_bars.jpg"), panel_b_bars, width = 9, height = 6.75, units = "in")

# --- Panel C: violin de G+C, secuenciados vs no --------------------------

panel_c <- prom_df %>%
  select(seq_id, g_c) %>%
  unique() %>%
  left_join(full_stats, by = "seq_id") %>%
  ggplot(aes(!is.na(sum_counts), g_c)) +
  geom_violin(draw_quantiles = 0.5, fill = "#AADAD4") +
  labs(x = "Promotores secuenciados", y = "Contenido de [G+C]") +
  theme_pubr() +
  theme(text = element_text(size = 25))
ggsave(file.path(out_dir, "panel_c_gc_bias_violin.jpg"), panel_c, width = 4.5, height = 6.75, units = "in")

message("Paneles guardados en ", out_dir, ". Panel B tiene 2 versiones alternativas (lineal/bars) - definir cual usar.")
