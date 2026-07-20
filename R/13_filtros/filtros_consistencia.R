# filtros_consistencia_replicas (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Panel A: correlacion de actividad media entre replicas, antes de
#   excluir promotores bimodales/inconsistentes (coloreado por si el
#   promotor termina en el subconjunto de alta confianza).
# Panel B: la misma correlacion, ya solo sobre el subconjunto filtrado.
# Panel C (SIN RESOLVER): el CSV menciona "subconjunto final de 12214
#   promotores de alta confianza" - hay dos candidatos igual de
#   plausibles en transcriptional_library/Tesis/tesis.R (n_per_rep.jpg y
#   venn_rep_postfilter.jpg), ambos generados abajo. Definir cual es.
#
# Requiere: data/processed/activity_stats_full.tsv,
# data/processed/activity_stats_highconf.tsv y data/processed/data_long.tsv
# (ver R/01_activity_stats/build_activity_stats.R) y data/raw/seq_data.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggExtra)
library(patchwork)
source("R/functions/fig_paths.R")

slug <- "filtros_consistencia_replicas"
out_dir <- fig_dir(slug)

full_stats <- read_tsv("data/processed/activity_stats_full.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)

# --- Panel A: correlacion entre replicas, pre-filtro (coloreada) ---------

panel_a_base <- full_stats %>%
  filter(sum_norm_counts > 1000) %>%
  select(seq_id, mean, rep) %>%
  pivot_wider(values_from = "mean", names_from = "rep") %>%
  mutate(Confianza = ifelse(seq_id %in% stats_highconf$seq_id, "Alta", "Baja") %>% fct_rev()) %>%
  arrange(desc(Confianza))

p <- panel_a_base %>%
  ggplot(aes(`Rep 1`, `Rep 2`, col = Confianza)) +
  geom_point(size = 0.5, alpha = 0.7) +
  labs(x = "Actividad media (Rep 1)", y = "Actividad media (Rep 2)") +
  scale_color_manual(values = c("Alta" = "#1B8C8E", "Baja" = "#0D2C54")) +
  theme_bw() +
  theme(text = element_text(size = 25)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
p <- ggExtra::ggMarginal(p, type = "density", margins = "both", col = "#0D2C54")
ggsave(file.path(out_dir, "panel_a_correlacion_prefiltro.jpg"), p, width = 9, height = 6.75, units = "in")

cor_prefiltro <- panel_a_base %>% summarise(corre = cor(`Rep 1`, `Rep 2`, use = "complete.obs", method = "pearson"))
message("Pearson pre-filtro: ", round(cor_prefiltro$corre, 3))

# --- Panel B: correlacion entre replicas, post-filtro ----------------------

panel_b_base <- stats_highconf %>%
  filter(sum_norm_counts > 1000) %>%
  select(seq_id, mean, rep) %>%
  pivot_wider(values_from = "mean", names_from = "rep")

p <- panel_b_base %>%
  ggplot(aes(`Rep 1`, `Rep 2`)) +
  geom_point(size = 0.5, alpha = 0.5, col = "#1B8C8E") +
  labs(x = "Actividad media (Rep 1)", y = "Actividad media (Rep 2)") +
  theme_bw() +
  theme(text = element_text(size = 25))
p <- ggExtra::ggMarginal(p, type = "density", margins = "both", col = "#0D2C54")
ggsave(file.path(out_dir, "panel_b_correlacion_postfiltro.jpg"), p, width = 7, height = 6.75, units = "in")

cor_postfiltro <- panel_b_base %>% summarise(corre = cor(`Rep 1`, `Rep 2`, use = "complete.obs", method = "pearson"))
message("Pearson post-filtro: ", round(cor_postfiltro$corre, 3))

# --- Panel C candidato 1: N de promotores detectados por replica ---------

data_long <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE)
n_lib <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE) %>% filter(!seq_id %in% dupl_id) %>% nrow()

count_by_rep_combo <- function(seq_id_rep_df, total_n) {
  tmp <- seq_id_rep_df %>%
    mutate(rep = factor("1|2")) %>%
    bind_rows(seq_id_rep_df) %>%
    distinct(seq_id, rep)
  tmp %>%
    count(seq_id) %>%
    filter(n == 3) %>%
    mutate(rep = factor("1&2")) %>%
    bind_rows(tmp %>% distinct(seq_id, rep)) %>%
    count(rep) %>%
    mutate(total_n = total_n)
}

p1 <- count_by_rep_combo(data_long %>% filter(counts > 0) %>% distinct(seq_id, rep), n_lib) %>%
  ggplot(aes(rep, n, fill = rep, label = n)) +
  geom_col(aes(y = total_n), fill = "#dfedf6") +
  geom_col(aes(y = n), fill = "#358AAA") +
  geom_label(nudge_y = 1800, aes(y = n), fill = "#AADAD4") +
  ggtitle("Totales observadas") +
  xlab("Réplica") +
  ylab("Frecuencia") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15))

p2 <- count_by_rep_combo(stats_highconf %>% distinct(seq_id, rep), n_lib) %>%
  ggplot(aes(rep, n, fill = rep, label = n)) +
  geom_col(aes(y = total_n), fill = "#dfedf6") +
  geom_col(aes(y = n), fill = "#358AAA") +
  geom_label(nudge_y = 1800, aes(y = n), fill = "#AADAD4") +
  ggtitle("Post-filtros") +
  xlab("Réplica") +
  ylab("Frecuencia") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15))

ggsave(file.path(out_dir, "panel_c_candidato1_n_por_replica.jpg"), p1 + p2, width = 9, height = 6.75, units = "in")

# --- Panel C candidato 2: Venn post-filtro ---------------------------------

spike_id <- c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")
all_seq <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE) %>% filter(!seq_id %in% dupl_id)

VennDiagram::venn.diagram(
  list(
    "Rep 1" = stats_highconf %>% filter(rep == "Rep 1", !seq_id %in% spike_id) %>% pull(seq_id) %>% unique(),
    "Rep 2" = stats_highconf %>% filter(rep == "Rep 2", !seq_id %in% spike_id) %>% pull(seq_id) %>% unique(),
    "Library" = all_seq$seq_id
  ),
  file.path(out_dir, "panel_c_candidato2_venn_postfiltro.png"),
  imagetype = "png",
  disable.logging = TRUE,
  alpha = c(.5, .5, .5),
  fill = c("#AADAD4", "#358AAA", "#dfedf6"),
  cat.cex = c(0, 0, 0)
)

message("Paneles guardados en ", out_dir, ". Panel C sin resolver - ver dos candidatos.")
