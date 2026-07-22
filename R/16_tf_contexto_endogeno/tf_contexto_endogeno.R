# tf_nfya_sp_contexto_endogeno (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Asociacion entre actividad transcripcional (bins de 100 promotores,
# ordenados por actividad media) y union de NFYA/SP1/SP2 segun ChIP-seq
# de ReMap2022: (1) union en cualquier linea celular, (2) union
# restringida a HEK293 para SP1/SP2, (3) union de NFYA consistente en
# las 3 lineas celulares con datos (interseccion K-562/HeLa-S3/Hep-G2,
# con diagrama de Venn de las 3 lineas).
#
# Requiere: data/processed/remap_tf_hits.tsv (ver
# R/00_prom_features/build_remap_tf_hits.R), data/processed/
# activity_stats_highconf.tsv y data/processed/prom_df.tsv (ver
# R/00_prom_features y R/01_activity_stats). Run from the TesisDoc repo
# root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("ReMap dataset", "SP1/2 only hek" y "NFY Remap intersect" - fixed a
# bug there: the TF binding columns were referenced as TF_NFYA/TF_SP1/
# TF_SP2, a name that only existed in a different, Analysis/Tables-based
# object (data_TF_remap.tsv) from a sibling script - here they're just
# NFYA/SP1/SP2, as build_remap_tf_hits.R produces them).

library(tidyverse)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "tf_nfya_sp_contexto_endogeno"
out_dir <- fig_dir(slug)

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

# Bins de 100 promotores segun rango de actividad media (por replica,
# se descarta el resto de la division para que los bins queden parejos).
excess <- data %>% group_by(rep) %>% summarise(excess = n() %% 100, .groups = "drop")
data <- data %>%
  group_split(rep) %>%
  map2(excess$excess, ~ .x %>%
    mutate(mean_rank_sw = row_number(mean)) %>%
    filter(mean_rank_sw > .y) %>%
    mutate(mean_rank_sw = row_number(mean), mean_sw = ceiling(mean_rank_sw / 100)) %>%
    ungroup()) %>%
  list_rbind()

tf_scatter <- function(df, feature, titulo) {
  df %>%
    group_by(rep, mean_sw) %>%
    summarise(prop = sum(.data[[feature]]) / n(), .groups = "drop") %>%
    ggplot(aes(mean_sw, prop)) +
    geom_point(col = thesis_clr) +
    geom_smooth(col = "#216869") +
    facet_wrap(~rep) +
    ylim(0, 1) +
    labs(x = "Actividad media (bins de 100, ordenados por rango)", y = "Proporción de promotores", title = titulo) +
    theme_bw(base_size = 16)
}

# --- Union en cualquier linea celular (NFYA, SP1, SP2) ---------------------

binary_df <- remap_hits %>%
  filter(TF %in% c("NFYA", "SP1", "SP2")) %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, mean, mean_sw) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(any_of(c("NFYA", "SP1", "SP2")), ~ replace_na(.x, FALSE)))

ggsave(file.path(out_dir, "NFYA_scatter.jpg"), tf_scatter(data_TF, "NFYA", "Unión a NFYA"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "SP1_scatter.jpg"), tf_scatter(data_TF, "SP1", "Unión a SP1"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "SP2_scatter.jpg"), tf_scatter(data_TF, "SP2", "Unión a SP2"), width = 9, height = 6.75, units = "in")

# --- Union restringida a HEK293 (SP1, SP2) ---------------------------------

hek_df <- remap_hits %>%
  filter(str_detect(sample, "HEK293"), TF %in% c("SP1", "SP2")) %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF_hek <- data %>%
  select(seq_id, rep, mean, mean_sw) %>%
  left_join(hek_df, by = "seq_id") %>%
  mutate(across(any_of(c("SP1", "SP2")), ~ replace_na(.x, FALSE)))

ggsave(file.path(out_dir, "SP1_scatter_hek.jpg"), tf_scatter(data_TF_hek, "SP1", "Unión a SP1 - HEK293"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "SP2_scatter_hek.jpg"), tf_scatter(data_TF_hek, "SP2", "Unión a SP2 - HEK293"), width = 9, height = 6.75, units = "in")

# --- NFYA: interseccion entre lineas celulares (venn + union consistente) --

nfya_samples <- remap_hits %>%
  filter(TF == "NFYA") %>%
  separate_rows(sample, sep = ",") %>%
  filter(sample != "GM12878")

venn_list <- nfya_samples %>%
  group_by(sample) %>%
  summarise(fps = list(unique(seq_id)), .groups = "drop") %>%
  deframe()

venn_fill <- colorRampPalette(c(thesis_clr, "#AADAD4", "#dfedf6"))(length(venn_list))
VennDiagram::venn.diagram(
  venn_list,
  file.path(out_dir, "NFYA_venn.png"),
  imagetype = "png",
  disable.logging = TRUE,
  alpha = rep(0.5, length(venn_list)),
  fill = venn_fill
)

nfya_intersect_proms <- Reduce(intersect, venn_list)
data_TF_nfya_intersect <- data %>%
  select(seq_id, rep, mean, mean_sw) %>%
  mutate(NFYA = seq_id %in% nfya_intersect_proms)

ggsave(
  file.path(out_dir, "NFYA_scatter_intersect.jpg"),
  tf_scatter(data_TF_nfya_intersect, "NFYA", paste0("Unión a NFYA - ", paste(names(venn_list), collapse = ", "))),
  width = 9, height = 6.75, units = "in"
)

message("Figuras guardadas en ", out_dir)
