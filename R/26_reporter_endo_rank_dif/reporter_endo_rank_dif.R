# reporter_endo_rank_dif (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Diferencia absoluta entre el rango de actividad del reportero y el
# rango de actividad endogena maxima (FANTOM5, todas las muestras),
# separando housekeeping vs. tissue-specific (tercil de sample_
# specificity_gini), contra un control permutado DENTRO de cada grupo.
#
# Los rangos se calculan sobre los dos terciles juntos (por replica). El
# control mezcla rank_endo entre los miembros de cada grupo dejando
# rank_reporter fijo (1000 permutaciones por replica y grupo), usando los
# rangos que el grupo efectivamente tiene: asi preserva donde esta ubicado
# cada grupo en ambos ordenamientos y solo destruye el apareamiento
# promotor a promotor. Un control global (pares de rangos sorteados sobre
# todo el conjunto) no sirve, porque los housekeeping estan concentrados
# arriba en los dos ordenamientos y eso ya achica su diferencia de rango
# por construccion.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv y data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("Tissue specificity" section, fantom_max_expresion / rank_dif, lineas
# 799-848).

library(tidyverse)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

set.seed(1)

slug <- "reporter_endo_rank_dif"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  mutate(max_tpm = replace_na(max_tpm, 0)) %>%
  filter(!is.na(sample_specificity_gini)) %>%
  group_by(rep) %>%
  mutate(specificity_tercile = cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) %>%
  filter(specificity_tercile != 2) %>%
  mutate(
    rank_reporter = row_number(mean),
    rank_endo = row_number(max_tpm),
    patron = ifelse(specificity_tercile == 1, "Housekeeping", "Tejido-específico")
  ) %>%
  ungroup()

observed <- data %>%
  transmute(rep, patron, value = abs(rank_endo - rank_reporter))

n_perm <- 1000

control <- data %>%
  group_by(rep, patron) %>%
  group_modify(~ tibble(value = abs(
    as.vector(replicate(n_perm, sample(.x$rank_endo))) - rep(.x$rank_reporter, n_perm)
  ))) %>%
  ungroup()

plot_data <- bind_rows(
  observed %>% mutate(tipo = "Observado"),
  control %>% mutate(tipo = "Permutado")
) %>%
  mutate(
    patron = factor(patron, levels = c("Housekeeping", "Tejido-específico")),
    tipo = factor(tipo, levels = c("Observado", "Permutado"))
  )

patron_colors <- c("Housekeeping" = "#1B8C8E", "Tejido-específico" = "#0D2C54")

# Relleno solo para las curvas observadas; las permutadas van sin relleno.
p <- ggplot(plot_data, aes(value, col = patron, linetype = tipo, group = interaction(patron, tipo))) +
  geom_density(
    data = ~ filter(.x, tipo == "Observado"), aes(fill = patron),
    alpha = 0.5, linewidth = 0.7, key_glyph = draw_key_path
  ) +
  geom_density(data = ~ filter(.x, tipo == "Permutado"), fill = NA, linewidth = 0.9, key_glyph = draw_key_path) +
  facet_wrap(~rep) +
  scale_color_manual(values = patron_colors) +
  scale_fill_manual(values = patron_colors, guide = "none") +
  scale_linetype_manual(values = c("Observado" = "solid", "Permutado" = "dashed")) +
  guides(
    col = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)),
    linetype = guide_legend(override.aes = list(colour = "black", linewidth = 1, alpha = 1), keywidth = unit(2.5, "lines"))
  ) +
  labs(
    x = "Diferencia absoluta de rango (reportero vs. actividad endógena)",
    y = "Densidad", col = NULL, linetype = NULL
  ) +
  theme_bw(base_size = 20) +
  theme(legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("reporter_endo_rank_dif.", ext)), p, width = 9, height = 6.75, units = "in")
}

message("Figuras guardadas en ", out_dir)
