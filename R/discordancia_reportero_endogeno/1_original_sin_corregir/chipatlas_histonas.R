# discordancia_reportero_endogeno / 1_original_sin_corregir / histonas
#
# VERSION ORIGINAL (sin corregir) del analisis de "diferencia de rango"
# aplicado a marcas de histona (ChIP-Atlas). Se mantiene aca (2026-07-30)
# como referencia de comparacion junto con remap_tf.R (mismo diagnostico
# de confusion aplica). No escribe via fig_dir()/mapping_figuras.csv,
# ver README.md de la carpeta padre.
#
# Analogo al panel de TFs de ReMap en remap_tf.R, pero con la union de
# cada marca de histona (ChIP-Atlas, union en cualquier experimento)
# como feature: para cada marca, dif_rank_rank ~ rank_reporter_scaled +
# marca + (1|replica), correccion BH sobre el efecto fijo. Sin GSEA (a
# diferencia del panel de TFs) - mismo criterio que R3.4 (chipatlas_
# histonas) frente a R3.3 (remap_tf_noise_gsea): no hay una ontologia
# funcional sensata para enriquecer sobre marcas de histona.
#
# Requiere: data/external/allPeaks_chipatlas_counted.tsv (ver
# R/00_prom_features/build_chipatlas_hits.R), data/processed/
# activity_stats_highconf.tsv, data/processed/prom_df.tsv y
# data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

out_dir <- "figures/discordancia_reportero_endogeno/1_original_sin_corregir"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data prep: diferencia de rango x union a cada marca de histona -------

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  mutate(max_tpm = replace_na(max_tpm, 0)) %>%
  group_by(rep) %>%
  mutate(
    rank_reporter = dense_rank(mean),
    rank_endo = row_number(max_tpm),
    dif_rank_rank = rank(abs(rank_endo - rank_reporter)),
    rank_reporter_scaled = scale(rank_reporter)[, 1]
  ) %>%
  ungroup()

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
  select(seq_id, rep, dif_rank_rank, rank_reporter_scaled) %>%
  left_join(binary_hist, by = "seq_id") %>%
  mutate(across(starts_with("Hist_"), ~ replace_na(.x, FALSE)))

# Marcas con >100 promotores marcados en AMBAS replicas.
nHist <- hist_data %>%
  group_by(rep) %>%
  summarise(across(starts_with("Hist_"), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(starts_with("Hist_")) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  distinct(name)

hist_data <- hist_data %>%
  mutate(
    rep = as.factor(rep),
    across(all_of(nHist$name), ~ factor(.x, levels = c("TRUE", "FALSE")))
  )

# --- LMM: efecto de cada marca sobre la diferencia de rango ---------------

resultados_lmm_hist <- map(nHist$name, function(var) {
  formula <- as.formula(paste(
    "dif_rank_rank ~ rank_reporter_scaled +",
    paste0("`", var, "`"), "+ (1 | rep)"
  ))
  fit <- lmerTest::lmer(formula, data = hist_data, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = str_remove(var, "^Hist_"))
}) %>%
  list_rbind() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(resultados_lmm_hist, file.path(out_dir, "chipatlas_histonas_rank_dif_lmm.tsv"))

# --- Plot: efecto de cada marca sobre la diferencia de rango --------------

p <- resultados_lmm_hist %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "FDR < 0.001",
      p_adj < 0.01 ~ "FDR < 0.01",
      p_adj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "ns"
    ),
    sig = factor(sig, levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns")),
    variable = fct_reorder(variable, estimate)
  ) %>%
  ggplot(aes(x = estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, alpha = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "FDR < 0.001" = "#c0392b",
    "FDR < 0.01" = "#e67e22",
    "FDR < 0.05" = "#2980b9",
    "ns" = "grey70"
  )) +
  labs(
    title = "[ORIGINAL, SIN CORREGIR] Efecto de marcas de histona (ChIP-Atlas) sobre la\ndiferencia de rango entre actividad endógena y del reportero",
    subtitle = "Diferencia de rango ~ marca de histona + actividad del reportero + (1|réplica)",
    x = "Estimado (β) con IC 95%", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "chipatlas_histonas_rank_dif_lmm.jpg"), p, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
