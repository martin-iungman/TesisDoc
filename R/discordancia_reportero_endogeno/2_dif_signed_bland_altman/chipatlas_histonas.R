# discordancia_reportero_endogeno / 2_dif_signed_bland_altman / histonas
#
# SACADO DEL PIPELINE NUMERADO 2026-07-30 (era R4.7/R/28_remap_rank_dif_
# lmm_gsea) y puesto en esta carpeta de staging junto con las otras 2
# variantes del analisis mientras el autor decide cual usar para la
# tesis - no escribe via fig_dir()/mapping_figuras.csv, ver README.md de
# la carpeta padre.
#
# Analogo al panel de TFs de ReMap en remap_tf.R, pero con la union de cada marca de histona (ChIP-Atlas,
# union en cualquier experimento) como feature: para cada marca,
# dif_signed ~ avg_rank + marca + (1|replica), correccion BH sobre el
# efecto fijo. Sin GSEA (a diferencia del panel de TFs) - mismo criterio
# que R3.4 (chipatlas_histonas) frente a R3.3 (remap_tf_noise_gsea): no
# hay una ontologia funcional sensata para enriquecer sobre marcas de
# histona.
#
# dif_signed = rr - re (rangos de reportero/endogena normalizados 0-1
# por replica; convencion elegida por el autor 2026-07-31): positivo =
# el reportero SOBREESTIMA respecto a la actividad endogena, negativo =
# el reportero SUBESTIMA.
#
# SIGNO DEL COEFICIENTE POR MARCA (`estimate`): cada marca es un
# factor(levels=c("TRUE","FALSE")), TRUE=referencia, asi que lmer/broom
# devuelven crudo el contraste FALSE-menos-TRUE (contraintuitivo para un
# grafico rotulado por nombre de marca). Por eso, justo despues de armar
# resultados_lmm_hist, invertimos el signo (y el IC) para que quede como
# el contraste directo TRUE-menos-FALSE: `estimate` positivo = promotores
# marcados tienen dif_signed mas alto = el reportero SOBREESTIMA para esa
# marca; negativo = SUBESTIMA. Ver remap_tf.R para la verificacion
# completa (sesion 2026-07-31).
#
# avg_rank (promedio de ambos rangos) es el control de nivel de
# actividad - ver R/27_rank_dif_lmm/rank_dif_lmm.R para la discusion
# completa de por que reemplaza a controlar por rank_reporter
# solo (REDISEÑADO 2026-07-29).
#
# Requiere: data/external/allPeaks_chipatlas_counted.tsv (ver
# R/00_prom_features/build_chipatlas_hits.R), data/processed/
# activity_stats_highconf.tsv, data/processed/prom_df.tsv y
# data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.
#
# No portado de un script del repo viejo. Nuevo analisis, combinando
# R/27_rank_dif_lmm/rank_dif_lmm.R (el modelo LMM) con la estructura de
# datos de R/22_chipatlas_histonas/chipatlas_histonas.R (ChIP-Atlas,
# filtro >100 promotores marcados en ambas replicas, bug del prefijo
# numerico y del duplicated() ya corregidos alli via
# R/00_prom_features/build_chipatlas_hits.R).

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Data prep: diferencia de rango x union a cada marca de histona -------

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance()

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
  select(seq_id, rep, dif_signed, avg_rank) %>%
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

# --- LMM: efecto de cada marca sobre la discordancia con signo -----------

resultados_lmm_hist <- map(nHist$name, function(var) {
  formula <- as.formula(paste(
    "dif_signed ~ avg_rank +",
    paste0("`", var, "`"), "+ (1 | rep)"
  ))
  fit <- lmerTest::lmer(formula, data = hist_data, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = str_remove(var, "^Hist_"))
}) %>%
  list_rbind() %>%
  mutate(
    conf.low_true = -conf.high,
    conf.high_true = -conf.low,
    estimate = -estimate,
    conf.low = conf.low_true,
    conf.high = conf.high_true
  ) %>%
  select(-conf.low_true, -conf.high_true) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(resultados_lmm_hist, file.path(out_dir, "chipatlas_histonas_rank_dif_lmm.tsv"))

# --- Plot: efecto de cada marca sobre la discordancia con signo ----------

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
    title = "Efecto de marcas de histona (ChIP-Atlas) sobre la discordancia (con signo)\nentre actividad del reportero y endógena",
    subtitle = "Discordancia (reportero - endo) ~ marca de histona + nivel de actividad + (1|réplica)",
    x = "Estimado (β) con IC 95% — positivo: el reportero sobreestima", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "chipatlas_histonas_rank_dif_lmm.jpg"), p, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
