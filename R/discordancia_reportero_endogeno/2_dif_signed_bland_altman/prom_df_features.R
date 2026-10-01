# discordancia_reportero_endogeno / 2_dif_signed_bland_altman / prom_df
#
# SACADO DEL PIPELINE NUMERADO 2026-07-30 (era R4.6/R/27_rank_dif_lmm) y
# puesto en esta carpeta de staging junto con las otras 2 variantes del
# analisis (1_original_sin_corregir, 3_doble_glm_dispersion) mientras el
# autor decide cual usar para la tesis - no escribe via fig_dir()/
# mapping_figuras.csv, ver README.md de la carpeta padre.
#
# Modelo mixto (LMM) del efecto de cada feature booleana del promotor
# sobre la discordancia CON SIGNO entre el rango de actividad del
# reportero y el rango de actividad endogena maxima (FANTOM5):
# dif_signed ~ avg_rank + feature + (1|replica), con correccion BH sobre
# el efecto fijo de la feature. dif_signed = rr - re (rangos de
# reportero/endogena normalizados 0-1 por replica; convencion elegida
# por el autor 2026-07-31): positivo = el reportero SOBREESTIMA respecto
# a la actividad endogena para promotores con esa feature, negativo =
# el reportero SUBESTIMA.
#
# SIGNO DEL COEFICIENTE POR FEATURE (`estimate`): cada feature es un
# factor(levels=c("TRUE","FALSE")), TRUE=referencia, asi que lmer/broom
# devuelven crudo el contraste FALSE-menos-TRUE (contraintuitivo para
# leer en un grafico rotulado por nombre de feature). Por eso, justo
# despues de armar resultados_lmm, invertimos el signo (y el IC) para que
# quede como el contraste directo TRUE-menos-FALSE: de ahi en mas,
# `estimate` positivo = promotores CON esa feature tienen dif_signed mas
# alto = el reportero SOBREESTIMA; negativo = SUBESTIMA. Verificado
# contra medias crudas por grupo y contra el modelo releveled a mano
# (sesion 2026-07-31). Ejemplo real: "Enhancers a 50kb" da estimate
# negativo (~-0.054) -> promotores con enhancer cercano SUBESTIMAN mas,
# consistente con la hipotesis de dependencia de contexto.
#
# REDISEÑADO 2026-07-29 (discutido con el autor, ver hilo de la sesion).
# La primera version de este script controlaba por rank_reporter_scaled
# unicamente (dif_rank_rank = rank(|rank_endo-rank_reporter|) ~
# rank_reporter_scaled + feature). Eso deja pasar el efecto de cada
# feature sobre rank_endo (nunca modelado) directamente al coeficiente
# de la feature: se confirmo empiricamente que ni un control mas
# flexible (spline df=4 sobre rank_reporter_scaled) ni matching exacto
# por estratos de ~100 promotores de rank_reporter resolvian esto (r
# entre 0.6 y 0.8 entre el efecto sobre dif_rank y el efecto de la misma
# feature sobre actividad cruda, medido en el analogo de TFs de ReMap,
# R/28_remap_rank_dif_lmm_gsea/). La solucion (analisis tipo
# Bland-Altman/mean-difference, add_rank_discordance() en
# R/functions/plot_helpers.R): descomponer (rank_reporter, rank_endo) en
# avg_rank (nivel de actividad, ortogonal por construccion a la
# discordancia - cor~-0.02 en estos datos) y dif_signed (discordancia).
# Controlando por avg_rank en vez de rank_reporter, la correlacion
# residual con el efecto de actividad cruda cae a ~0.27 para la
# discordancia CON SIGNO - la version de MAGNITUD (|dif_signed|) sigue
# fuertemente confundida (r~-0.85) incluso con este control, probable-
# mente porque es una pregunta de varianza/heterocedasticidad y no de
# nivel medio, que el diseño Bland-Altman no resuelve. Por eso este
# script reporta la discordancia CON SIGNO (direccion del sesgo
# reportero-vs-endogeno), no la magnitud como la version original. El
# "control negativo" de sensibilidad de la version anterior (alta
# actividad del reportero) se elimino: con avg_rank como control ese
# control quedaria algebraicamente equivalente a dif_signed (rr = avg_
# rank - dif_signed/2), tautologico - ya no es un chequeo valido.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv y data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/rank_dif_lmm.R -
# solo la parte que produce la figura de la tesis (Plots/Review/
# endo_reporter_rank_dif_lmm.pdf); se omiten los graficos exploratorios
# iniciales (hex plots) y el chequeo de supuestos del modelo al final,
# que no son figuras de la tesis. tidy_data usa build_tidy_features()
# (R/functions/plot_helpers.R, compartido con R7/R3.2) en vez del
# "tidy_data" del original (no definido en ningun script del repo viejo).

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance()

# --- Tabla de features booleanas (todas en espanol) -----------------------

tidy_data <- build_tidy_features(data, keep = c("seq_id", "dif_signed", "avg_rank"))
vars_dicotomicas <- tidy_data %>% select(-seq_id, -rep, -dif_signed, -avg_rank) %>% names()

# --- LMM: efecto de cada feature sobre la discordancia con signo ---------

resultados_lmm <- map(vars_dicotomicas, function(var) {
  formula <- as.formula(paste(
    "dif_signed ~ avg_rank +",
    paste0("`", var, "`"), "+ (1 | rep)"
  ))
  fit <- lmerTest::lmer(formula, data = tidy_data, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = var)
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

stopifnot(round(resultados_lmm$estimate[resultados_lmm$variable == "TATA-box"], 4) == -0.0459)

# --- Extra: promoter orientation (unidirectional vs. bidirectional) -------
# Separate fit, NOT via build_tidy_features()/vars_dicotomicas above: this
# feature comes from PRO-seq (Core & Martins Orientation Index, see
# ../transcriptional_library/Analysis/scripts/promoter_orientation.R), not
# from prom_df's regular columns, and "not_detected" (no PRO-seq signal)
# must be excluded from the tested universe for this feature only, rather
# than folded into TRUE/FALSE like the rest of build_tidy_features() does.
# Read from transcriptional_library/Analysis/Tables/prom_df.tsv (not this
# repo's data/processed/prom_df.tsv, which doesn't have this column yet)
# per author instruction 2026-09-23.

orientation_df <- read_tsv("../transcriptional_library/Analysis/Tables/prom_df.tsv", show_col_types = FALSE) %>%
  select(seq_id, orientation) %>%
  distinct()

tidy_data_orient <- data %>%
  select(seq_id, rep, dif_signed, avg_rank) %>%
  inner_join(orientation_df, by = "seq_id") %>%
  filter(orientation %in% c("unidirectional", "bidirectional")) %>%
  mutate(`Promotor unidireccional` = factor(
    ifelse(orientation == "unidirectional", "TRUE", "FALSE"),
    levels = c("TRUE", "FALSE")
  ))

message("Promotor unidireccional: n = ", nrow(tidy_data_orient), " (excluding not_detected/NA orientation)")

fit_orient <- lmerTest::lmer(dif_signed ~ avg_rank + `Promotor unidireccional` + (1 | rep), data = tidy_data_orient, REML = FALSE)
resultado_orient <- broom.mixed::tidy(fit_orient, effects = "fixed", conf.int = TRUE) %>%
  rename_with(~"p_value", matches("p[._]value")) %>%
  filter(grepl("Promotor unidireccional", term, fixed = TRUE)) %>%
  mutate(variable = "Promotor unidireccional") %>%
  mutate(
    conf.low_true = -conf.high,
    conf.high_true = -conf.low,
    estimate = -estimate,
    conf.low = conf.low_true,
    conf.high = conf.high_true
  ) %>%
  select(-conf.low_true, -conf.high_true)

resultados_lmm <- bind_rows(resultados_lmm, resultado_orient) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(resultados_lmm, file.path(out_dir, "rank_dif_lmm.tsv"))

# --- Grafico: estimado (beta) +- IC 95%, coloreado por significancia (BH) -

p <- resultados_lmm %>%
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
    title = "Efecto de features del promotor sobre la discordancia (con signo)\nentre actividad del reportero y endógena",
    subtitle = "Discordancia (reportero - endo) ~ feature del promotor + nivel de actividad + (1|réplica)",
    x = "Estimado (β) con IC 95% — positivo: el reportero sobreestima", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "rank_dif_lmm.jpg"), p, width = 9, height = 6.75, units = "in")

# --- Comparacion con el efecto sobre actividad cruda (R7) ------------------
# Mismo diagnostico que para TFs (remap_tf.R, misma carpeta): el diseño
# corregido (avg_rank+dif_signed) deberia mostrar mucha menos correlacion
# residual con el efecto de actividad cruda que el diseño original
# (1_original_sin_corregir/prom_df_features.R, que da r=0.285, aunque con
# solo 31 features esa r ya no era significativa - ver discusion de la
# sesion 2026-07-31).

activity_summary <- read_tsv("data/processed/activity_summary.tsv", show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(activity_estimate = mean(estimate), .groups = "drop")

comparacion <- resultados_lmm %>% inner_join(activity_summary, by = c("variable" = "feature"))
r_val <- cor(comparacion$activity_estimate, comparacion$estimate)
write_tsv(comparacion, file.path(out_dir, "rank_dif_lmm_vs_activity.tsv"))

p_comparacion <- comparacion %>%
  ggplot(aes(activity_estimate, estimate)) +
  geom_point(alpha = 0.6, col = "#358AAA") +
  geom_smooth(method = "lm", col = "#216869") +
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, label = paste0("r = ", round(r_val, 2)), size = 6, label.size = 0) +
  labs(
    x = "Efecto sobre la actividad cruda (Wilcoxon, R7)",
    y = "Efecto sobre la discordancia con signo (avg_rank + dif_signed)",
    title = "¿El diseño corregido sigue reflejando actividad cruda? (features curadas)"
  ) +
  theme_bw(base_size = 14)
ggsave(file.path(out_dir, "rank_dif_lmm_vs_activity.jpg"), p_comparacion, width = 9, height = 6.75, units = "in")
message("Correlacion entre efecto de actividad (R7) y efecto de discordancia con signo: r=", round(r_val, 3))

message("Figuras guardadas en ", out_dir)
