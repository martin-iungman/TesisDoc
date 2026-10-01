# discordancia_reportero_endogeno / 1_original_sin_corregir / prom_df
#
# VERSION ORIGINAL (sin corregir) del analisis de "diferencia de rango"
# entre actividad del reportero y actividad endogena maxima (FANTOM5),
# portado de transcriptional_library/Analysis/scripts/rank_dif_lmm.R.
# Se mantiene aca (2026-07-30, junto con 2_dif_signed_bland_altman y
# 3_doble_glm_dispersion) como referencia de comparacion, PORQUE ESTE
# DISEÑO ESTA CONFUNDIDO: controla por rank_reporter_scaled unicamente,
# lo cual deja pasar el efecto de cada feature sobre rank_endo (nunca
# modelado) directamente al coeficiente de la feature. Ver
# remap_tf.R en esta misma carpeta para el diagnostico cuantitativo
# (correlacion r~0.8 entre el efecto sobre "diferencia de rango" y el
# efecto sobre actividad cruda de R2.6/R7 - el analisis practicamente no
# aporta nada mas alla de re-detectar que features predicen actividad).
# No escribe via fig_dir()/mapping_figuras.csv (fuera del pipeline
# numerado mientras se decide que version usar), ver README.md de la
# carpeta padre.
#
# Modelo mixto (LMM) del efecto de cada feature booleana del promotor
# sobre la diferencia de rango (magnitud, sin signo dentro del calculo
# pero con signo perdido al usar dense_rank/row_number crudos - ver
# comentario en dif_rank_rank) entre actividad del reportero y actividad
# endogena maxima (FANTOM5): dif_rank_rank ~ rank_reporter_scaled +
# feature + (1|replica), correccion BH por feature, con un control
# negativo (alta actividad del reportero) para validar sensibilidad del
# modelo.
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
library(ggrepel)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/1_original_sin_corregir"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
    rank_reporter_scaled = scale(rank_reporter)[, 1],
    `Control (alta actividad del reportero)` = rank_reporter > median(rank_reporter)
  ) %>%
  ungroup()

# --- Tabla de features booleanas (todas en espanol), + el control -------

tidy_data <- build_tidy_features(
  data,
  keep = c("seq_id", "dif_rank_rank", "rank_reporter_scaled", "Control (alta actividad del reportero)")
)
vars_dicotomicas <- tidy_data %>% select(-seq_id, -rep, -dif_rank_rank, -rank_reporter_scaled) %>% names()

# --- LMM: efecto de cada feature sobre la diferencia de rango -------------

resultados_lmm <- map(vars_dicotomicas, function(var) {
  formula <- as.formula(paste(
    "dif_rank_rank ~ rank_reporter_scaled +",
    paste0("`", var, "`"), "+ (1 | rep)"
  ))
  fit <- lmerTest::lmer(formula, data = tidy_data, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = var)
}) %>%
  list_rbind() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

control_check <- resultados_lmm %>%
  filter(variable == "Control (alta actividad del reportero)") %>%
  select(variable, estimate, conf.low, conf.high, p_value, p_adj)
message("Control (alta actividad del reportero): estimate=", round(control_check$estimate, 3), ", p_adj=", signif(control_check$p_adj, 3))

write_tsv(resultados_lmm, file.path(out_dir, "rank_dif_lmm.tsv"))

# --- Grafico: estimado (beta) +- IC 95%, coloreado por significancia (BH) -

p <- resultados_lmm %>%
  filter(variable != "Control (alta actividad del reportero)") %>%
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
    title = "[ORIGINAL, SIN CORREGIR] Efecto de features del promotor sobre la\ndiferencia de rango entre actividad endógena y del reportero",
    subtitle = "Diferencia de rango ~ feature del promotor + actividad del reportero + (1|réplica)",
    x = "Estimado (β) con IC 95%", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "rank_dif_lmm.jpg"), p, width = 9, height = 6.75, units = "in")

# --- Comparacion con el efecto sobre actividad cruda (R7) ------------------
# Mismo diagnostico que en remap_tf.R (misma carpeta): si el analisis de
# "diferencia de rango" aportara algo mas alla de actividad cruda, no
# deberia correlacionar fuerte con el efecto de Wilcoxon de R7
# (summary_features.R, data/processed/activity_summary.tsv).

activity_summary <- read_tsv("data/processed/activity_summary.tsv", show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(activity_estimate = mean(estimate), .groups = "drop")

comparacion <- resultados_lmm %>%
  filter(variable != "Control (alta actividad del reportero)") %>%
  inner_join(activity_summary, by = c("variable" = "feature"))
r_val <- cor(comparacion$activity_estimate, comparacion$estimate)

# Outliers: features cuyo efecto sobre "diferencia de rango" se aparta
# mas de lo que predice la tendencia general con actividad cruda (mismo
# criterio de residuo que en 3_doble_glm_dispersion/prom_df_features.R) -
# son las mas interesantes, ya que el resto es en gran parte redundante
# con actividad.
m <- lm(estimate ~ activity_estimate, data = comparacion)
comparacion <- comparacion %>%
  mutate(resid = residuals(m), rank_resid = rank(-abs(resid)), destacado = rank_resid <= 5)
write_tsv(comparacion, file.path(out_dir, "rank_dif_lmm_vs_activity.tsv"))

p_comparacion <- comparacion %>%
  ggplot(aes(activity_estimate, estimate)) +
  geom_smooth(method = "lm", col = "grey50", se = TRUE, linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(size = destacado), col = "#358AAA", alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = comparacion %>% filter(destacado),
    aes(label = variable), size = 4, fontface = "bold", color = "#AD343E",
    box.padding = 0.6, max.overlaps = Inf, min.segment.length = 0
  ) +
  scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.8), guide = "none") +
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, label = paste0("r = ", round(r_val, 2)), size = 6, label.size = 0) +
  labs(
    x = "Efecto sobre la actividad cruda (Wilcoxon, R7)",
    y = "Efecto sobre la diferencia de rango (LMM, sin corregir)",
    title = "[ORIGINAL, SIN CORREGIR] ¿Aporta algo mas alla de actividad cruda?",
    subtitle = "Etiquetadas las 5 mayores excepciones a la tendencia general (mayor |residuo|)"
  ) +
  theme_bw(base_size = 14)
ggsave(file.path(out_dir, "rank_dif_lmm_vs_activity.jpg"), p_comparacion, width = 10, height = 7.5, units = "in")
message("Correlacion entre efecto de actividad (R7) y efecto de diferencia de rango (SIN CORREGIR): r=", round(r_val, 3))

message("Figuras guardadas en ", out_dir)
