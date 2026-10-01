# discordancia_reportero_endogeno / 3_doble_glm_dispersion / prom_df
#
# Tercera variante del analisis de discordancia reportero/endogeno
# (2026-07-30, junto con 1_original_sin_corregir y 2_dif_signed_bland_
# altman en esta misma carpeta de staging - no escribe via fig_dir()/
# mapping_figuras.csv, ver README.md de la carpeta padre).
#
# A diferencia de 2_dif_signed_bland_altman (que solo modela la MEDIA de
# dif_signed), este script usa un modelo doble-GLM (glmmTMB con
# dispformula) que ajusta simultaneamente:
#   - la MEDIA de dif_signed  (~ avg_rank + rep + feature)   - direccion
#     del sesgo reportero-vs-endogeno, igual que en 2_dif_signed_
#     bland_altman (deberia dar resultados muy similares, es basicamente
#     el mismo submodelo).
#   - la DISPERSION de dif_signed (~ avg_rank + feature)     - "cuanta
#     incertidumbre/inconsistencia genera esta feature, en cualquier
#     sentido", la pregunta que motivo probar este metodo (ver hilo de
#     la sesion 2026-07-30): ni un control lineal/spline/residuo/
#     matching exacto por rank_reporter, ni avg_rank + profundidad de
#     lectura, lograron desconfundir la MAGNITUD de la discordancia
#     (|dif_signed|) del efecto de actividad cruda (r~-0.85 con todos).
#     Este modelo de dispersion (mas riguroso que una regresion cruda
#     sobre |dif_signed|) da el mismo resultado (r~-0.8 en el analogo de
#     TFs, ver remap_tf.R) - confirma que no es un problema de metodo,
#     sino una propiedad real de los datos.
#
# El valor de este analisis no esta en "arreglar" la magnitud (no se
# pudo) sino en encontrar EXCEPCIONES a la tendencia general actividad-
# vs-dispersion: features cuyo efecto de dispersion se aparta mucho de
# lo que su efecto de actividad predeciria (residuo de dispersion~
# actividad). TATA-box parece ser la mayor excepcion (mas discordante de
# lo esperado por su actividad) - consistente con R3.1 (ruido_cgi_tata,
# TATA predice mas ruido/varianza con un diseño totalmente distinto).
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv, data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R) y data/processed/
# activity_summary.tsv (ver R/15_summary_features/summary_features.R,
# efecto de cada feature sobre actividad cruda - R7). Run from the
# TesisDoc repo root.

library(tidyverse)
library(glmmTMB)
library(ggrepel)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/3_doble_glm_dispersion"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance()

tidy_data <- build_tidy_features(data, keep = c("seq_id", "dif_signed", "avg_rank")) %>%
  mutate(rep = as.factor(rep))
vars <- tidy_data %>% select(-seq_id, -rep, -dif_signed, -avg_rank) %>% names()

# --- Doble-GLM: media y dispersion de dif_signed, por feature -------------

fit_double_glm <- function(var) {
  df <- tidy_data
  df$feat_test <- df[[var]]
  fit <- tryCatch(
    glmmTMB(dif_signed ~ avg_rank + rep + feat_test, dispformula = ~ avg_rank + feat_test, data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(variable = var, mean_estimate = NA_real_, mean_pval = NA_real_, disp_estimate = NA_real_, disp_pval = NA_real_))
  }
  s <- summary(fit)
  cc_mean <- s$coefficients$cond
  cc_disp <- s$coefficients$disp
  tibble(
    variable = var,
    mean_estimate = cc_mean["feat_testFALSE", "Estimate"], mean_pval = cc_mean["feat_testFALSE", "Pr(>|z|)"],
    disp_estimate = cc_disp["feat_testFALSE", "Estimate"], disp_pval = cc_disp["feat_testFALSE", "Pr(>|z|)"]
  )
}

message("Ajustando doble-GLM para ", length(vars), " features...")
resultados <- map(vars, fit_double_glm) %>%
  list_rbind() %>%
  mutate(mean_padj = p.adjust(mean_pval, method = "BH"), disp_padj = p.adjust(disp_pval, method = "BH"))
write_tsv(resultados, file.path(out_dir, "prom_df_doble_glm.tsv"))

# --- Plot: efecto sobre la MEDIA (direccion del sesgo) --------------------

p_mean <- resultados %>%
  mutate(
    sig = case_when(mean_padj < 0.001 ~ "FDR < 0.001", mean_padj < 0.01 ~ "FDR < 0.01", mean_padj < 0.05 ~ "FDR < 0.05", TRUE ~ "ns") %>%
      factor(levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns")),
    variable = fct_reorder(variable, mean_estimate)
  ) %>%
  ggplot(aes(x = mean_estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  scale_color_manual(values = c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")) +
  labs(title = "Doble-GLM: efecto sobre la MEDIA de la discordancia (dirección del sesgo)", x = "Estimado (β) — media", y = NULL, color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "prom_df_media.jpg"), p_mean, width = 9, height = 6.75, units = "in")

# --- Plot: efecto sobre la DISPERSION (incertidumbre, sin importar signo) -

p_disp <- resultados %>%
  mutate(
    sig = case_when(disp_padj < 0.001 ~ "FDR < 0.001", disp_padj < 0.01 ~ "FDR < 0.01", disp_padj < 0.05 ~ "FDR < 0.05", TRUE ~ "ns") %>%
      factor(levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns")),
    variable = fct_reorder(variable, disp_estimate)
  ) %>%
  ggplot(aes(x = disp_estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  scale_color_manual(values = c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")) +
  labs(title = "Doble-GLM: efecto sobre la DISPERSIÓN de la discordancia (incertidumbre)", x = "Estimado (β, escala log) — dispersión", y = NULL, color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "prom_df_dispersion.jpg"), p_disp, width = 9, height = 6.75, units = "in")

# --- Comparacion: dispersion vs. efecto sobre actividad cruda (R7), con
# las 5 mayores excepciones a la tendencia general etiquetadas ----------

activity_summary <- read_tsv("data/processed/activity_summary.tsv", show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(activity_estimate = mean(estimate), .groups = "drop")

comparacion <- resultados %>% inner_join(activity_summary, by = c("variable" = "feature"))
m <- lm(disp_estimate ~ activity_estimate, data = comparacion)
r_val <- cor(comparacion$disp_estimate, comparacion$activity_estimate)
comparacion <- comparacion %>%
  mutate(
    resid = residuals(m),
    rank_resid = rank(-abs(resid)),
    destacado = rank_resid <= 5,
    sig = ifelse(disp_padj < 0.05, "Dispersión: significativo (BH)", "Dispersión: ns")
  )
write_tsv(comparacion, file.path(out_dir, "prom_df_dispersion_vs_actividad.tsv"))

p_comp <- comparacion %>%
  ggplot(aes(activity_estimate, disp_estimate)) +
  geom_smooth(method = "lm", col = "grey50", se = TRUE, linetype = "dashed", linewidth = 0.6) +
  geom_point(aes(color = sig, size = destacado), alpha = 0.8) +
  geom_text_repel(
    data = comparacion %>% filter(destacado),
    aes(label = variable), size = 4, fontface = "bold", color = "#AD343E",
    box.padding = 0.6, max.overlaps = Inf, min.segment.length = 0
  ) +
  scale_color_manual(values = c("Dispersión: significativo (BH)" = "#358AAA", "Dispersión: ns" = "grey70")) +
  scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.8), guide = "none") +
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, label = paste0("r = ", round(r_val, 2)), size = 6, label.size = 0) +
  labs(
    x = "Efecto sobre la actividad cruda (Wilcoxon, R7)",
    y = "Efecto sobre la dispersión de la discordancia\n(doble-GLM, controlando avg_rank)",
    color = NULL,
    title = "Features curadas: actividad vs. dispersión reportero/endógeno",
    subtitle = "Etiquetadas las 5 mayores excepciones a la tendencia general (mayor |residuo|)"
  ) +
  theme_bw(base_size = 14) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "prom_df_dispersion_vs_actividad.jpg"), p_comp, width = 10, height = 7.5, units = "in")

message("r (dispersion vs actividad): ", round(r_val, 3))
message("Figuras guardadas en ", out_dir)
