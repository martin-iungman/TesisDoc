# ruido_summary (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Resumen del efecto de cada feature booleana del promotor sobre el
# ruido transcripcional: AUC de la feature prediciendo var_rank_sw alto
# (rango de varianza dentro de cada bin de actividad media - ver
# R/19_ruido_cgi_tata), separado en features de secuencia (seq) y de
# contexto endogeno (endo) igual que summary_features_secuencia /
# summary_features_endogeno, pero los dos paneles forman una sola figura
# (A = secuencia, B = contexto endogeno).
#
# Significancia: p-valor de Wald usando la varianza de DeLong del AUC
# (z = (AUC-0.5)/SE_delong), corregido por multiples tests con BH -
# resync 2026-09-30, reemplaza al criterio anterior de "IC bootstrap que
# no cruza 0.5". Se verifico que el metodo anterior en realidad ya
# calculaba el IC via DeLong (pROC::ci.auc usa DeLong por default; el
# argumento boot.n se ignoraba silenciosamente sin method="bootstrap"),
# asi que el cambio no es metodologico sino que agrega el p-valor
# explicito para poder corregir por tests multiples, cosa que un IC no
# permite. Se comparo explicitamente contra bootstrap real (boot.n=2000,
# resampleo estratificado) para las 33 features x 2 replicas: acuerdo
# practicamente perfecto (Spearman en -log10(p) = 0.99), con la unica
# discrepancia de significancia post-BH en un caso limite (p~0.05), y el
# bootstrap con piso de resolucion en 1/boot.n=0.0005 (da p=0 exacto en
# 25/66 tests donde el efecto es muy fuerte, cosa que DeLong no tiene).
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y
# R/01_activity_stats). Corre en segundos (DeLong es una formula cerrada,
# no hace falta resamplear). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (bloque de resumen de la seccion de ruido, plot_auc2). build_tidy_
# features()/feature_groups (identicas a las de R/15_summary_features,
# el resumen de actividad) viven en R/functions/plot_helpers.R. Incluye
# "Promotor unidireccional" (ver orientacion_promotor) como cualquier otra feature de
# build_tidy_features(), sin tratamiento especial.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(pROC)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "ruido_summary"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data <- add_noise_rank(data)

tidy_data <- build_tidy_features(data, keep = "var_rank_sw")

# --- AUC (p-valor DeLong, corregido por BH) de cada feature prediciendo --
# --- ruido alto ------------------------------------------------------------

vbles <- tidy_data %>% select(-rep, -var_rank_sw) %>% names()
repname <- unique(tidy_data$rep)

auc_df <- map(repname, function(r) {
  df <- tidy_data %>% filter(rep == r)
  map(vbles, function(f) {
    resp <- df[[f]]
    ok <- !is.na(resp)
    roc_obj <- pROC::roc(resp[ok], df[["var_rank_sw"]][ok], ci = TRUE, direction = ">", quiet = TRUE)
    auc_val <- as.numeric(roc_obj$ci[2])
    se_delong <- sqrt(pROC::var(roc_obj, method = "delong"))
    pval <- 2 * pnorm(-abs((auc_val - 0.5) / se_delong))
    tibble(feature = f, rep = r, AUC = auc_val, ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]), pval = pval)
  }) %>% list_rbind()
}) %>% list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, method = "BH"))

write_tsv(auc_df, file.path(out_dir, "ruido_summary_auc.tsv"))

# --- Grafico resumen --------------------------------------------------
# Features con pval_corr<0.05 Y mismo sentido (ruido alto/bajo) en ambas
# replicas (mismo criterio de consistencia que summary_features_secuencia). plot_noise_auc_
# summary() en plot_helpers.R (compartida con remap_tf_noise_gsea, la version de esto
# para TFs de ReMap, que todavia usa el criterio de IC ya que no calcula
# pval_corr).

# Mismo eje x en los dos paneles para que las barras sean comparables.
auc_signif <- auc_df %>% filter(pval_corr < 0.05)
x_lim <- range(c(auc_signif$ci2.5, auc_signif$ci97.5, 0.5)) - 0.5

panel_noise <- function(grupo, y_label) {
  plot_noise_auc_summary(auc_df %>% filter(feature %in% feature_groups$feature[feature_groups$group == grupo])) +
    coord_cartesian(xlim = x_lim) +
    labs(y = y_label)
}
panel_seq <- panel_noise("seq", "Elementos del promotor")
panel_endo <- panel_noise("endo", "Características del promotor")

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}
save_jpg_pdf(panel_seq, "summary_noise_seq", 13.5, 6.75)
save_jpg_pdf(panel_endo, "summary_noise_endo", 13.5, 6.75)

# Figura completa: alto de cada panel proporcional a su numero de features
# (+3 por los ejes, si no el titulo del eje y de endo no entra).
n_feat <- map_int(list(panel_seq, panel_endo), ~ n_distinct(.x$data$feature))
figura <- (panel_seq / panel_endo) +
  plot_layout(heights = n_feat + 3, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top", plot.tag = element_text(size = 28, face = "bold"))
save_jpg_pdf(figura, "summary_noise", 13.5, 2 + 0.45 * sum(n_feat + 3))

# --- Modelo por feature: ruido ~ feature * actividad -----------------------
# Todas las features juntas (sin separar seq/endo). Para cada feature y
# replica: lm(var_rank_sw ~ feature * pct_act), con pct_act = percentil de
# la actividad media dentro de la replica (0-1) centrado en 0.5. Asi:
#   - efecto principal = diferencia en var_rank_sw (posicion 1-100 dentro
#     del bin de actividad) entre promotores con y sin la feature, a
#     actividad mediana;
#   - interaccion = cuanto cambia ese efecto entre el percentil 0 y el 100
#     de actividad.
# Complementario, no va a la tesis: se guarda en
# figures/exploratorio_ruido_feature_x_actividad/ (solo jpg).
# BH sobre todos los tests (features x replicas x terminos). Se muestran las
# features con algun termino significativo y del mismo signo en ambas
# replicas; los terminos que no lo cumplen van en gris.

expl_dir <- "figures/exploratorio_ruido_feature_x_actividad"
dir.create(expl_dir, showWarnings = FALSE, recursive = TRUE)

tidy_mod <- data %>%
  group_by(rep) %>%
  mutate(pct_act = percent_rank(mean) - 0.5) %>%
  ungroup() %>%
  build_tidy_features(keep = c("var_rank_sw", "pct_act"))

terminos <- c(feat = "Efecto principal
(a actividad mediana)", int = "Interacción con la actividad
(cambio del efecto del percentil 0 al 100)")
lm_df <- expand_grid(rep = repname, feature = vbles) %>%
  pmap(function(rep, feature) {
    df <- tidy_mod[tidy_mod$rep == rep, ]
    df$f <- as.numeric(as.logical(df[[feature]]))
    fit <- lm(var_rank_sw ~ f * pct_act, data = df)
    co <- summary(fit)$coefficients[c("f", "f:pct_act"), , drop = FALSE]
    ci <- confint(fit)[c("f", "f:pct_act"), , drop = FALSE]
    tibble(
      feature = feature, rep = rep, termino = c("feat", "int"), n_feature = sum(df$f, na.rm = TRUE),
      estimate = co[, 1], ci_low = ci[, 1], ci_high = ci[, 2], pval = co[, 4]
    )
  }) %>%
  list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, method = "BH"))
write_tsv(lm_df, file.path(expl_dir, "ruido_feature_x_actividad_lm.tsv"))

lm_plot <- lm_df %>%
  group_by(feature, termino) %>%
  mutate(consistente = all(pval_corr < 0.05) & n_distinct(sign(estimate)) == 1) %>%
  group_by(feature) %>%
  filter(any(consistente)) %>%
  ungroup() %>%
  mutate(
    efecto = case_when(!consistente ~ "n.s.", estimate > 0 ~ "Más ruido", TRUE ~ "Menos ruido"),
    termino = factor(terminos[termino], levels = terminos)
  )
orden <- lm_plot %>% filter(termino == terminos["feat"]) %>% group_by(feature) %>% summarise(m = mean(estimate)) %>% arrange(m) %>% pull(feature)

panel_lm <- lm_plot %>%
  mutate(feature = factor(feature, levels = orden), rep = factor(rep, levels = c("Rep 2", "Rep 1"))) %>%
  ggplot(aes(x = estimate, y = feature, group = rep)) +
  geom_col(orientation = "y", position = "dodge", aes(fill = efecto, alpha = rep)) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), position = position_dodge(0.9), height = 0.05, col = "#777777", linewidth = 1) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  facet_wrap(~termino, scales = "free_x") +
  scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7), breaks = c("Rep 1", "Rep 2")) +
  scale_fill_manual(values = c("Más ruido" = "#D6741F", "Menos ruido" = "#7FB800", "n.s." = "grey75"), breaks = c("Más ruido", "Menos ruido", "n.s.")) +
  theme_pubr(base_size = 20) +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(size = 20), panel.spacing = unit(2, "lines")) +
  labs(fill = "Efecto", x = "Coeficiente (rango de varianza dentro del bin, 1-100)", y = "Características del promotor", alpha = "")
ggsave(file.path(expl_dir, "ruido_feature_x_actividad.jpg"), panel_lm, width = 18, height = 2.5 + 0.55 * length(orden), units = "in")

# Efecto predicho de cada feature a lo largo de la actividad
# (efecto principal + interaccion * (percentil - 0.5)), con su banda de
# confianza de 95%, para las mismas features.
efecto_pred <- expand_grid(rep = repname, feature = orden) %>%
  pmap(function(rep, feature) {
    df <- tidy_mod[tidy_mod$rep == rep, ]
    df$f <- as.numeric(as.logical(df[[feature]]))
    fit <- lm(var_rank_sw ~ f * pct_act, data = df)
    grid <- seq(-0.5, 0.5, by = 0.05)
    L <- cbind(0, 1, 0, grid)
    est <- drop(L %*% coef(fit))
    se <- sqrt(rowSums((L %*% vcov(fit)) * L))
    tibble(feature = feature, rep = rep, pct = grid + 0.5, est = est, lo = est - 1.96 * se, hi = est + 1.96 * se)
  }) %>%
  list_rbind()

panel_pred <- efecto_pred %>%
  mutate(feature = factor(feature, levels = rev(orden))) %>%
  ggplot(aes(pct * 100, est, col = rep, fill = rep)) +
  geom_hline(yintercept = 0, linewidth = 0.4, col = "grey50") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, col = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~feature, labeller = label_wrap_gen(25)) +
  scale_color_manual(values = c("Rep 1" = "black", "Rep 2" = "grey55")) +
  scale_fill_manual(values = c("Rep 1" = "black", "Rep 2" = "grey55")) +
  theme_pubr(base_size = 16) +
  theme(legend.position = "top", strip.background = element_blank()) +
  labs(x = "Percentil de actividad media", y = "Efecto de la característica sobre el
rango de varianza (con - sin)", col = "", fill = "")
ggsave(file.path(expl_dir, "ruido_feature_x_actividad_curvas.jpg"), panel_pred, width = 16, height = 3 + 2.6 * ceiling(length(orden) / 4), units = "in")

message("Figuras guardadas en ", out_dir)
