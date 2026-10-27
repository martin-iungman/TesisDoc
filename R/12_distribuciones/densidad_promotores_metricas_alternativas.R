# Exploratorio: variantes de densidad_individual_promotores (extra de
# distribuciones_expresion, R/12_distribuciones/densidad_promotores_individuales.R) que usan
# todos los eventos de citometria.
#
# Con scale_x_log10() se descartan los eventos con senal compensada <= 0
# (10% en BTG1_1, ~45% en LSM1_1/PPP1R14B_3, 48% en el control US), y lo
# que queda positivo tiene casi la misma forma en todas las muestras, asi
# que los promotores parecen superponerse con el control. Aca:
#   1. densidad_asinh: misma figura pero en escala asinh(x / cofactor),
#      lineal cerca de 0 y logaritmica lejos (transformacion estandar en
#      citometria), sin descartar eventos.
#   2. metricas_por_pocillo: por cada pocillo (3 por promotor, 1 para cada
#      control) la mediana de senal menos la del control US, y el % de
#      celulas por encima del percentil 99 del control (% EGFP+).
#
# Salida en figures/exploratorio_densidad_promotores/ (no es una figura
# numerada). Requiere los mismos datos de citometria que
# densidad_promotores_individuales.R. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

out_dir <- "figures/exploratorio_densidad_promotores"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- load_citometry_stable_validation(path_citometry_stable_validation)

cmv <- "CMV (control positivo)"
prom_colors <- c(
  KIAA0753_1 = "#AD343E", TMEM87A_1 = "#FFB400", BTG1_1 = "#7FB800", METAP2_1 = "#D6741F",
  PPP1R14B_3 = "#1B8C8E", ZKSCAN2_1 = "#9EDBCB", LSM1_1 = "#0D2C54", ETS1_1 = "#25afe9",
  setNames("black", cmv)
)
prom_order <- c("KIAA0753_1", "TMEM87A_1", sort(setdiff(names(prom_colors), c("KIAA0753_1", "TMEM87A_1", cmv))), cmv)

control <- df %>% filter(name == "Control") %>% pull(Comp_FL2_A)
proms <- df %>%
  filter(name != "Control") %>%
  mutate(name = factor(recode(name, Strong = cmv), levels = prom_order))

# --- 1. Densidad en escala asinh -------------------------------------------

# Cofactor ~ ancho del ruido del control (su IQR es ~ -20 a 24): por debajo
# la escala es ~lineal, por encima ~log.
cofactor <- 20
asinh_trans <- scales::trans_new(
  "asinh", function(x) asinh(x / cofactor), function(y) sinh(y) * cofactor
)

p_asinh <- proms %>%
  ggplot() +
  geom_density(data = tibble(Comp_FL2_A = control), aes(Comp_FL2_A), col = "grey") +
  geom_density(aes(Comp_FL2_A, col = name)) +
  scale_x_continuous(trans = asinh_trans, breaks = c(-100, 0, 100, 1000, 10000), labels = scales::label_number(big.mark = "")) +
  scale_color_manual(values = prom_colors) +
  facet_wrap(~name) +
  xlab("Señal de EGFP\n(unidades de fluorescencia relativa, escala asinh)") +
  ylab("Densidad") +
  theme_pubclean() +
  theme(legend.position = "none", axis.title = element_text(size = 16), plot.margin = margin(5.5, 15, 5.5, 5.5))
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("densidad_asinh.", ext)), p_asinh, width = 9, height = 6.75, units = "in")
}

# --- 2. Metricas resumen por pocillo ---------------------------------------

umbral_pos <- quantile(control, 0.99)
mediana_control <- median(control)

por_pocillo <- proms %>%
  group_by(name, sample_name) %>%
  summarise(
    n = n(),
    frac_le0 = mean(Comp_FL2_A <= 0),
    mediana = median(Comp_FL2_A),
    pct_pos = 100 * mean(Comp_FL2_A > umbral_pos),
    .groups = "drop"
  ) %>%
  mutate(delta_mediana = mediana - mediana_control)
write_tsv(por_pocillo, file.path(out_dir, "metricas_por_pocillo.tsv"))
message("Umbral EGFP+ (p99 del control) = ", round(umbral_pos, 1), "; mediana del control = ", round(mediana_control, 2))

metric_panel <- function(var, ylab_txt) {
  por_pocillo %>%
    ggplot(aes(name, .data[[var]], col = name)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, col = "grey40") +
    geom_point(size = 2.5) +
    scale_y_log10() +
    scale_color_manual(values = prom_colors) +
    scale_x_discrete(labels = function(x) str_replace(x, " \\(", "\n(")) +
    labs(x = NULL, y = ylab_txt) +
    theme_pubclean() +
    theme(legend.position = "none", axis.title = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1))
}
p_metricas <- cowplot::plot_grid(
  metric_panel("delta_mediana", "Mediana de EGFP − mediana del control\n(escala log)"),
  metric_panel("pct_pos", "% de células EGFP+\n(> percentil 99 del control, escala log)"),
  nrow = 1, labels = c("A", "B")
)
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("metricas_por_pocillo.", ext)), p_metricas, width = 12, height = 6, units = "in", bg = "white")
}

# --- 3. Correlacion library vs citometria (distribuciones_expresion C) -----

# Citometria por pocillo: media en log10 como en distribuciones_expresion
# (descarta los eventos <= 0) y media/SD en asinh (usa todos).
cito_pocillo <- proms %>%
  filter(name != cmv) %>%
  mutate(name = as.character(name), x_asinh = asinh(Comp_FL2_A / cofactor)) %>%
  group_by(name, sample_name) %>%
  summarise(
    media_log = mean(log10(Comp_FL2_A[Comp_FL2_A > 0])),
    media_asinh = mean(x_asinh),
    sd_asinh = sd(x_asinh),
    .groups = "drop"
  )
# Ruido de citometria desacoplado de la media: residuo de SD ~ media entre
# los 8 promotores (analogo, muy rudimentario, a var_rank_sw).
cito_pocillo$sd_asinh_resid <- residuals(lm(sd_asinh ~ media_asinh, data = cito_pocillo))
cito <- cito_pocillo %>%
  group_by(name) %>%
  summarise(across(c(media_log, media_asinh, sd_asinh, sd_asinh_resid), list(m = mean, sd = sd)), .groups = "drop")

# Library: media y varianza entre gates (todas las secuencias, igual que
# distribuciones_expresion) y var_rank_sw (set highconf, como metrica_ruido),
# promediadas entre replicas.
lib_rep <- read_tsv("data/processed/activity_stats_full.tsv", show_col_types = FALSE) %>%
  filter(name %in% cito$name) %>%
  select(name, rep, mean, var)
rank_rep <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE) %>%
  add_noise_rank() %>%
  filter(name %in% cito$name) %>%
  select(name, rep, var_rank_sw)
lib <- lib_rep %>%
  left_join(rank_rep, by = c("name", "rep")) %>%
  group_by(name) %>%
  summarise(
    across(c(mean, var), list(m = mean, sd = ~ sd(.x) %>% replace_na(0))),
    var_rank_m = mean(var_rank_sw, na.rm = TRUE),
    var_rank_sd = sd(var_rank_sw, na.rm = TRUE) %>% replace_na(0),
    n_rep_rank = sum(!is.na(var_rank_sw)),
    .groups = "drop"
  ) %>%
  mutate(var_rank_m = ifelse(n_rep_rank == 0, NA, var_rank_m))
comp <- inner_join(lib, cito, by = "name")
write_tsv(comp, file.path(out_dir, "library_vs_citometria.tsv"))

cor_lab <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  sprintf(
    "n = %d; Pearson r = %.2f; Spearman rho = %.2f",
    sum(ok), cor(x[ok], y[ok]), cor(x[ok], y[ok], method = "spearman")
  )
}
scatter_panel <- function(xv, yv, xlab_txt, ylab_txt, titulo) {
  d <- comp %>% filter(!is.na(.data[[paste0(xv, "_m")]]))
  d %>%
    ggplot(aes(.data[[paste0(xv, "_m")]], .data[[paste0(yv, "_m")]], col = name)) +
    geom_errorbar(aes(ymin = .data[[paste0(yv, "_m")]] - .data[[paste0(yv, "_sd")]], ymax = .data[[paste0(yv, "_m")]] + .data[[paste0(yv, "_sd")]]), width = 0) +
    geom_errorbarh(aes(xmin = .data[[paste0(xv, "_m")]] - .data[[paste0(xv, "_sd")]], xmax = .data[[paste0(xv, "_m")]] + .data[[paste0(xv, "_sd")]]), height = 0) +
    geom_point(size = 3.5) +
    # limites fijos para que los paneles compartan una sola leyenda aunque
    # a alguno le falte un promotor
    scale_color_manual(values = prom_colors, limits = setdiff(prom_order, cmv)) +
    labs(
      x = xlab_txt, y = ylab_txt, col = "Promotor", title = titulo,
      subtitle = cor_lab(d[[paste0(xv, "_m")]], d[[paste0(yv, "_m")]])
    ) +
    theme_bw() +
    theme(text = element_text(size = 14), plot.title = element_text(size = 14))
}

p_media <- (scatter_panel("mean", "media_log", "Media (ensayo high-throughput)", "Media de log10(EGFP)\n(solo eventos > 0)", "Actual (log10)") +
  scatter_panel("mean", "media_asinh", "Media (ensayo high-throughput)", "Media de asinh(EGFP / 20)\n(todos los eventos)", "Asinh")) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("correlacion_media_log_vs_asinh.", ext)), p_media, width = 13, height = 5.5, units = "in")
}

# --- 4. Lo mismo para ruido -------------------------------------------------

p_ruido <- (scatter_panel("var", "sd_asinh", "Varianza entre gates (ensayo high-throughput)", "SD de asinh(EGFP / 20)", "Ruido crudo") +
  scatter_panel("var_rank", "sd_asinh_resid", "Rango de varianza dentro del bin de actividad\n(var_rank_sw, set highconf)", "Residuo de SD ~ media\n(citometría, asinh)", "Ruido desacoplado de la actividad") + theme(legend.position = "none")) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A")
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("correlacion_ruido_asinh.", ext)), p_ruido, width = 13, height = 5.5, units = "in")
}

message("Exploratorio guardado en ", out_dir)
