# puffin_selectivity_predscore (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Validacion contra el modelo de secuencia PUFFIN (Dudnyk et al. 2024):
# actividad del reportero (bins por rango) vs. dos scores de PUFFIN -
#   - selectivity_dudnik: selectividad tisular predicha por Puffin-D
#     sobre el gen endogeno correspondiente (dato externo publicado,
#     tabla suplementaria del paper).
#   - teoTSS_score: score de PUFFIN en la posicion del TSS (pos=16 de
#     252pb) al correr el modelo sobre nuestras propias secuencias de
#     la library (EXCEPCION, ver mas abajo).
# scatter (bins de 100) + boxplot (20 bins de rango) para cada score, y
# selectividad por patron de expresion (housekeeping vs. tejido-especifico).
# Figura completa puffin_selectivity_predscore.jpg/.pdf (2x2, ver abajo).
#
# Requiere: data/processed/activity_stats_highconf.tsv,
# data/processed/prom_df.tsv, data/external/dudnyk2024_sup.xlsx,
# data/external/EPD/human38_epdnew.bed. Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "#PUFFIN") y seq_features.qmd (join selectivity_dudnik contra
# EPD). EXCEPCION: teoTSS_score se lee de Analysis/Tables/
# Dudnyk_puffin_prediction_summ.tsv (path_puffin_pred en
# R/00_prom_features/analysis_tables_exceptions.R) - el pipeline que
# genera esa tabla desde cero (correr el modelo puffin.py sobre nuestra
# library) esta portado pero NO ejecutado en
# R/25_puffin/puffin_prediction_model.R, ver ese archivo para el motivo.

library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(patchwork)

# GenomicRanges/S4Vectors define S4 generics that shadow dplyr's
# select/filter/rename; pin them to the dplyr versions (same guard as
# R/00_prom_features/build_prom_features.R).
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

source("R/00_prom_features/analysis_tables_exceptions.R")
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "puffin_selectivity_predscore"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

# --- selectivity_dudnik: join externo contra la tabla suplementaria de
# Dudnyk et al. 2024 (Puffin-D), por posicion de TSS contra EPD ------------

epd <- import.bed("data/external/EPD/human38_epdnew.bed")
selectivity_puffin <- readxl::read_xlsx("data/external/dudnyk2024_sup.xlsx", sheet = 5, skip = 1) %>%
  rename(
    start = TSS,
    selectivity_dudnik = `Selectivity (Puffin-D)`,
    motif_selectivity_dudnik = `Motif Selectivity`,
    dispersion_tissue_dudnik = `Dispersion index (log10)`,
    mean_expr_log10_dudnik = `Mean expression (log10)`
  ) %>%
  select(-geneID, -geneName) %>%
  distinct()
selectivity_puffin$end <- selectivity_puffin$start
gr_selectivity_puffin <- GRanges(selectivity_puffin)
selectivity_puffin_lib <- plyranges::join_overlap_intersect_directed(
  promoters(epd, upstream = 10, downstream = 10), gr_selectivity_puffin
) %>%
  values() %>%
  as_tibble() %>%
  add_count(name) %>%
  filter(n == 1) %>%
  select(-c(n, score))
prom_df <- left_join(prom_df, selectivity_puffin_lib, by = "name")

# --- teoTSS_score: prediccion de PUFFIN sobre nuestra library (EXCEPCION) --

df_pred <- read_tsv(path_puffin_pred, show_col_types = FALSE) %>% select(seq_id, teoTSS_score)
prom_df <- left_join(prom_df, df_pred, by = "seq_id")

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))
data_sw <- add_mean_sw_bins(data)

save_jpg_pdf <- function(plot, name, width = 9, height = 6.75) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}

# IC 95% de rho de Spearman por bootstrap (igual que especificidad_tisular_fantom)
set.seed(1)
rho_boot_ci <- function(x, y, n_iter = 1000) {
  r <- replicate(n_iter, {
    i <- sample(length(x), replace = TRUE)
    cor(x[i], y[i], method = "spearman")
  })
  list(lo = unname(quantile(r, 0.025)), hi = unname(quantile(r, 0.975)))
}

# rho de Spearman entre actividad del reportero (mean) y el score, sobre
# promotores individuales (no bins), por replica y grupo
rho_labels <- function(df, score_col, ...) {
  df %>%
    filter(!is.na(.data[[score_col]])) %>%
    group_by(rep, ...) %>%
    summarise(
      rho = cor(mean, .data[[score_col]], method = "spearman"),
      ci = list(rho_boot_ci(mean, .data[[score_col]])),
      n = n(),
      .groups = "drop"
    ) %>%
    unnest_wider(ci) %>%
    mutate(label = sprintf("rho == \"%.2f [%.2f, %.2f]; n = %d\"", rho, lo, hi, n))
}

lab_x_bins <- "Bins de actividad del promotor reportero"
lab_sel <- "Índice de selectividad (Puffin-D)"
lab_pred <- "Puntaje de Puffin en el TSS"
tit_sel <- "Selectividad predicha del gen endógeno"
tit_pred <- "Predicción de Puffin sobre la secuencia del fragmento"

# --- Scatter (bins de 100 por rango): mediana del score por bin ------------

puffin_scatter <- function(data_sw, score_col, ylab, titulo) {
  rho <- rho_labels(data_sw, score_col)
  data_sw %>%
    filter(!is.na(.data[[score_col]])) %>%
    group_by(rep, mean_sw) %>%
    summarise(puffin = median(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(mean_sw, puffin)) +
    geom_point(col = "#14AFB2") +
    geom_smooth(col = "#216869", method = "lm") +
    geom_text(data = rho, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.5, size = 4, parse = TRUE, inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    facet_wrap(~rep) +
    labs(x = paste(lab_x_bins, "(bins de 100)"), y = ylab, title = titulo) +
    theme_bw(base_size = 15)
}

p_sel_scatter <- puffin_scatter(data_sw, "selectivity_dudnik", lab_sel, tit_sel)
p_pred_scatter <- puffin_scatter(data_sw, "teoTSS_score", paste(lab_pred, "(mediana)"), tit_pred)
save_jpg_pdf(p_sel_scatter, "PUFFIN_selectivity_scatter")
save_jpg_pdf(p_pred_scatter, "PUFFIN_predscore_scatter")

# --- Boxplot (20 bins de rango, sin outliers) ------------------------------

puffin_boxplot <- function(data, score_col, ylab, titulo) {
  rho <- rho_labels(data, score_col)
  data %>%
    filter(!is.na(.data[[score_col]])) %>%
    group_by(rep) %>%
    mutate(bin = row_number(mean) %>% cut_number(20) %>% as.numeric() %>% as_factor()) %>%
    ggplot(aes(bin, .data[[score_col]])) +
    geom_boxplot(alpha = 0.6, fill = "#216869", outliers = FALSE) +
    geom_text(data = rho, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.5, size = 4, parse = TRUE, inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    facet_wrap(~rep) +
    scale_x_discrete(breaks = c("1", "5", "10", "15", "20")) +
    labs(x = paste(lab_x_bins, "(por rango)"), y = ylab, title = titulo) +
    theme_bw(base_size = 15) +
    theme(panel.spacing = unit(1.5, "lines"))
}

p_sel_box <- puffin_boxplot(data, "selectivity_dudnik", lab_sel, tit_sel)
p_pred_box <- puffin_boxplot(data, "teoTSS_score", lab_pred, tit_pred)
save_jpg_pdf(p_sel_box, "PUFFIN_selectivity_boxplot")
save_jpg_pdf(p_pred_box, "PUFFIN_predscore_boxplot")

# --- Selectividad por patron de expresion (housekeeping vs. tejido-
# especifico): terciles de sample_specificity_gini por replica sobre todos
# los promotores con Gini, igual que especificidad_tisular_fantom (los cortes
# no dependen de tener selectividad de Puffin-D); actividad reordenada dentro
# de cada grupo y agrupada en 30 bins ----------------------------------------

patron_colors <- c("Housekeeping" = "#1B8C8E", "Tejido-específico" = "#0D2C54")
data_patron <- data %>%
  filter(!is.na(sample_specificity_gini)) %>%
  group_by(rep) %>%
  mutate(specificity_tercile = cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) %>%
  ungroup() %>%
  filter(specificity_tercile != 2, !is.na(selectivity_dudnik)) %>%
  mutate(patron = ifelse(specificity_tercile == 1, "Housekeeping", "Tejido-específico"))

rho_patron <- rho_labels(data_patron, "selectivity_dudnik", patron) %>%
  mutate(vjust = ifelse(patron == "Housekeeping", 1.5, 3))

p_selectivity_pattern <- data_patron %>%
  group_by(rep, patron) %>%
  mutate(activity_bin = row_number(mean) %>% cut_number(30) %>% as.numeric()) %>%
  group_by(rep, patron, activity_bin) %>%
  summarise(selectivity_dudnik = median(selectivity_dudnik), .groups = "drop") %>%
  ggplot(aes(activity_bin, selectivity_dudnik, col = patron, fill = patron)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm") +
  geom_text(
    data = rho_patron, aes(x = -Inf, y = Inf, label = label, colour = patron, vjust = vjust),
    hjust = -0.05, size = 4, parse = TRUE, inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  facet_wrap(~rep) +
  scale_color_manual(values = patron_colors) +
  scale_fill_manual(values = patron_colors) +
  labs(
    x = paste(lab_x_bins, "(30 bins por grupo)"), y = lab_sel,
    title = "Selectividad predicha según patrón de expresión",
    col = "Patrón de expresión", fill = "Patrón de expresión"
  ) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")
save_jpg_pdf(p_selectivity_pattern, "PUFFIN_selectivity_by_pattern")

# --- Mapa de calor: score por decil de actividad x decil de Gini -----------
# Mediana del score por celda (numero = promotores en la celda); deciles por
# replica sobre los promotores con Gini y score. En el titulo de cada faceta,
# los coeficientes de lm(score ~ rango de actividad * Gini), todo
# estandarizado (rango de actividad = percent_rank de mean); tabla en
# puffin_heatmap_lm.tsv.

puffin_heatmap <- function(data, score_col, fill_lab, titulo) {
  d <- data %>%
    filter(!is.na(sample_specificity_gini), !is.na(.data[[score_col]])) %>%
    group_by(rep) %>%
    mutate(
      rango_act = percent_rank(mean),
      dec_act = ntile(mean, 10),
      dec_gini = ntile(sample_specificity_gini, 10)
    ) %>%
    ungroup()

  coefs <- d %>%
    group_by(rep) %>%
    group_modify(~ {
      m <- lm(scale(.x[[score_col]]) ~ scale(rango_act) * scale(sample_specificity_gini), data = .x)
      co <- summary(m)$coefficients[-1, ]
      tibble(termino = c("actividad", "gini", "interaccion"), beta = co[, 1], se = co[, 2], pval = co[, 4], n = nobs(m), r2 = summary(m)$r.squared)
    }) %>%
    ungroup() %>%
    mutate(score = score_col)
  rep_lab <- coefs %>%
    group_by(rep) %>%
    summarise(lab = sprintf(
      "%s\nβ actividad = %.2f; β Gini = %.2f\nβ interacción = %.2f",
      rep[1], beta[termino == "actividad"], beta[termino == "gini"], beta[termino == "interaccion"]
    )) %>%
    deframe()
  # mediana de Gini de cada decil (pooled) para rotular el eje y
  gini_lab <- d %>%
    group_by(dec_gini) %>%
    summarise(lab = sprintf("%d (%.2f)", dec_gini[1], median(sample_specificity_gini))) %>%
    pull(lab)

  p <- d %>%
    group_by(rep, dec_act, dec_gini) %>%
    summarise(valor = median(.data[[score_col]]), n = n(), .groups = "drop") %>%
    ggplot(aes(dec_act, dec_gini, fill = valor)) +
    geom_tile() +
    geom_text(aes(label = n), size = 2.6, colour = "white") +
    facet_wrap(~rep, labeller = as_labeller(rep_lab)) +
    scale_fill_viridis_c(name = fill_lab) +
    scale_x_continuous(breaks = 1:10, expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:10, labels = gini_lab, expand = c(0, 0)) +
    labs(x = "Decil de actividad del promotor reportero", y = "Decil de Gini (mediana)", title = titulo) +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank())
  list(plot = p, coefs = coefs)
}

hm_sel <- puffin_heatmap(data, "selectivity_dudnik", "Selectividad\n(mediana)", paste(tit_sel, "según especificidad tisular"))
hm_pred <- puffin_heatmap(data, "teoTSS_score", "Puntaje\n(mediana)", "Predicción de Puffin según especificidad tisular")
save_jpg_pdf(hm_sel$plot, "PUFFIN_selectivity_heatmap_gini")
save_jpg_pdf(hm_pred$plot, "PUFFIN_predscore_heatmap_gini")
write_tsv(bind_rows(hm_sel$coefs, hm_pred$coefs), file.path(out_dir, "puffin_heatmap_lm.tsv"))
print(as.data.frame(bind_rows(hm_sel$coefs, hm_pred$coefs)))

# Figura completa: arriba prediccion de Puffin sobre nuestras secuencias,
# boxplot (A) y mapa de calor por especificidad tisular (B); abajo
# selectividad Puffin-D del gen endogeno, boxplot (C) y mapa de calor (D). Los
# scatters y la selectividad por patron de expresion quedan como paneles
# sueltos.
figura <- p_pred_box + hm_pred$plot + p_sel_box + hm_sel$plot +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
save_jpg_pdf(figura, slug, 18, 14)

message("Figuras guardadas en ", out_dir)
