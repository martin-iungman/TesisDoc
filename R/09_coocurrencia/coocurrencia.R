# coocurrencia_motivos (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Co-ocurrencia entre features de promotores: frecuencia de cada feature
# sobre su propio universo no-NA (el universo difiere entre features,
# plot_feature_frequency()), heatmap de coeficiente phi (correlacion de
# Pearson sobre variables binarias - contempla co-presencia y co-ausencia,
# apropiado para prevalencias moderadas) y heatmap de log2(lift)
# (observado/esperado bajo independencia, P(A∩B)/(P(A)*P(B)) - a
# diferencia de phi no tiene techo cuando las prevalencias estan lejos de
# 50/50, revela relaciones fuertes entre features raras que phi subestima;
# pseudocount de Haldane-Anscombe para evitar lift=0 en pares mutuamente
# excluyentes), cada uno con todas las features y solo con las
# estadisticamente significativas (misma direccion de efecto en ambas
# replicas y ambos limites del IC, ver R7/activity_summary.tsv) - 6
# figuras en total. Jaccard y overlap (Szymkiewicz-Simpson) se evaluaron
# tambien (ver plot_cooccurrence_jaccard()/plot_cooccurrence_overlap() en
# R/functions/plot_helpers.R) pero no se incluyen en esta figura final.
#
# build_cooccurrence_features()/plot_cooccurrence_phi()/
# plot_cooccurrence_lift()/plot_feature_frequency() en
# R/functions/plot_helpers.R, compartidos con R5.7 (promalt_coocurrencia).
# Incluye "Promotor unidireccional" (orientacion de PRO-cap, ver orientacion_promotor
# / R/00_prom_features/build_promoter_orientation.R) como feature mas de
# build_cooccurrence_features() - resync 2026-09-30: dejo de agregarse
# localmente para estar disponible en todos los analisis de features
# (R7/activity_summary, R3.2/ruido_summary, R5.7/promalt), no solo aca.
#
# Requiere: data/processed/prom_df.tsv (ver R/00_prom_features) y
# data/processed/activity_summary.tsv (ver
# R/15_summary_features/summary_features.R). Run from the TesisDoc repo root.

library(tidyverse)
library(fastDummies)
library(patchwork)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")

# Dos figuras: coocurrencia_motivos (frecuencia de cada feature, un panel) y
# coocurrencia_motivos_heatmaps (phi y log2 lift).
slug <- "coocurrencia_motivos"
slug_heat <- "coocurrencia_motivos_heatmaps"
out_dir <- fig_dir(slug)
out_dir_heat <- fig_dir(slug_heat)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
act_summary <- read_tsv("data/processed/activity_summary.tsv", show_col_types = FALSE)
act_signif_features <- act_summary %>%
  mutate(sign = sign(estimate)) %>%
  filter(pval_corr < 0.05) %>%
  count(feature, sign) %>%
  filter(n == 6) %>%
  pull(feature)

# Correcciones de rotulos y de una definicion SOLO para estas dos figuras
# (build_cooccurrence_features() es compartida con R7, R3.2, etc., que
# siguen con los nombres viejos - decision 2026-11-20):
# - "Sin actividad en ratón" de la funcion compartida convierte NA de
#   turnover en FALSE (universo = toda la library). Aca se recalcula con
#   universo = los promotores clasificados por Young (turnover no-NA) y
#   se rotula segun lo que incluye: expression-turnover + mouse-diminished.
# - TSS fuerte = INR_strong_TSS (patron YCASW) -> "Inr fuerte".
# - Ventanas de PhyloP con la geometria actual (TSS = 0 en la posicion 237,
#   fragmento -236 a +15): los cortes de build_prom_features.R caen en
#   -51/+15, -151/-52 y -236/-152.
# - "Insertado en humanos" -> "Inserción en humanos", como en
#   conservacion_evolucion_library.
RATON <- "Actividad perdida o reducida en ratón"
renombres <- c(
  "Inr fuerte" = "TSS fuerte",
  "Inserción en humanos" = "Insertado en humanos",
  "Alta conservación (−51 a +15 pb)" = "Alta conservación (16 a -50pb)",
  "Alta conservación (−151 a −52 pb)" = "Alta conservación (-50 a -150)",
  "Alta conservación (−236 a −152 pb)" = "Alta conservación (-150 a -235)",
  setNames("Sin actividad en ratón", RATON)
)

raton <- prom_df %>%
  transmute(seq_id, raton = ifelse(is.na(turnover), NA, turnover %in% c("expression-turnover", "mouse-diminished")))

tidy_data <- build_cooccurrence_features(prom_df) %>%
  select(-`Sin actividad en ratón`) %>%
  left_join(raton, by = "seq_id") %>%
  rename(!!RATON := raton) %>%
  rename(any_of(renombres[names(renombres) != RATON]))
stopifnot(!any(renombres %in% names(tidy_data)))

# activity_summary.tsv usa los nombres viejos
act_signif_features <- recode(act_signif_features, !!!setNames(names(renombres), renombres))

save_fig <- function(p, name, width, height, dir = out_dir) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(dir, paste0(name, ".", ext)), p, width = width, height = height, units = "in")
  }
}

save_pair <- function(data_for_matrix, suffix) {
  feat_mat <- data_for_matrix %>%
    select(-seq_id) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  p_phi <- plot_cooccurrence_phi(feat_mat)
  p_lift <- plot_cooccurrence_lift(feat_mat)
  p_freq <- plot_feature_frequency(feat_mat)
  freq_height <- pmax(4, ncol(feat_mat) * 0.28 + 1.5)
  save_fig(p_phi, paste0("cooccurrence_phi", suffix), 9, 6.75, out_dir_heat)
  save_fig(p_lift, paste0("cooccurrence_lift", suffix), 9, 6.75, out_dir_heat)
  save_fig(p_freq, paste0("feature_frequency", suffix), 8, freq_height)
  list(freq = p_freq, phi = p_phi, lift = p_lift, mat = feat_mat, freq_height = freq_height)
}

# --- 6 figuras: todas las features / solo las significativas -------------

todas <- save_pair(tidy_data, "")
invisible(save_pair(tidy_data %>% select(seq_id, any_of(act_signif_features)), "_signif"))

# --- Figuras finales (todas las features), sin titulo de panel ------------

sin_titulo <- labs(title = NULL)

# coocurrencia_motivos: un solo panel (sin letra)
save_fig(todas$freq + sin_titulo, slug, 8, todas$freq_height)

# coocurrencia_motivos_heatmaps: A phi arriba, B lift abajo, en vertical para
# que el texto de los ejes quede legible a ancho de pagina
texto_heat <- theme(axis.text = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(size = 13))
heatmaps <- (todas$phi + sin_titulo + texto_heat) / (todas$lift + sin_titulo + texto_heat) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))
save_fig(heatmaps, slug_heat, 11, 19, out_dir_heat)

# Estadisticos citados en el mapping (phi de Islas CpG con otras features)
phi_mat <- suppressWarnings(cor(todas$mat, use = "pairwise.complete.obs"))
diag(phi_mat) <- NA
print(round(sort(phi_mat["Islas CpG", ], decreasing = TRUE)[1:4], 2))
message("Maximo phi global: ", round(max(phi_mat, na.rm = TRUE), 2))

message("Figuras guardadas en ", out_dir, " y ", out_dir_heat)
