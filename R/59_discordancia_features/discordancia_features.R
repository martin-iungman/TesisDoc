# discordancia_features (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Efecto de cada caracteristica booleana del promotor sobre la discordancia
# con signo entre el reportero y la actividad endogena maxima (FANTOM5), en
# promotores activos en HEK293 (hek_tpm > 0):
#
#   dif_signed ~ avg_rank + rep + caracteristica
#
# con errores estandar robustos agrupados por promotor, porque las dos replicas
# de un promotor no son independientes (ajustar_discordancia_robusto() en
# R/functions/plot_helpers.R, que explica por que no se usa (1 | seq_id)).
# Correccion BH sobre las caracteristicas de la figura.
#
# CONVENCION DE SIGNO: dif_signed = rr - re (rangos de reportero y de
# max_tpm, normalizados por replica; add_rank_discordance() en
# R/functions/plot_helpers.R). Positivo = el reportero SOBREESTIMA respecto
# a la actividad endogena; negativo = la SUBESTIMA. El coeficiente es el
# contraste con caracteristica - sin caracteristica.
#
# DISENO (analisis tipo Bland-Altman): avg_rank = (rr + re) / 2 es el nivel
# de actividad, casi ortogonal a dif_signed. Controlar solo por el rango del
# reportero dejaba pasar al coeficiente el efecto de la caracteristica sobre
# la actividad cruda. La magnitud |dif_signed| sigue confundida con la
# actividad aun con este control, por eso se reporta la discordancia con signo.
#
# ORDEN DE LOS PASOS: los rangos (rr, re) se calculan dentro de los activos en
# HEK293. Las caracteristicas se calculan ANTES de filtrar, sobre todos los
# promotores de alta confianza, asi que los terciles (G+C, conservacion,
# especificidad, forma, DNase) usan los mismos cortes que summary_features_
# secuencia, summary_features_endogeno y especificidad_tisular_fantom.
#
# CARACTERISTICAS FUERA DE LA FIGURA (siguen en la tabla, con su motivo):
# - Menos de min_promotores promotores activos en HEK293 con la caracteristica:
#   grupo demasiado chico para reportar ("No detectado (FANTOM5)", "Sin modulo
#   cis-regulatorio", "LINE" y "Tracto poli(dA:dT) proximal").
# - "Alta actividad en HEK293": no es una tautologia algebraica con max_tpm,
#   pero funciona como un indicador de re. A avg_rank fijo, un re mas alto
#   implica un dif_signed mas bajo, asi que su coeficiente negativo sale de
#   la construccion. El nulo de la ultima seccion lo mide: permutando hek_tpm
#   entre promotores de max_tpm similar (conserva solo el acople con re) el
#   coeficiente es tan o mas negativo que el observado.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv y data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION, ver
# R/00_prom_features/analysis_tables_exceptions.R). Para el diagnostico final,
# data/processed/activity_summary.tsv (R/15_summary_features). Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/rank_dif_lmm.R.

library(tidyverse)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "discordancia_features"
out_dir <- fig_dir(slug)

min_promotores <- 100
feature_mecanica <- "Alta actividad en HEK293"

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data_highconf <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id")

data <- data_highconf %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  add_rank_discordance()

message(
  "n = ", n_distinct(data$seq_id), " promotores activos en HEK293 de ",
  n_distinct(data_highconf$seq_id), " de alta confianza (", nrow(data), " mediciones)"
)

# --- Modelo por caracteristica ------------------------------------------------

features_highconf <- build_tidy_features(data_highconf, keep = "seq_id")
caracteristicas <- setdiff(names(features_highconf), c("seq_id", "rep"))

tidy_data <- data %>%
  transmute(seq_id, rep = factor(rep, levels = levels(features_highconf$rep)), dif_signed, avg_rank) %>%
  inner_join(features_highconf, by = c("seq_id", "rep"), relationship = "one-to-one")
stopifnot(nrow(tidy_data) == nrow(data))

ajustar_caracteristica <- function(var) {
  datos <- tidy_data %>%
    transmute(seq_id, rep, dif_signed, avg_rank, presente = as.character(.data[[var]])) %>%
    filter(!is.na(presente)) %>%
    mutate(presente = factor(presente, levels = c("FALSE", "TRUE")))
  ajustar_discordancia_robusto(datos) %>%
    transmute(
      caracteristica = var, estimate, conf.low, conf.high, std.error, statistic, df, p_value,
      n_promotores = n_distinct(datos$seq_id[datos$presente == "TRUE"]),
      n_mediciones = nrow(datos)
    )
}

resultados <- map(caracteristicas, ajustar_caracteristica) %>%
  list_rbind() %>%
  mutate(
    motivo_exclusion = case_when(
      caracteristica == feature_mecanica ~ "acoplada a re por construccion",
      n_promotores < min_promotores ~ paste0("menos de ", min_promotores, " promotores"),
      TRUE ~ NA_character_
    ),
    en_figura = is.na(motivo_exclusion)
  )

resultados <- resultados %>%
  left_join(
    resultados %>% filter(en_figura) %>% transmute(caracteristica, p_adj = p.adjust(p_value, method = "BH")),
    by = "caracteristica"
  ) %>%
  arrange(desc(en_figura), desc(estimate))

write_tsv(resultados, "tables/discordancia_features_lmm.tsv")

# --- Figura ---------------------------------------------------------------------

colores_fdr <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")

p <- resultados %>%
  filter(en_figura) %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "FDR < 0.001",
      p_adj < 0.01 ~ "FDR < 0.01",
      p_adj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "ns"
    ) %>% factor(levels = names(colores_fdr)),
    caracteristica = fct_reorder(caracteristica, estimate, .desc = TRUE)
  ) %>%
  ggplot(aes(x = estimate, y = caracteristica, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.3, alpha = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = colores_fdr) +
  labs(x = "Efecto (β) e IC 95% (positivo: el reportero sobreestima)", y = NULL, color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

for (ext in c("jpg", "pdf")) {
  ggsave(
    file.path(out_dir, paste0(slug, ".", ext)), p,
    width = 9, height = 6.75, units = "in",
    device = if (ext == "pdf") cairo_pdf else "jpeg"
  )
}

# --- Diagnostico: correlacion con el efecto sobre la actividad cruda ------------

activity_summary <- read_tsv("data/processed/activity_summary.tsv", show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(efecto_actividad = mean(estimate), .groups = "drop")

comparacion <- resultados %>%
  filter(en_figura) %>%
  inner_join(activity_summary, by = c("caracteristica" = "feature"))
message(
  "Correlacion con el efecto sobre la actividad cruda: r = ",
  round(cor(comparacion$estimate, comparacion$efecto_actividad), 2),
  " (", nrow(comparacion), " caracteristicas)"
)

# --- Nulo mecanico de la caracteristica excluida --------------------------------
# Permuta hek_tpm entre promotores del mismo percentil de max_tpm: conserva el
# acople de la caracteristica con re y rompe cualquier otra relacion con el
# reportero. Solo importa el coeficiente, asi que alcanza con lm.

coef_alta_hek <- function(df) {
  unname(coef(lm(dif_signed ~ avg_rank + rep + alta, data = df))["altaTRUE"])
}
mediana_hek <- median(data$hek_tpm)
por_promotor <- data %>%
  distinct(seq_id, hek_tpm, max_tpm) %>%
  mutate(percentil_max = ntile(max_tpm, 100))

set.seed(1)
nulo <- map_dbl(1:200, function(i) {
  permutado <- por_promotor %>%
    group_by(percentil_max) %>%
    mutate(hek_perm = sample(hek_tpm)) %>%
    ungroup() %>%
    select(seq_id, hek_perm)
  data %>%
    left_join(permutado, by = "seq_id") %>%
    mutate(alta = hek_perm > mediana_hek) %>%
    coef_alta_hek()
})

observado <- data %>% mutate(alta = hek_tpm > mediana_hek) %>% coef_alta_hek()
write_tsv(
  tibble(caracteristica = feature_mecanica, observado, nulo_media = mean(nulo), nulo_sd = sd(nulo), n_permutaciones = length(nulo)),
  "tables/discordancia_features_nulo_alta_hek.tsv"
)
message(feature_mecanica, ": observado = ", round(observado, 3), ", nulo = ", round(mean(nulo), 3), " (DE ", round(sd(nulo), 3), ")")

message("Figuras guardadas en ", out_dir)
