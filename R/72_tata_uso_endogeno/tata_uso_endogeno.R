# tata_uso_endogeno (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
#
# Los promotores con TATA fuerte que NO son activos en el reportero: que les
# pasa en su contexto endogeno. Version de figura de p7 (parte inferior), p12
# y p13 del exploratorio R/44_exploratorio_tata_ruido (ver su handoff.md); el
# script recalcula todo desde los datos.
#
# Paneles:
# (A) Percentil de actividad en el reportero (media, IC 95% de t) por quintil
#     del indice de Gini de FANTOM5 (sample_specificity_gini; quintiles sobre
#     todos los promotores con Gini de cada replica, x = mediana del quintil),
#     TATA fuertes vs. resto, n de las TATA fuertes en cada punto. Las TATA
#     fuertes son mas activas que el resto en todos los quintiles: el motivo
#     sigue aportando actividad en los promotores tejido-especificos.
# (B) Solo TATA fuertes, actividad baja vs. alta en el reportero: TPM maximo
#     entre las muestras de tejido + celula primaria de FANTOM5 (max_tpm) y TPM
#     en CAGE de HEK293 (hek_tpm), en escala log (TPM + 1), con n y mediana.
#     Las de baja actividad alcanzan niveles altos en algun tejido pero estan
#     mayormente apagadas en HEK293.
#
# Estadisticos para el texto (tsv, mismos modelos que p12 y p13 de R/44),
# sobre todos los promotores, TATA fuertes vs. resto:
#   - lm(log2(TPM + 1) ~ act_alta * tata_fuerte) para max_tpm y hek_tpm:
#     efecto de la TATA fuerte en baja y en alta actividad e interaccion (IC 95%).
#   - glm(detectado ~ act_alta * tata_fuerte): FANTOM5 detectado = max_tpm > 1
#     (umbral de active_fantom5_samples), HEK293 detectado = hek_tpm > 0.
#   - Resumen por celda TATA fuerte x actividad (n, mediana, cuartiles, sin
#     senal, > 100 TPM).
#
# Definiciones (identicas a R/44_exploratorio_tata_ruido y tata_ruido_especificidad):
#   - Mismos promotores: add_noise_rank() (descarta el resto de n/100 de cada
#     replica) y fragmentos sin N, para que los n coincidan con R3.
#   - TATA fuerte: canonica (Bucher + pesos W de EPD, corte -8.17, -31 a -26
#     misma hebra) con puntaje >= la mediana de las canonicas.
#   - Actividad alta: percentil 75 de mean dentro de cada replica.
#
# Salidas: fig_dir("tata_uso_endogeno") (figura completa y paneles, jpg y pdf;
# un tsv por panel y uno por modelo).
# Requiere data/raw/library.bed, BSgenome.Hsapiens.UCSC.hg38,
# data/processed/activity_stats_highconf.tsv, data/processed/prom_df.tsv y
# data/processed/fantom_endo_activity_summary.tsv (EXCEPCION, ver
# R/00_prom_features/build_fantom_endo_activity_summary.R).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/functions/tata_scan.R")

select <- dplyr::select
filter <- dplyr::filter

out_dir <- fig_dir("tata_uso_endogeno")

CANON <- TATA_CANON
col_tata <- "#216869"
col_hek <- "#14AFB2"
clr_grupo <- c("TATA fuertes" = col_tata, "Resto" = "grey55")
clr_medida <- c("Máximo en FANTOM5" = col_tata, "HEK293" = col_hek)
th <- theme_pubr(base_size = 16) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")
tag_theme <- theme(plot.tag = element_text(size = 24, face = "bold"))
num <- function(x) format(x, decimal.mark = ".", big.mark = "", scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}
save_tsv <- function(x, name) write_tsv(x, file.path(out_dir, name))

# --- Datos ----------------------------------------------------------------------

lib <- tata_load_library()
code <- lib$code[!lib$has_n, ] # 1 fragmento con N fuera
sc <- tata_scores(code)
hits <- tata_hits(sc$misma, sc$opuesta, TATA_CORTE_EPD)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_endo <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  filter(seq_id %in% rownames(code)) %>%
  group_by(rep) %>%
  mutate(act_pct = 100 * percent_rank(mean), act_alta = act_pct >= 75) %>%
  ungroup() %>%
  inner_join(tata_classify(hits, rownames(code)), by = "seq_id") %>%
  mutate(score_canon = tata_best_in_window(sc$misma, CANON[1], CANON[2])[seq_id])
stopifnot(all((data$clase == "Canónica") == data$TATA_EPD))

corte_puntaje <- round(median(data$score_canon[data$clase == "Canónica"]), 2)
message("Corte de puntaje (mediana de las canonicas): ", corte_puntaje)

dat <- data %>%
  left_join(fantom_endo %>% select(seq_id, max_tpm), by = "seq_id") %>%
  mutate(
    tata_fuerte = clase == "Canónica" & score_canon >= corte_puntaje - TATA_TOL,
    grupo = factor(if_else(tata_fuerte, "TATA fuertes", "Resto"), levels = c("TATA fuertes", "Resto")),
    actividad = factor(if_else(act_alta, "Alta", "Baja"), levels = c("Baja", "Alta"))
  ) %>%
  select(seq_id, rep, act_pct, act_alta, actividad, tata_fuerte, grupo, gini = sample_specificity_gini, max_tpm, hek_tpm)
stopifnot(!anyNA(dat$max_tpm), !anyNA(dat$hek_tpm))
print(count(dat, rep, grupo, actividad))

# --- Panel A: actividad por quintil de Gini -------------------------------------

panel_a_df <- dat %>%
  filter(!is.na(gini)) %>%
  group_by(rep) %>%
  mutate(quintil = ntile(gini, 5)) %>%
  group_by(rep, quintil) %>%
  mutate(gini_mediana = median(gini)) %>%
  group_by(rep, grupo, quintil, gini_mediana) %>%
  summarise(n = n(), act_media = mean(act_pct), h = qt(0.975, n() - 1) * sd(act_pct) / sqrt(n()), prop_alta = mean(act_alta), .groups = "drop") %>%
  mutate(conf.low = act_media - h, conf.high = act_media + h) %>%
  select(-h)
save_tsv(panel_a_df, "panel_a_actividad_quintil_gini.tsv")
print(panel_a_df, n = Inf, width = Inf)

dodge_a <- position_dodge(width = 0.02)
panel_a <- ggplot(panel_a_df, aes(gini_mediana, act_media, color = grupo)) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, linewidth = 0.9, position = dodge_a) +
  geom_line(position = dodge_a) +
  geom_point(size = 3, position = dodge_a) +
  geom_text(data = filter(panel_a_df, grupo == "TATA fuertes"), aes(y = conf.low, label = paste0("n=", n)), vjust = 1.6, size = 4.2, show.legend = FALSE) +
  scale_color_manual(values = clr_grupo) +
  scale_x_continuous(labels = num) +
  facet_wrap(~rep) +
  th +
  labs(x = "Especificidad tisular (índice de Gini en FANTOM5, mediana del quintil)", y = "Actividad en el reportero\n(percentil medio, IC 95%)", color = NULL)
save_jpg_pdf(panel_a, "panel_a_actividad_quintil_gini", 12, 5.5)

# --- Panel B: uso endogeno de las TATA fuertes ----------------------------------

medidas <- c(max_tpm = "Máximo en FANTOM5", hek_tpm = "HEK293")
largo <- dat %>%
  pivot_longer(c(max_tpm, hek_tpm), names_to = "variable", values_to = "tpm") %>%
  mutate(medida = factor(medidas[variable], levels = medidas), umbral = if_else(variable == "max_tpm", 1, 0))

panel_b_celdas <- largo %>%
  group_by(medida, rep, grupo, actividad) %>%
  summarise(
    n = n(), mediana = median(tpm), q1 = quantile(tpm, 0.25), q3 = quantile(tpm, 0.75),
    n_cero = sum(tpm == 0), n_no_detectado = sum(tpm <= umbral), prop_no_detectado = mean(tpm <= umbral),
    mediana_detectados = median(tpm[tpm > umbral]), prop_mayor_100 = mean(tpm > 100), .groups = "drop"
  )
save_tsv(panel_b_celdas, "panel_b_uso_endogeno_por_celda.tsv")
print(panel_b_celdas, n = Inf, width = Inf)

panel_b_df <- filter(largo, tata_fuerte)
etiq_b <- panel_b_celdas %>%
  filter(grupo == "TATA fuertes") %>%
  mutate(lab = paste0("n=", n, "\n", num(round(mediana, 1))))
dodge_b <- position_dodge(width = 0.85)
panel_b <- ggplot(panel_b_df, aes(actividad, tpm + 1, fill = medida)) +
  geom_violin(alpha = 0.5, color = NA, scale = "width", position = dodge_b) +
  geom_boxplot(width = 0.15, outlier.shape = NA, position = dodge_b, show.legend = FALSE) +
  geom_text(data = etiq_b, aes(actividad, 4 * max(panel_b_df$tpm + 1), label = lab, group = medida), position = dodge_b, size = 4.2, lineheight = 0.9) +
  scale_y_log10(breaks = 10^(0:5), labels = num) +
  scale_fill_manual(values = clr_medida) +
  facet_wrap(~rep) +
  th +
  labs(x = "Actividad en el reportero (TATA fuertes)", y = "Actividad endógena\n(TPM de CAGE + 1, escala log)", fill = NULL)
save_jpg_pdf(panel_b, "panel_b_uso_endogeno", 12, 6)

# --- Modelos para el texto ------------------------------------------------------

lm_interaccion <- function(d, f) {
  m <- lm(f, data = d)
  b <- coef(m)
  V <- vcov(m)
  t_int <- "act_altaTRUE:tata_fuerteTRUE"
  l <- setNames(rep(0, length(b)), names(b))
  l[c("tata_fuerteTRUE", t_int)] <- 1
  est_alta <- sum(l * b)
  se_alta <- sqrt(drop(t(l) %*% V %*% l))
  cf <- coef(summary(m))
  tibble(
    n = nobs(m),
    tata_baja = b[["tata_fuerteTRUE"]], tata_baja_low = tata_baja - 1.96 * cf["tata_fuerteTRUE", 2], tata_baja_high = tata_baja + 1.96 * cf["tata_fuerteTRUE", 2], p_tata_baja = cf["tata_fuerteTRUE", 4],
    tata_alta = est_alta, tata_alta_low = est_alta - 1.96 * se_alta, tata_alta_high = est_alta + 1.96 * se_alta,
    interaccion = b[[t_int]], int_low = interaccion - 1.96 * cf[t_int, 2], int_high = interaccion + 1.96 * cf[t_int, 2], p_interaccion = cf[t_int, 4]
  )
}
modelos_lm <- largo %>%
  group_by(medida, rep) %>%
  group_modify(~ lm_interaccion(.x, log2(tpm + 1) ~ act_alta * tata_fuerte)) %>%
  ungroup() %>%
  mutate(across(c(tata_baja, tata_baja_low, tata_baja_high, tata_alta, tata_alta_low, tata_alta_high, interaccion, int_low, int_high), ~ 2^.x, .names = "{.col}_veces"))
save_tsv(modelos_lm, "modelo_log2_tpm_interaccion.tsv")
print(modelos_lm, n = Inf, width = Inf)

modelos_det <- largo %>%
  group_by(medida, rep) %>%
  group_modify(function(d, k) {
    m <- glm(I(tpm > umbral) ~ act_alta * tata_fuerte, family = binomial, data = d)
    cf <- coef(summary(m))
    t_int <- "act_altaTRUE:tata_fuerteTRUE"
    tibble(
      umbral_tpm = d$umbral[1],
      OR_tata_baja = exp(cf["tata_fuerteTRUE", 1]), OR_tata_baja_low = exp(cf["tata_fuerteTRUE", 1] - 1.96 * cf["tata_fuerteTRUE", 2]),
      OR_tata_baja_high = exp(cf["tata_fuerteTRUE", 1] + 1.96 * cf["tata_fuerteTRUE", 2]), p_tata_baja = cf["tata_fuerteTRUE", 4],
      OR_interaccion = exp(cf[t_int, 1]), p_interaccion = cf[t_int, 4]
    )
  }) %>%
  ungroup()
save_tsv(modelos_det, "modelo_deteccion_interaccion.tsv")
print(modelos_det, width = Inf)

# --- Figura completa ------------------------------------------------------------

fig <- (panel_a / panel_b) + plot_annotation(tag_levels = "A") & tag_theme
save_jpg_pdf(fig, "tata_uso_endogeno", 12, 12)

message("Salidas en ", out_dir)
