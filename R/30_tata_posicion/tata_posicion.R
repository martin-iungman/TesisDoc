# tata_posicion (ver docs/mapping_figuras.csv para el numero de figura
# vigente) - El efecto de la TATA-box sobre la actividad depende de su
# posicion canonica respecto del TSS. Version de figura del analisis
# exploratorio R/30_tata_posicion/exploratorio_tata_posicion.qmd; la busqueda
# (regla de EPD, puntajes, clases, perfil por ventana) vive en
# R/functions/tata_scan.R y la comparten ambos.
#
# Regla: matriz de Bucher (1990) con los pesos W de EPD, misma hebra, columna
# de referencia (columna 4, la segunda T de TATA) entre -31 y -26 del TSS,
# corte -8.17 - reproduce TATA_EPD exactamente. Solo promotores FP, cruzados
# con activity_stats_highconf.tsv; cada replica por separado. Referencia
# "sin motivo" = sin ninguna coincidencia >= corte en todo el fragmento, en
# ninguna hebra (salvo en tata_puntaje_corte, ver abajo).
#
# Genera dos figuras (cada panel en .jpg y .pdf):
#   tata_posicion (perfil por posicion + controles de composicion; figura
#   completa en tata_posicion.jpg/.pdf: A perfil, B modelo A/T, C permutadas):
#   - A, tata_posicion_perfil: HL (IC 95%) por ventana de 10 pb (alineadas
#      para que -33 a -24 contenga la ventana canonica), misma hebra / hebra
#      opuesta, por replica, puntos rellenos si BH < 0.05; abajo, proporcion
#      de promotores con coincidencia por ventana. Franja gris = -31 a -26.
#   - B, tata_controles_modelo_AT: coeficientes de cada clase vs. "sin
#      motivo" (puntos de percentil) sin covariable / + A/T (-50 a -1) / + A/T
#      proximal y del fragmento.
#   - C, tata_controles_permutadas: perfil de la misma hebra con la matriz
#      real (color) vs. 10 matrices con las columnas permutadas (gris; misma
#      composicion).
#   tata_puntaje (efecto segun el puntaje; el panel por bins fijos de puntaje
#   lo genera R/30_tata_posicion/tata_puntaje_bins.R):
#   - tata_puntaje_dosis.jpg: percentil de actividad medio (IC 95%) por bin de
#      1 unidad del mejor puntaje W (bins de los extremos unidos hasta tener
#      >= 30 promotores en cada replica), ventana canonica vs. ventana control
#      (-131 a -126). Linea vertical en -8.17.
#   - tata_puntaje_corte.jpg: HL de la TATA canonica para cortes de -10 a -4
#      cada 0.5 vs. una referencia FIJA (sin ninguna coincidencia >= -10 en
#      todo el fragmento, en ninguna hebra), para que el aumento del efecto
#      con el corte no venga de que la referencia cambie.
# Tablas (.tsv, en la carpeta de la figura que corresponde): ver la seccion
# "Tablas" al final.
#
# Requiere: data/raw/library.bed, BSgenome.Hsapiens.UCSC.hg38,
# data/processed/activity_stats_highconf.tsv y data/processed/prom_df.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork) # perfil: efecto arriba, proporcion abajo
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
source("R/functions/fig_paths.R")
source("R/functions/tata_scan.R")

select <- dplyr::select
filter <- dplyr::filter

dir_posicion <- fig_dir("tata_posicion")
dir_puntaje <- fig_dir("tata_puntaje")

CORTE <- TATA_CORTE_EPD
CANON <- TATA_CANON
CONTROL <- TATA_CANON - 100 # ventana control de la dosis-respuesta, mismo tamaño
CORTES_BARRIDO <- seq(-10, -4, by = 0.5)
CORTE_REF_BARRIDO <- min(CORTES_BARRIDO)
TOL <- TATA_TOL

col_misma <- "#216869"
col_opuesta <- "#14AFB2"
col_alto <- "#D6741F"
col_bajo <- "#7FB800"
th <- theme_pubclean(base_size = 18) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")
minus <- function(x) sub("^-", "−", format(x, trim = TRUE))
canon_lab <- paste0(minus(CANON[1]), " a ", minus(CANON[2]))
control_lab <- paste0(minus(CONTROL[1]), " a ", minus(CONTROL[2]))
franja <- annotate("rect", xmin = CANON[1], xmax = CANON[2], ymin = -Inf, ymax = Inf, fill = "grey88")
save_tsv <- function(x, dir, name) write_tsv(x, file.path(dir, name))
save_jpg_pdf <- function(plot, dir, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}

# --- Busqueda ---------------------------------------------------------------

lib <- tata_load_library()
code <- lib$code[!lib$has_n, ] # 1 fragmento con N fuera
sc <- tata_scores(code)
hits_10 <- tata_hits(sc$misma, sc$opuesta, CORTE_REF_BARRIDO) # coincidencias >= -10 (referencia fija del barrido de cortes)
hits <- filter(hits_10, score >= CORTE - TOL) # regla de EPD

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE) %>%
  inner_join(prom_df, by = c("seq_id", "name")) %>%
  filter(seq_id %in% rownames(code)) %>%
  group_by(rep) %>%
  mutate(act_pct = 100 * percent_rank(mean)) %>%
  ungroup()

clases <- tata_classify(hits, unique(data$seq_id))
data_clases <- inner_join(data, clases, by = "seq_id")

# --- perfil del efecto por posicion ----------------------------------------

col_hebra <- scale_color_manual(values = c("Misma hebra" = col_misma, "Hebra opuesta" = col_opuesta))
perfil <- tata_profile_effect(data, hits)

dodge <- position_dodge(width = 3)
panel_perfil_efecto <- ggplot(perfil, aes(win_center, estimate, color = hebra)) +
  franja +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = dodge, alpha = 0.7) +
  geom_line(aes(group = seg), position = dodge, alpha = 0.6) +
  geom_point(aes(shape = sig), position = dodge, size = 2.2, fill = "white") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "BH < 0,05", `FALSE` = "n.s."), name = NULL) +
  col_hebra +
  facet_grid(rep ~ .) +
  scale_x_continuous(breaks = seq(-240, 0, 20), labels = minus) +
  th +
  labs(x = NULL, y = "Diferencia de actividad vs.\nsin motivo (HL, IC 95%)", color = NULL)
panel_perfil_prop <- ggplot(perfil, aes(win_center, 100 * prop, color = hebra)) +
  franja +
  geom_line(aes(group = seg)) +
  geom_point(size = 1.3) +
  col_hebra +
  facet_grid(. ~ rep) +
  scale_x_continuous(breaks = seq(-200, 0, 50), labels = minus) +
  th +
  theme(legend.position = "none") +
  labs(x = "Posición de la coincidencia respecto del TSS (ventanas de 10 pb)", y = "% de promotores\ncon coincidencia")
panel_perfil <- panel_perfil_efecto / panel_perfil_prop + plot_layout(heights = c(2.2, 1))
save_jpg_pdf(panel_perfil, dir_posicion, "tata_posicion_perfil", 14, 10)

# --- matrices permutadas (misma hebra) -----------------------------------

perms <- tata_permutations(10, seed = 2024)
perfil_perm <- imap_dfr(perms, function(p, nm) {
  s <- tata_scores(code, TATA_W[, p])
  h <- tata_hits(s$misma, s$opuesta, CORTE)
  # referencia: sin coincidencias de la matriz permutada en ninguna hebra
  tata_profile_effect(data, filter(h, hebra == "Misma hebra"), ref_hits = h) %>%
    mutate(perm = nm, orden = paste(p, collapse = ","))
})
perfil_real_misma <- filter(perfil, hebra == "Misma hebra")

panel_permutadas <- ggplot(perfil_perm, aes(win_center, estimate)) +
  franja +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_line(aes(group = paste(perm, seg)), color = "grey65", alpha = 0.7) +
  geom_line(data = perfil_real_misma, aes(group = seg), color = col_misma, linewidth = 1.1) +
  geom_point(data = perfil_real_misma, color = col_misma, size = 2) +
  facet_wrap(~rep, ncol = 1) +
  scale_x_continuous(breaks = seq(-200, 0, 50), labels = minus) +
  th +
  labs(x = "Posición respecto del TSS (ventanas de 10 pb)", y = "Diferencia de actividad vs.\nsin motivo (HL)")
save_jpg_pdf(panel_permutadas, dir_posicion, "tata_controles_permutadas", 10, 8)

# --- dosis-respuesta por bins fijos de puntaje ----------------------------

# Bins de 1 unidad de puntaje W; se fusiona el bin con menos promotores (en la
# replica con menos) con su vecino mas chico hasta que todos tengan >= min_n
# en ambas replicas - en la practica une los extremos.
merge_bins <- function(bin, rep, min_n = 30) {
  repeat {
    tab <- table(factor(bin), rep)
    m <- apply(tab, 1, min)
    if (all(m >= min_n) || length(m) == 1) break
    lv <- as.numeric(names(m))
    i <- which.min(m)
    j <- if (i == 1) 2 else if (i == length(m)) i - 1 else if (m[i - 1] <= m[i + 1]) i - 1 else i + 1
    bin[bin == lv[i]] <- lv[j]
  }
  bin
}
best <- tibble(
  seq_id = rownames(code),
  canonica = tata_best_in_window(sc$misma, CANON[1], CANON[2]),
  control = tata_best_in_window(sc$misma, CONTROL[1], CONTROL[2])
) %>%
  pivot_longer(-seq_id, names_to = "ventana", values_to = "score") %>%
  mutate(ventana = recode(ventana, canonica = paste0("Canónica (", canon_lab, ")"), control = paste0("Control (", control_lab, ")")))
dosis <- data %>%
  select(rep, seq_id, mean, act_pct) %>%
  inner_join(best, by = "seq_id", relationship = "many-to-many") %>%
  group_by(ventana) %>%
  mutate(bin = merge_bins(floor(score), rep)) %>%
  ungroup()
dosis_bins <- dosis %>%
  group_by(rep, ventana, bin) %>%
  summarise(
    puntaje_min = min(score), puntaje_max = max(score), puntaje_medio = mean(score), n = n(),
    act_pct_media = mean(act_pct), se = sd(act_pct) / sqrt(n()), .groups = "drop"
  ) %>%
  mutate(conf.low = act_pct_media - 1.96 * se, conf.high = act_pct_media + 1.96 * se)
stopifnot(all(dosis_bins$n >= 30))

panel_dosis <- ggplot(dosis_bins, aes(puntaje_medio, act_pct_media, color = ventana)) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = CORTE, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_color_manual(values = setNames(c(col_misma, "grey55"), c(paste0("Canónica (", canon_lab, ")"), paste0("Control (", control_lab, ")")))) +
  scale_x_continuous(labels = minus) +
  facet_wrap(~rep, ncol = 1) +
  th +
  labs(x = "Mejor puntaje en la ventana (suma de pesos W; bins de 1 unidad)", y = "Percentil de actividad\n(media, IC 95%)", color = NULL)
save_jpg_pdf(panel_dosis, dir_puntaje, "tata_puntaje_dosis", 10, 8)

# --- control por A/T (modelo lineal) --------------------------------------

modelos_lab <- c("Sin covariable", "+ A/T (−50 a −1)", "+ A/T proximal y del fragmento")
at_prox <- tata_at_prox(code, -50, -1)
modelos <- data_clases %>%
  mutate(at_prox = at_prox[seq_id], at_frag = 1 - g_c, clase = relevel(clase, ref = "Sin motivo")) %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    bind_rows(
      broom::tidy(lm(act_pct ~ clase, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[1]),
      broom::tidy(lm(act_pct ~ clase + at_prox, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[2]),
      broom::tidy(lm(act_pct ~ clase + at_prox + at_frag, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[3])
    )
  }) %>%
  ungroup() %>%
  filter(str_detect(term, "^clase")) %>%
  mutate(
    clase = factor(str_remove(term, "^clase"), levels = rev(setdiff(TATA_CLASES, "Sin motivo"))),
    modelo = factor(modelo, levels = modelos_lab)
  )

panel_modelo_at <- ggplot(modelos, aes(estimate, clase, color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.6), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  scale_color_manual(values = setNames(c("grey55", col_opuesta, col_misma), modelos_lab)) +
  facet_wrap(~rep, ncol = 1) +
  th +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Diferencia de percentil de actividad vs. sin motivo", y = NULL, color = NULL)
save_jpg_pdf(panel_modelo_at, dir_posicion, "tata_controles_modelo_AT", 10, 8)

# Figura completa de tata_posicion: perfil (A) arriba; modelo A/T (B) y
# matrices permutadas (C) abajo
# (fila vacia intermedia para separar A de B y C)
figura_posicion <- wrap_elements(full = panel_perfil) / plot_spacer() /
  (wrap_elements(full = panel_modelo_at) | wrap_elements(full = panel_permutadas)) +
  plot_layout(heights = c(10, 0.8, 8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 26, face = "bold"))
save_jpg_pdf(figura_posicion, dir_posicion, "tata_posicion", 20, 18.8)

# --- efecto canonico segun el corte, referencia fija -----------------------

best_canon <- tata_best_in_window(sc$misma, CANON[1], CANON[2])
ref_fija <- setdiff(rownames(code), hits_10$seq_id)
barrido <- expand_grid(rep = sort(unique(data$rep)), corte = CORTES_BARRIDO) %>%
  mutate(res = map2(rep, corte, function(r, ct) {
    d <- filter(data, rep == r)
    x <- d$mean[d$seq_id %in% names(best_canon)[best_canon >= ct - TOL]]
    y <- d$mean[d$seq_id %in% ref_fija]
    tata_hl_vs(x, y) %>% mutate(n = length(x), n_ref = length(y))
  })) %>%
  unnest(res)
n_ref_corte <- barrido %>% distinct(rep, n_ref)

panel_corte <- ggplot(barrido, aes(corte, estimate, color = rep)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_vline(xintercept = CORTE, linetype = "dashed", color = "grey40") +
  annotate("text", x = CORTE, y = Inf, label = "−8,17 (EPD)", hjust = -0.08, vjust = 1.5, size = 5, color = "grey30") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.2)) +
  geom_text(aes(y = -Inf, label = n, vjust = if_else(rep == "Rep 1", -2.6, -1.1)), size = 3.8, show.legend = FALSE) +
  annotate("text", x = min(CORTES_BARRIDO) - 0.45, y = -Inf, label = "n canónicos", hjust = 0, vjust = -4.3, size = 3.8, color = "grey30") +
  scale_color_manual(values = c("Rep 1" = col_opuesta, "Rep 2" = col_misma)) +
  scale_x_continuous(breaks = seq(-10, -4, 1), labels = minus, expand = expansion(add = c(0.6, 0.3))) +
  scale_y_continuous(limits = c(-0.15, NA), expand = expansion(mult = c(0.02, 0.08))) +
  th +
  labs(
    x = "Corte (suma de pesos W)", y = "Diferencia de actividad vs.\nreferencia fija (HL, IC 95%)", color = NULL,
    subtitle = paste0("TATA canónica según el corte; referencia: sin coincidencias ≥ −10\n(n = ",
      paste0(n_ref_corte$n_ref, " en ", n_ref_corte$rep, collapse = ", "), ")")
  )
save_jpg_pdf(panel_corte, dir_puntaje, "tata_puntaje_corte", 10, 7)

# --- Tablas ----------------------------------------------------------------

# Datos de cada panel
save_tsv(perfil, dir_posicion, "efecto_por_ventana.tsv")
save_tsv(bind_rows(perfil_real_misma %>% mutate(perm = "Matriz real", orden = paste(seq_len(ncol(TATA_W)), collapse = ",")), perfil_perm), dir_posicion, "efecto_por_ventana_permutadas.tsv")
save_tsv(dosis_bins, dir_puntaje, "dosis_por_bin.tsv")
save_tsv(modelos, dir_posicion, "modelo_clase_AT.tsv")
save_tsv(barrido, dir_puntaje, "efecto_por_corte.tsv")

# clases x replica (vs. sin motivo)
res_clases <- tata_compare_classes(data_clases, "rep")
save_tsv(res_clases, dir_posicion, "clases_vs_sin_motivo.tsv")

# Spearman del mejor puntaje con la actividad, total y entre los que superan el corte
spearman <- data %>%
  select(rep, seq_id, mean) %>%
  inner_join(best, by = "seq_id", relationship = "many-to-many") %>%
  group_by(rep, ventana) %>%
  summarise(
    n = n(), rho = cor(score, mean, method = "spearman"),
    p = cor.test(score, mean, method = "spearman", exact = FALSE)$p.value,
    n_sobre_corte = sum(score >= CORTE - TOL),
    rho_sobre_corte = cor(score[score >= CORTE - TOL], mean[score >= CORTE - TOL], method = "spearman"),
    p_sobre_corte = cor.test(score[score >= CORTE - TOL], mean[score >= CORTE - TOL], method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  )
save_tsv(spearman, dir_puntaje, "spearman_puntaje_actividad.tsv")

# efecto de las clases por estrato: todos, angostos (tercil inferior de
# interquantile_width dentro de cada replica), con / sin isla CpG (CGI_frag).
# La referencia "sin motivo" es la del mismo estrato.
estratos <- bind_rows(
  data_clases %>% mutate(estrato = "Todos"),
  data_clases %>%
    filter(!is.na(interquantile_width)) %>%
    group_by(rep) %>%
    filter(ntile(interquantile_width, 3) == 1) %>%
    ungroup() %>%
    mutate(estrato = "Angostos (tercil inferior de interquantile_width)"),
  data_clases %>% filter(!is.na(CGI_frag)) %>% mutate(estrato = if_else(CGI_frag, "Con isla CpG (CGI_frag)", "Sin isla CpG (CGI_frag)"))
)
res_estratos <- tata_compare_classes(estratos, c("estrato", "rep"))
save_tsv(res_estratos, dir_posicion, "clases_por_estrato.tsv")

# % de promotores con alguna coincidencia >= corte en cualquier posicion y hebra
con_hit <- unique(hits$seq_id)
pct_hit <- bind_rows(
  tibble(universo = "Library (promotores FP)", n = nrow(code), n_con_coincidencia = sum(rownames(code) %in% con_hit)),
  data %>% group_by(universo = paste0("Con actividad (", rep, ")")) %>%
    summarise(n = n(), n_con_coincidencia = sum(seq_id %in% con_hit), .groups = "drop")
) %>% mutate(pct = round(100 * n_con_coincidencia / n, 1))
save_tsv(pct_hit, dir_posicion, "pct_promotores_con_coincidencia.tsv")

# pico de la hebra opuesta en la ventana canonica: cuantos tienen tambien la
# canonica en la misma hebra, y efecto de los que no (vs. sin motivo)
canon_ids <- hits %>% filter(hebra == "Misma hebra", between(pos, CANON[1], CANON[2])) %>% pull(seq_id) %>% unique()
op_ids <- hits %>% filter(hebra == "Hebra opuesta", between(pos, CANON[1], CANON[2])) %>% pull(seq_id) %>% unique()
op_solo <- setdiff(op_ids, canon_ids)
sin_motivo <- setdiff(rownames(code), con_hit)
opuesta <- bind_rows(
  tibble(rep = "Library (promotores FP)", n_opuesta = length(op_ids), n_tambien_canonica = sum(op_ids %in% canon_ids)),
  data %>% group_by(rep) %>%
    summarise(n_opuesta = sum(seq_id %in% op_ids), n_tambien_canonica = sum(seq_id %in% intersect(op_ids, canon_ids)), .groups = "drop")
) %>%
  mutate(prop_tambien_canonica = round(n_tambien_canonica / n_opuesta, 3)) %>%
  left_join(
    data %>% group_by(rep) %>%
      summarise(tata_hl_vs(mean[seq_id %in% op_solo], mean[seq_id %in% sin_motivo]),
        n_opuesta_sin_canonica = sum(seq_id %in% op_solo), n_ref = sum(seq_id %in% sin_motivo), .groups = "drop"),
    by = "rep"
  )
save_tsv(opuesta, dir_posicion, "hebra_opuesta_ventana_canonica.tsv")

message("Figuras y tablas guardadas en ", paste(c(dir_posicion, dir_puntaje), collapse = ", "))
