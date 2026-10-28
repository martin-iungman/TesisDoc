# EXPLORATORIO - el efecto de la TATA-box sobre el ruido (var_rank_sw)
# depende de su posicion canonica, de la secuencia del motivo (y no de la
# composicion A/T) y de la calidad del motivo (puntaje)? En particular:
# la TATA-box se asocia a mas ruido solo en promotores de alta actividad;
# eso es porque el efecto depende del nivel de actividad (lectura A) o
# porque las TATA de los promotores poco activos son mas debiles (lectura B)?
#
# Metodo (decidido en el Paso 0, ver mensaje de commit):
#   - var_rank_sw con add_noise_rank() sobre el set completo de cada replica
#     (igual que ruido_cgi): row_number(var) dentro de bins de 100 promotores
#     por actividad media. Despues se restringe a los promotores FP sin N
#     (los que tienen busqueda del motivo).
#   - AUC = P(ruido de los casos > ruido de la referencia), exacta
#     (Mann-Whitney, pROC) con IC 95% y p-valor de Wald de DeLong. NO se usa
#     noise_auc() (trapecio sobre umbrales 1-99): subestima ~0.01 porque
#     omite el tramo final de la curva y falla en estratos chicos (celdas
#     vacias -> NA -> tramos descartados). noise_auc() solo se reporta para
#     la TATA_EPD y CGI_frag globales, para comparar con ruido_cgi/ruido_tata.
#   - Actividad alta = percentil 75 de mean dentro de cada replica (no el
#     corte de noise_roc_by_expr(), que usa 0.75 * max(mean_sw) con las dos
#     replicas juntas).
#   - Busqueda del motivo: R/functions/tata_scan.R, igual que tata_posicion /
#     tata_puntaje (matriz de Bucher con pesos W de EPD, ambas hebras, corte
#     -8.17, ventanas de 10 pb alineadas en -33 a -24, clases, 10 matrices
#     permutadas con seed 2024, mejor puntaje en -31 a -26). La busqueda no
#     esta guardada en data/processed/ y tarda segundos, asi que se recalcula.
#   - Referencia "sin motivo" = sin ninguna coincidencia >= -8.17 en todo el
#     fragmento, en ninguna hebra (como en tata_posicion).
#
# Puntos (prefijo de los archivos de salida):
#   p0  AUC de referencia (TATA_EPD global, por actividad, por isla CpG),
#       con pROC y con noise_auc().
#   p1  AUC por ventana de posicion y hebra (formato de tata_posicion_perfil)
#       y por clase (canonica / solo no canonica / solo hebra opuesta; mas la
#       hebra opuesta en la ventana canonica sin la canonica).
#   p2  mismo perfil con las 10 matrices permutadas; modelo lineal del rango
#       de varianza con clase + A/T proximal (-50 a -1) + A/T del fragmento
#       (como tata_controles_modelo_AT).
#   p3  AUC por bin de puntaje del mejor match en -31 a -26 (sin corte,
#       bins de tata_puntaje_bins) vs. el resto de los promotores, ventana
#       canonica vs. control (-131 a -126).
#   p4  puntaje x actividad: canonicas con puntaje alto/bajo x actividad
#       alta/baja vs. sin motivo del mismo estrato de actividad;
#       distribucion del puntaje de las canonicas por actividad; Spearman
#       puntaje-ruido dentro de las canonicas por estrato.
#   p5  p4 dentro de los promotores con isla CpG (CGI_frag) y sin isla.
#   p6  que distingue a las TATA fuertes (canonicas de puntaje alto, p4) de
#       actividad baja de las de actividad alta: Fisher / Wilcoxon entre
#       ellas y, como control, interaccion actividad x TATA fuerte sobre
#       todos los promotores (variables del panorama de caracteristicas,
#       Polycomb en ReMap, elementos de secuencia y composicion; H1-H4).
#   p7  especificidad tisular continua (Gini de FANTOM5 y numero de muestras
#       con actividad): TATA fuertes vs. resto por actividad, y percentil de
#       actividad por quintil de Gini.
#   p8  AUC del ruido de las TATA fuertes vs. sin motivo por tercil de Gini
#       y estrato de actividad; ruido vs. Gini continuo.
#   p9  p1 y p2 (perfil por ventana, clases, permutadas, modelo A/T) dentro
#       de cada estrato de actividad (percentil 75 por replica).
#   p10 var_rank_sw ~ clase * actividad_alta (+ A/T): si el ajuste por A/T
#       reduce el efecto de la canonica por igual en ambos estratos.
#
# Union de Polycomb: data/processed/remap_tf_hits.tsv.
#
# Salidas en figures/exploratorio_tata_ruido/ (sin fig_dir: no es figura de
# la tesis). Requiere: data/raw/library.bed, BSgenome.Hsapiens.UCSC.hg38,
# data/processed/activity_stats_highconf.tsv y data/processed/prom_df.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
library(pROC)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
source("R/functions/plot_helpers.R")
source("R/functions/tata_scan.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename # S4Vectors la enmascara (build_tidy_features la usa)

out_dir <- "figures/exploratorio_tata_ruido"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

CORTE <- TATA_CORTE_EPD
CANON <- TATA_CANON
CONTROL <- TATA_CANON - 100
UMBRAL_ALTO <- -5 # puntaje "alto" de p4 si da n razonables (ver abajo)
LIMITES <- c(-Inf, -16, -14, -12, -10, CORTE, -7, -6, -5, Inf)
ETIQUETAS <- c("< −16", "−16 a −14", "−14 a −12", "−12 a −10", "−10 a −8.17", "−8.17 a −7", "−7 a −6", "−6 a −5", "≥ −5")
N_MIN <- 20

col_misma <- "#216869"
col_opuesta <- "#14AFB2"
clr_rep <- c("Rep 1" = "#0D2C54", "Rep 2" = "#AD343E")
th <- theme_pubclean(base_size = 18) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")
minus <- function(x) sub("^-", "−", format(x, trim = TRUE))
franja <- annotate("rect", xmin = CANON[1], xmax = CANON[2], ymin = -Inf, ymax = Inf, fill = "grey88")
save_tsv <- function(x, name) write_tsv(x, file.path(out_dir, name))
save_jpg <- function(plot, name, width, height) ggsave(file.path(out_dir, paste0(name, ".jpg")), plot, width = width, height = height, units = "in")
options(width = 220)
show <- function(x, titulo) {
  cat("\n==", titulo, "==\n")
  print(x %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)
}

# AUC de x (casos) vs. y (referencia): P(x > y), exacta, IC y p de DeLong
auc_vs <- function(x, y, min_n = 10) {
  if (length(x) < min_n || length(y) < min_n) {
    return(tibble(n = length(x), n_ref = length(y), AUC = NA_real_, conf.low = NA_real_, conf.high = NA_real_, se = NA_real_, p = NA_real_))
  }
  r <- pROC::roc(controls = y, cases = x, direction = "<", ci = TRUE, quiet = TRUE)
  auc <- as.numeric(r$auc)
  se <- sqrt(pROC::var(r, method = "delong"))
  tibble(n = length(x), n_ref = length(y), AUC = auc, conf.low = r$ci[1], conf.high = r$ci[3], se = se, p = 2 * pnorm(-abs((auc - 0.5) / se)))
}
# diferencia de dos AUC de muestras independientes (z de DeLong)
auc_diff_p <- function(a1, se1, a2, se2) 2 * pnorm(-abs((a1 - a2) / sqrt(se1^2 + se2^2)))

# --- Datos y busqueda --------------------------------------------------------

lib <- tata_load_library()
code <- lib$code[!lib$has_n, ] # 1 fragmento con N fuera
sc <- tata_scores(code)
hits <- tata_hits(sc$misma, sc$opuesta, CORTE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  filter(seq_id %in% rownames(code)) %>%
  group_by(rep) %>%
  mutate(act_pct = 100 * percent_rank(mean), actividad = if_else(act_pct >= 75, "Alta", "Baja")) %>%
  ungroup() %>%
  inner_join(tata_classify(hits, rownames(code)), by = "seq_id") %>%
  mutate(
    score_canon = tata_best_in_window(sc$misma, CANON[1], CANON[2])[seq_id],
    score_control = tata_best_in_window(sc$misma, CONTROL[1], CONTROL[2])[seq_id]
  )
stopifnot(all((data$clase == "Canónica") == data$TATA_EPD))
show(count(data, rep, actividad, clase) %>% pivot_wider(names_from = clase, values_from = n), "n por clase y actividad")

sin_motivo <- function(d) d$var_rank_sw[d$clase == "Sin motivo"]

# --- p0: valores de referencia ------------------------------------------------

auc_feature <- function(d, feature) {
  d <- d[!is.na(d[[feature]]), ]
  auc_vs(d$var_rank_sw[d[[feature]]], d$var_rank_sw[!d[[feature]]]) %>%
    mutate(AUC_noise_auc = noise_auc(roc_curve(d, feature)))
}
p0 <- bind_rows(
  data %>% group_by(rep) %>% group_modify(~ auc_feature(.x, "CGI_frag")) %>% mutate(feature = "CGI_frag", estrato = "Todos"),
  data %>% group_by(rep) %>% group_modify(~ auc_feature(.x, "TATA_EPD")) %>% mutate(feature = "TATA_EPD", estrato = "Todos"),
  data %>% group_by(rep, estrato = paste("Actividad", tolower(actividad))) %>% group_modify(~ auc_feature(.x, "TATA_EPD")) %>% mutate(feature = "TATA_EPD"),
  data %>% filter(!is.na(CGI_frag)) %>% group_by(rep, estrato = if_else(CGI_frag, "Con isla", "Sin isla")) %>%
    group_modify(~ auc_feature(.x, "TATA_EPD")) %>% mutate(feature = "TATA_EPD"),
  data %>% filter(!is.na(CGI_frag)) %>% group_by(rep, estrato = paste(if_else(CGI_frag, "Con isla", "Sin isla"), "- actividad", tolower(actividad))) %>%
    group_modify(~ auc_feature(.x, "TATA_EPD")) %>% mutate(feature = "TATA_EPD")
) %>%
  ungroup() %>%
  select(feature, estrato, rep, everything())
save_tsv(p0, "p0_auc_referencia.tsv")
show(p0, "p0 AUC de referencia (TATA_EPD vs. el resto)")

# --- p1: perfil por posicion ---------------------------------------------------

# Como tata_profile_effect(), con AUC del rango de varianza en vez de HL
profile_auc <- function(d_all, hits_u, ref_hits = hits_u, min_n = N_MIN) {
  hit_ids <- unique(ref_hits$seq_id)
  wins <- hits_u %>% distinct(seq_id, hebra, win)
  d_all %>%
    group_by(rep) %>%
    group_modify(function(d, k) {
      ref <- d$var_rank_sw[!d$seq_id %in% hit_ids]
      n_tot <- nrow(d)
      d %>%
        inner_join(wins, by = "seq_id", relationship = "many-to-many") %>%
        group_by(hebra, win) %>%
        filter(n() >= min_n) %>%
        summarise(auc_vs(var_rank_sw, ref), prop = n() / n_tot, .groups = "drop") %>%
        mutate(p_adj = p.adjust(p, "BH")) %>%
        arrange(hebra, win) %>%
        group_by(hebra) %>%
        mutate(seg = paste(hebra, cumsum(c(1, diff(win) != 10)))) %>%
        ungroup()
    }) %>%
    ungroup() %>%
    mutate(win_center = win + 4.5, sig = p_adj < 0.05)
}

perfil <- profile_auc(data, hits)
save_tsv(perfil, "p1_auc_por_ventana.tsv")
show(perfil %>% filter(between(win, -63, 6)) %>% select(rep, hebra, win, n, AUC, conf.low, conf.high, p_adj), "p1 AUC por ventana (-63 a +6)")

col_hebra <- scale_color_manual(values = c("Misma hebra" = col_misma, "Hebra opuesta" = col_opuesta))
dodge <- position_dodge(width = 3)
panel_perfil_auc <- ggplot(perfil, aes(win_center, AUC, color = hebra)) +
  franja +
  geom_hline(yintercept = 0.5, color = "grey50") +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = dodge, alpha = 0.7) +
  geom_line(aes(group = seg), position = dodge, alpha = 0.6) +
  geom_point(aes(shape = sig), position = dodge, size = 2.2, fill = "white") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "BH < 0.05", `FALSE` = "n.s."), name = NULL) +
  col_hebra +
  facet_grid(rep ~ .) +
  scale_x_continuous(breaks = seq(-240, 0, 20), labels = minus) +
  th +
  labs(x = NULL, y = "AUC del ruido vs.\nsin motivo (IC 95%)", color = NULL)
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
save_jpg(panel_perfil_auc / panel_perfil_prop + plot_layout(heights = c(2.2, 1)), "p1_perfil_posicion_auc", 14, 10)

# Clases vs. sin motivo, y la hebra opuesta en la ventana canonica sin la
# canonica (el pico de la hebra opuesta en -33 a -24 es en su mayoria la
# misma TATA leida en complemento reverso)
canon_ids <- data$seq_id[data$clase == "Canónica"] %>% unique()
op_canon_solo <- hits %>%
  filter(hebra == "Hebra opuesta", between(pos, CANON[1], CANON[2])) %>%
  pull(seq_id) %>% unique() %>% setdiff(canon_ids)
clases_auc <- data %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    ref <- sin_motivo(d)
    bind_rows(
      map_dfr(setdiff(TATA_CLASES, "Sin motivo"), ~ auc_vs(d$var_rank_sw[d$clase == .x], ref) %>% mutate(clase = .x)),
      auc_vs(d$var_rank_sw[d$seq_id %in% op_canon_solo], ref) %>% mutate(clase = "Hebra opuesta en −31 a −26, sin canónica")
    ) %>% mutate(p_adj = p.adjust(p, "BH"))
  }) %>%
  ungroup() %>%
  mutate(clase = fct_rev(fct_inorder(clase)))
save_tsv(clases_auc, "p1_auc_por_clase.tsv")
show(clases_auc %>% select(rep, clase, n, n_ref, AUC, conf.low, conf.high, p_adj), "p1 AUC por clase vs. sin motivo")

forest <- function(df, y, color = "rep", xlab = "AUC del ruido vs. sin motivo (IC 95%)") {
  ggplot(df, aes(AUC, .data[[y]], color = .data[[color]])) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.6), linewidth = 1) +
    geom_point(position = position_dodge(width = 0.6), size = 2.8) +
    scale_color_manual(values = clr_rep) +
    th +
    labs(x = xlab, y = NULL, color = NULL)
}
save_jpg(forest(clases_auc, "clase"), "p1_clases_auc", 10, 5)

# --- p2: matrices permutadas y modelo A/T ------------------------------------

perms <- tata_permutations(10, seed = 2024)
perfil_perm <- imap_dfr(perms, function(p, nm) {
  s <- tata_scores(code, TATA_W[, p])
  h <- tata_hits(s$misma, s$opuesta, CORTE)
  profile_auc(data, filter(h, hebra == "Misma hebra"), ref_hits = h) %>%
    mutate(perm = nm, orden = paste(p, collapse = ","))
})
perfil_real_misma <- filter(perfil, hebra == "Misma hebra")
save_tsv(bind_rows(perfil_real_misma %>% mutate(perm = "Matriz real"), perfil_perm), "p2_auc_por_ventana_permutadas.tsv")
perm_canon <- bind_rows(perfil_real_misma %>% mutate(perm = "Matriz real"), perfil_perm) %>%
  filter(win == TATA_WIN_ORIGIN) %>%
  select(rep, perm, n, AUC, conf.low, conf.high, p_adj)
save_tsv(perm_canon, "p2_permutadas_ventana_canonica.tsv")
show(perm_canon, "p2 ventana -33 a -24 (misma hebra), matriz real vs. permutadas")
show(perfil_perm %>% group_by(rep) %>% summarise(max_AUC = max(AUC, na.rm = TRUE), n_ventanas_BH_auc_mayor = sum(sig & AUC > 0.5, na.rm = TRUE), n_ventanas = n()), "p2 permutadas: maximo AUC en cualquier ventana")

panel_perm <- ggplot(perfil_perm, aes(win_center, AUC)) +
  franja +
  geom_hline(yintercept = 0.5, color = "grey50") +
  geom_line(aes(group = paste(perm, seg)), color = "grey65", alpha = 0.7) +
  geom_line(data = perfil_real_misma, aes(group = seg), color = col_misma, linewidth = 1.1) +
  geom_point(data = perfil_real_misma, color = col_misma, size = 2) +
  facet_wrap(~rep, ncol = 1) +
  scale_x_continuous(breaks = seq(-200, 0, 50), labels = minus) +
  th +
  labs(x = "Posición respecto del TSS (ventanas de 10 pb)", y = "AUC del ruido vs.\nsin motivo")
save_jpg(panel_perm, "p2_permutadas_auc", 10, 8)

modelos_lab <- c("Sin covariable", "+ A/T (−50 a −1)", "+ A/T proximal y del fragmento")
at_prox <- tata_at_prox(code, -50, -1)
modelos <- data %>%
  mutate(at_prox = at_prox[seq_id], at_frag = 1 - g_c, clase = relevel(clase, ref = "Sin motivo")) %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    bind_rows(
      broom::tidy(lm(var_rank_sw ~ clase, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[1]),
      broom::tidy(lm(var_rank_sw ~ clase + at_prox, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[2]),
      broom::tidy(lm(var_rank_sw ~ clase + at_prox + at_frag, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[3])
    )
  }) %>%
  ungroup() %>%
  filter(str_detect(term, "^clase|^at_")) %>%
  mutate(modelo = factor(modelo, levels = modelos_lab))
save_tsv(modelos, "p2_modelo_clase_AT.tsv")
show(modelos %>% select(rep, modelo, term, estimate, conf.low, conf.high, p.value), "p2 modelo lineal var_rank_sw ~ clase + A/T")

panel_modelo <- modelos %>%
  filter(str_detect(term, "^clase")) %>%
  mutate(clase = factor(str_remove(term, "^clase"), levels = rev(setdiff(TATA_CLASES, "Sin motivo")))) %>%
  ggplot(aes(estimate, clase, color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.6), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  scale_color_manual(values = setNames(c("grey55", col_opuesta, col_misma), modelos_lab)) +
  facet_wrap(~rep, ncol = 1) +
  th +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Diferencia de rango de varianza vs. sin motivo", y = NULL, color = NULL)
save_jpg(panel_modelo, "p2_modelo_AT", 10, 8)

# --- p3: AUC segun el puntaje -------------------------------------------------

# Cada bin vs. el resto de los promotores (no hay referencia "sin motivo"
# limpia por debajo del corte: un promotor con puntaje canonico bajo puede
# tener coincidencias en otra parte). Mismos bins que tata_puntaje_bins.
puntaje <- data %>%
  select(rep, seq_id, var_rank_sw, canonica = score_canon, control = score_control) %>%
  pivot_longer(c(canonica, control), names_to = "ventana", values_to = "score") %>%
  mutate(
    ventana = if_else(ventana == "canonica", "Canónica (−31 a −26)", "Control (−131 a −126)"),
    bin = cut(score + TATA_TOL, breaks = LIMITES, labels = ETIQUETAS, right = FALSE)
  )
puntaje_bins <- puntaje %>%
  group_by(rep, ventana) %>%
  group_modify(function(d, k) {
    map_dfr(levels(d$bin), function(b) {
      en <- d$bin == b
      auc_vs(d$var_rank_sw[en], d$var_rank_sw[!en]) %>%
        mutate(bin = b, puntaje_mediana = median(d$score[en]), ruido_medio = mean(d$var_rank_sw[en]))
    })
  }) %>%
  ungroup() %>%
  mutate(bin = factor(bin, levels = ETIQUETAS))
save_tsv(puntaje_bins, "p3_auc_por_bin_puntaje.tsv")
show(puntaje_bins %>% select(rep, ventana, bin, n, puntaje_mediana, ruido_medio, AUC, conf.low, conf.high, p), "p3 AUC por bin de puntaje (vs. resto)")
show(puntaje %>% group_by(rep, ventana) %>% summarise(rho = cor(score, var_rank_sw, method = "spearman"), rho_sobre_corte = cor(score[score >= CORTE - TATA_TOL], var_rank_sw[score >= CORTE - TATA_TOL], method = "spearman"), n_sobre_corte = sum(score >= CORTE - TATA_TOL), .groups = "drop"), "p3 Spearman puntaje-ruido")

panel_puntaje <- ggplot(filter(puntaje_bins, n >= N_MIN), aes(puntaje_mediana, AUC, color = ventana)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = CORTE, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  geom_line() +
  geom_point(size = 3) +
  annotate("text", x = CORTE, y = Inf, label = "−8.17 (EPD)", hjust = 1.05, vjust = 1.4, size = 5, color = "grey30") +
  scale_color_manual(values = c("Canónica (−31 a −26)" = col_misma, "Control (−131 a −126)" = "grey55")) +
  scale_x_continuous(labels = function(x) sub("^-", "−", x), expand = expansion(add = c(1, 0.8))) +
  facet_wrap(~rep) +
  th +
  labs(x = "Mejor puntaje TATA en la ventana (mediana del bin)", y = "AUC del ruido\nvs. resto (IC 95%)", color = NULL)
save_jpg(panel_puntaje, "p3_puntaje_bins_auc", 14, 6.5)

# --- p4/p5: puntaje x actividad ------------------------------------------------

canon <- filter(data, clase == "Canónica")
mediana_puntaje <- canon %>% group_by(rep) %>% summarise(mediana = median(score_canon), n_ge_umbral = sum(score_canon >= UMBRAL_ALTO - TATA_TOL), n = n())
show(mediana_puntaje, "p4 mediana del puntaje de las canonicas y n >= -5")
show(canon %>% count(rep, actividad, alto = score_canon >= UMBRAL_ALTO - TATA_TOL), "p4 n por celda con el umbral -5")

# Distribucion del puntaje de las canonicas por actividad
dist_puntaje <- canon %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    w <- wilcox.test(score_canon ~ actividad, data = d, conf.int = TRUE, exact = FALSE)
    d %>%
      group_by(actividad) %>%
      summarise(n = n(), mediana = median(score_canon), q1 = quantile(score_canon, 0.25), q3 = quantile(score_canon, 0.75), prop_ge_umbral = mean(score_canon >= UMBRAL_ALTO - TATA_TOL), .groups = "drop") %>%
      mutate(HL_alta_menos_baja = unname(w$estimate), p_wilcoxon = w$p.value)
  }) %>%
  ungroup()
save_tsv(dist_puntaje, "p4_distribucion_puntaje_por_actividad.tsv")
show(dist_puntaje, "p4 puntaje de las canonicas por actividad")

panel_dist <- ggplot(canon, aes(actividad, score_canon, fill = actividad)) +
  geom_violin(alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = UMBRAL_ALTO, linetype = "dashed", color = "grey40") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.35, size = 5) +
  scale_fill_manual(values = c(Alta = col_misma, Baja = "grey60")) +
  scale_y_continuous(labels = function(x) sub("^-", "−", x)) +
  facet_wrap(~rep) +
  th +
  theme(legend.position = "none") +
  labs(x = "Actividad (percentil 75 por réplica)", y = "Puntaje de la TATA canónica")
save_jpg(panel_dist, "p4_distribucion_puntaje", 9, 6)

# AUC de cada celda (puntaje x actividad) vs. sin motivo del mismo estrato de
# actividad; ademas la canonica completa por estrato, la diferencia alto - bajo
# dentro de cada estrato (z de DeLong) y Spearman puntaje-ruido dentro de las
# canonicas de cada estrato (sin dicotomizar).
celdas_auc <- function(d_all, corte_puntaje, grupo = "Todos") {
  d_all %>%
    group_by(rep, actividad) %>%
    group_modify(function(d, k) {
      ref <- sin_motivo(d)
      cn <- d[d$clase == "Canónica", ]
      alto <- cn$score_canon >= corte_puntaje - TATA_TOL
      res <- bind_rows(
        auc_vs(cn$var_rank_sw, ref) %>% mutate(puntaje = "Todas"),
        auc_vs(cn$var_rank_sw[alto], ref) %>% mutate(puntaje = "Alto"),
        auc_vs(cn$var_rank_sw[!alto], ref) %>% mutate(puntaje = "Bajo")
      )
      rho <- if (nrow(cn) >= 10) cor.test(cn$score_canon, cn$var_rank_sw, method = "spearman", exact = FALSE) else list(estimate = NA, p.value = NA)
      res %>% mutate(
        p_alto_vs_bajo = auc_diff_p(AUC[puntaje == "Alto"], se[puntaje == "Alto"], AUC[puntaje == "Bajo"], se[puntaje == "Bajo"]),
        rho_puntaje_ruido = unname(rho$estimate), p_rho = rho$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(grupo = grupo, corte_puntaje = corte_puntaje, puntaje = factor(puntaje, levels = c("Todas", "Bajo", "Alto")))
}

# Corte del puntaje: mediana de las canonicas (pooled de ambas replicas,
# redondeada; las dos medianas son practicamente iguales) - el umbral -5 deja
# pocas canonicas de puntaje alto en baja actividad (ver la tabla de n arriba);
# -5 se reporta igual como sensibilidad.
CORTE_MEDIANA <- round(median(canon$score_canon), 2)
cat("\nCorte de puntaje (mediana de las canonicas):", CORTE_MEDIANA, "\n")

p4 <- bind_rows(
  celdas_auc(data, CORTE_MEDIANA),
  celdas_auc(data, UMBRAL_ALTO),
  celdas_auc(filter(data, CGI_frag %in% TRUE), CORTE_MEDIANA, "Con isla CpG"),
  celdas_auc(filter(data, CGI_frag %in% FALSE), CORTE_MEDIANA, "Sin isla CpG"),
  celdas_auc(filter(data, CGI_frag %in% TRUE), UMBRAL_ALTO, "Con isla CpG"),
  celdas_auc(filter(data, CGI_frag %in% FALSE), UMBRAL_ALTO, "Sin isla CpG")
) %>%
  mutate(actividad = factor(actividad, levels = c("Baja", "Alta")))
save_tsv(p4, "p4_p5_auc_puntaje_x_actividad.tsv")
show(p4 %>% select(grupo, corte_puntaje, rep, actividad, puntaje, n, n_ref, AUC, conf.low, conf.high, p, p_alto_vs_bajo, rho_puntaje_ruido, p_rho), "p4/p5 AUC puntaje x actividad vs. sin motivo del estrato")

# Sensibilidad: referencia = todos los no canonicos del estrato (como TATA_EPD vs. resto)
p4_resto <- data %>%
  group_by(rep, actividad) %>%
  group_modify(function(d, k) {
    ref <- d$var_rank_sw[d$clase != "Canónica"]
    cn <- d[d$clase == "Canónica", ]
    alto <- cn$score_canon >= CORTE_MEDIANA - TATA_TOL
    bind_rows(auc_vs(cn$var_rank_sw[alto], ref) %>% mutate(puntaje = "Alto"), auc_vs(cn$var_rank_sw[!alto], ref) %>% mutate(puntaje = "Bajo"))
  }) %>%
  ungroup()
save_tsv(p4_resto, "p4_auc_puntaje_x_actividad_ref_resto.tsv")
show(p4_resto, "p4 sensibilidad: referencia = no canonicos del estrato")

panel_celdas <- function(df) {
  df %>%
    mutate(lab = paste0("n=", n)) %>%
    ggplot(aes(actividad, AUC, color = puntaje)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5), linewidth = 1) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_text(aes(y = conf.high, label = lab), vjust = -0.6, size = 4, position = position_dodge(width = 0.5), show.legend = FALSE) +
    scale_color_manual(values = c(Todas = "grey55", Bajo = col_opuesta, Alto = col_misma)) +
    facet_grid(grupo ~ rep) +
    th +
    labs(x = "Actividad (percentil 75 por réplica)", y = "AUC del ruido vs. sin motivo\ndel mismo estrato (IC 95%)", color = "Puntaje de la TATA canónica")
}
save_jpg(panel_celdas(filter(p4, corte_puntaje == CORTE_MEDIANA, grupo == "Todos")), "p4_puntaje_x_actividad", 10, 6.5)
save_jpg(panel_celdas(filter(p4, corte_puntaje == UMBRAL_ALTO, grupo == "Todos")), "p4_puntaje_x_actividad_corte_menos5", 10, 6.5)
save_jpg(panel_celdas(filter(p4, corte_puntaje == CORTE_MEDIANA, grupo != "Todos")), "p5_puntaje_x_actividad_por_isla", 10, 10)

# --- p6: que distingue a las TATA fuertes de baja actividad ---------------------

# TATA fuerte = canonica con puntaje >= mediana de las canonicas (las "Alto"
# de p4). Para cada variable: TATA fuertes de actividad alta vs. baja
# (Fisher / Wilcoxon, BH sobre todas las variables dentro de cada replica y
# conjunto) y, como control, variable ~ actividad_alta * tata_fuerte sobre
# todos los promotores (logistica / lineal): OR (o diferencia de medias) de
# alta vs. baja en las TATA fuertes y en el resto, y p de la interaccion.
# Variables: definiciones de build_tidy_features() (panorama de
# caracteristicas; terciles calculados sobre todos los promotores de cada
# replica), columnas de prom_df y union de Polycomb en ReMap (cualquier
# muestra). "Polycomb (cualquiera)" = RING1, CBX8, JARID2, PHF19, EZH2 o
# SUZ12; RNF2 (RING1B) se agrega aparte, fuera de "cualquiera".

PC_TFS <- c("RING1", "CBX8", "JARID2", "PHF19", "EZH2", "SUZ12")
remap_pc <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE) %>%
  filter(TF %in% c(PC_TFS, "RNF2")) %>%
  distinct(seq_id, TF)

hipotesis_vars <- tribble(
  ~variable, ~hipotesis,
  "Alta especificidad tisular", "H1", "Baja especificidad tisular", "H1",
  "Promotores angostos", "H1", "Promotores anchos", "H1",
  "No detectado (FANTOM5)", "H1", "Alta actividad en HEK293", "H1",
  "Alta accesibilidad de cromatina (DNase-seq)", "H1",
  "Alta conservación (16 a -50pb)", "H1", "Alta conservación (-50 a -150)", "H1", "Alta conservación (-150 a -235)", "H1",
  "Sin actividad en ratón", "H1", "Promotor unidireccional", "H1", "Sin módulo cis-regulatorio", "H1",
  "Polycomb (cualquiera)", "H1",
  "ReMap: RING1", "H1", "ReMap: CBX8", "H1", "ReMap: JARID2", "H1", "ReMap: PHF19", "H1",
  "ReMap: EZH2", "H1", "ReMap: SUZ12", "H1", "ReMap: RNF2", "H1",
  "TSS no canónico", "H2", "E-box CAGGTR", "H2", "Tracto poli(dA:dT) proximal", "H2", "Elementos transponibles (TE_overlap > 0)", "H2",
  "PyPu en -1/0", "H3", "Inr (EPD)", "H3", "TSS fuerte (Inr fuerte)", "H3",
  "Isla CpG (CGI_frag)", "H4", "Contenido G+C", "H4", "CpG o/e", "H4",
  "CCAAT-box", "Otros", "GC-box", "Otros", "YY1", "Otros", "CRE (CREB/ATF)", "Otros"
) %>% mutate(hipotesis = factor(hipotesis, levels = c("H1", "H2", "H3", "H4", "Otros")))
vars_cont <- c("Contenido G+C", "CpG o/e")

endo_vars <- hipotesis_vars$variable[hipotesis_vars$hipotesis == "H1" & !str_detect(hipotesis_vars$variable, "Polycomb|ReMap")]
tidy_endo <- build_tidy_features(data, keep = "seq_id") %>%
  mutate(rep = as.character(rep), across(where(is.factor), ~ .x == "TRUE")) %>%
  select(seq_id, rep, all_of(endo_vars))

feat <- data %>%
  transmute(
    seq_id, rep, act_alta = actividad == "Alta",
    tata_fuerte = clase == "Canónica" & score_canon >= CORTE_MEDIANA - TATA_TOL,
    score_canon,
    `Contenido G+C` = g_c, `CpG o/e` = cpg_oe,
    `Isla CpG (CGI_frag)` = as.logical(CGI_frag),
    `Elementos transponibles (TE_overlap > 0)` = replace_na(TE_overlap > 0, FALSE),
    `E-box CAGGTR` = as.logical(Ebox_CAGGTR), `Tracto poli(dA:dT) proximal` = as.logical(polyAT_prox),
    `TSS no canónico` = as.logical(other_TSS), `PyPu en -1/0` = as.logical(PyPu_TSS),
    `Inr (EPD)` = as.logical(INR_EPD), `TSS fuerte (Inr fuerte)` = as.logical(INR_strong_TSS),
    `CCAAT-box` = as.logical(CCAAT_EPD), `GC-box` = as.logical(GCbox_EPD),
    YY1 = as.logical(YY1_motif), `CRE (CREB/ATF)` = as.logical(CRE_motif),
    `Polycomb (cualquiera)` = seq_id %in% remap_pc$seq_id[remap_pc$TF %in% PC_TFS]
  )
for (tf in c(PC_TFS, "RNF2")) feat[[paste("ReMap:", tf)]] <- feat$seq_id %in% remap_pc$seq_id[remap_pc$TF == tf]
feat <- left_join(feat, tidy_endo, by = c("seq_id", "rep"))

# Concordancia del grupo de actividad de las TATA fuertes entre replicas
conc <- feat %>%
  filter(tata_fuerte) %>%
  select(seq_id, rep, act_alta) %>%
  pivot_wider(names_from = rep, values_from = act_alta)
conc_tab <- conc %>%
  summarise(
    n_rep1 = sum(!is.na(`Rep 1`)), n_rep2 = sum(!is.na(`Rep 2`)), n_ambas = sum(!is.na(`Rep 1`) & !is.na(`Rep 2`)),
    concordantes_alta = sum(`Rep 1` & `Rep 2`, na.rm = TRUE), concordantes_baja = sum(!`Rep 1` & !`Rep 2`, na.rm = TRUE),
    discordantes = sum(`Rep 1` != `Rep 2`, na.rm = TRUE)
  ) %>%
  mutate(prop_concordantes = (concordantes_alta + concordantes_baja) / n_ambas)
save_tsv(conc_tab, "p6_concordancia_actividad_tata_fuertes.tsv")
show(conc_tab, "p6 concordancia del grupo de actividad de las TATA fuertes")
show(count(feat, rep, tata_fuerte, act_alta), "p6 n por grupo")

# Promotores (todos, no solo TATA) con el mismo grupo de actividad en ambas replicas
ids_conc <- feat %>%
  select(seq_id, rep, act_alta) %>%
  pivot_wider(names_from = rep, values_from = act_alta) %>%
  filter(!is.na(`Rep 1`), !is.na(`Rep 2`), `Rep 1` == `Rep 2`) %>%
  pull(seq_id)

or_glm <- function(dd) {
  tab <- table(factor(dd$act_alta, c(TRUE, FALSE)), factor(dd$valor, c(TRUE, FALSE)))
  if (any(tab == 0)) {
    return(tibble(OR = NA_real_, OR_low = NA_real_, OR_high = NA_real_))
  }
  b <- coef(summary(glm(valor ~ act_alta, family = binomial, data = dd)))["act_altaTRUE", ]
  tibble(OR = exp(b[1]), OR_low = exp(b[1] - 1.96 * b[2]), OR_high = exp(b[1] + 1.96 * b[2]))
}
compara_variable <- function(d, v) {
  d <- d %>% mutate(valor = .data[[v]]) %>% filter(!is.na(valor))
  s <- filter(d, tata_fuerte)
  r <- filter(d, !tata_fuerte)
  if (v %in% vars_cont) {
    w <- wilcox.test(s$valor[s$act_alta], s$valor[!s$act_alta], exact = FALSE)
    lm_dif <- function(dd) coef(lm(valor ~ act_alta, data = dd))[["act_altaTRUE"]]
    tibble(
      tipo = "continua",
      fuerte_alta = median(s$valor[s$act_alta]), fuerte_baja = median(s$valor[!s$act_alta]),
      dif_medianas_fuertes = fuerte_alta - fuerte_baja, p = w$p.value,
      dif_medias_fuertes = lm_dif(s), dif_medias_resto = lm_dif(r),
      p_interaccion = coef(summary(lm(valor ~ act_alta * tata_fuerte, data = d)))["act_altaTRUE:tata_fuerteTRUE", 4]
    )
  } else {
    ft <- fisher.test(table(factor(s$act_alta, c(TRUE, FALSE)), factor(s$valor, c(TRUE, FALSE))))
    m <- glm(valor ~ act_alta * tata_fuerte, family = binomial, data = d)
    cf <- coef(summary(m))
    tibble(
      tipo = "binaria",
      fuerte_alta = mean(s$valor[s$act_alta]), fuerte_baja = mean(s$valor[!s$act_alta]),
      resto_alta = mean(r$valor[r$act_alta]), resto_baja = mean(r$valor[!r$act_alta]),
      OR_fisher = unname(ft$estimate), OR_fisher_low = ft$conf.int[1], OR_fisher_high = ft$conf.int[2], p = ft$p.value,
      p_interaccion = if ("act_altaTRUE:tata_fuerteTRUE" %in% rownames(cf)) cf["act_altaTRUE:tata_fuerteTRUE", 4] else NA_real_
    ) %>%
      bind_cols(or_glm(s) %>% rename_with(~ paste0(.x, "_fuertes"))) %>%
      bind_cols(or_glm(r) %>% rename_with(~ paste0(.x, "_resto")))
  }
}
p6 <- bind_rows(
  feat %>% mutate(conjunto = "Todos"),
  feat %>% filter(seq_id %in% ids_conc) %>% mutate(conjunto = "Actividad concordante entre réplicas")
) %>%
  group_by(conjunto, rep) %>%
  group_modify(function(d, k) {
    map_dfr(hipotesis_vars$variable, ~ compara_variable(d, .x) %>% mutate(variable = .x)) %>%
      mutate(
        n_fuerte_alta = sum(d$tata_fuerte & d$act_alta), n_fuerte_baja = sum(d$tata_fuerte & !d$act_alta),
        p_adj = p.adjust(p, "BH"), p_interaccion_adj = p.adjust(p_interaccion, "BH")
      )
  }) %>%
  ungroup() %>%
  left_join(hipotesis_vars, by = "variable") %>%
  select(conjunto, rep, hipotesis, variable, tipo, n_fuerte_alta, n_fuerte_baja, everything())

# Puntaje de la TATA (solo existe en las TATA fuertes): alta vs. baja
p6_puntaje <- bind_rows(
  feat %>% mutate(conjunto = "Todos"),
  feat %>% filter(seq_id %in% ids_conc) %>% mutate(conjunto = "Actividad concordante entre réplicas")
) %>%
  filter(tata_fuerte) %>%
  group_by(conjunto, rep) %>%
  summarise(
    mediana_alta = median(score_canon[act_alta]), mediana_baja = median(score_canon[!act_alta]),
    p = wilcox.test(score_canon[act_alta], score_canon[!act_alta], exact = FALSE)$p.value, .groups = "drop"
  ) %>%
  mutate(dif_medianas = mediana_alta - mediana_baja)
save_tsv(p6, "p6_variables_tata_fuertes_alta_vs_baja.tsv")
save_tsv(p6_puntaje, "p6_puntaje_tata_fuertes_alta_vs_baja.tsv")
show(p6_puntaje, "p6 puntaje de las TATA fuertes, alta vs. baja actividad")
show(p6 %>% filter(conjunto == "Todos", tipo == "binaria") %>%
  select(rep, hipotesis, variable, fuerte_alta, fuerte_baja, OR_fisher, p_adj, resto_alta, resto_baja, OR_resto, p_interaccion, p_interaccion_adj), "p6 binarias (Todos)")
show(p6 %>% filter(tipo == "continua") %>%
  select(conjunto, rep, variable, fuerte_alta, fuerte_baja, dif_medianas_fuertes, p_adj, dif_medias_fuertes, dif_medias_resto, p_interaccion), "p6 continuas")
show(p6 %>% filter(conjunto != "Todos", tipo == "binaria") %>%
  select(rep, variable, n_fuerte_alta, n_fuerte_baja, fuerte_alta, fuerte_baja, OR_fisher, p_adj, OR_resto, p_interaccion), "p6 binarias (concordantes)")

forest_p6 <- function(df) {
  df %>%
    filter(tipo == "binaria") %>%
    select(rep, hipotesis, variable, p_interaccion_adj, starts_with("OR_") & (ends_with("_fuertes") | ends_with("_resto"))) %>%
    pivot_longer(starts_with("OR_"), names_to = c(".value", "grupo"), names_pattern = "(OR|OR_low|OR_high)_(fuertes|resto)") %>%
    mutate(
      grupo = recode(grupo, fuertes = "TATA fuertes", resto = "Resto"),
      variable = factor(variable, levels = rev(hipotesis_vars$variable)),
      interaccion = if_else(p_interaccion_adj < 0.05, "Interacción BH < 0.05", "Interacción n.s.", missing = "Interacción n.s.")
    ) %>%
    filter(!is.na(OR)) %>%
    ggplot(aes(OR, variable, color = grupo)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_linerange(aes(xmin = OR_low, xmax = OR_high), position = position_dodge(width = 0.6), linewidth = 0.9) +
    geom_point(aes(shape = interaccion), position = position_dodge(width = 0.6), size = 2.6, fill = "white") +
    scale_x_log10() +
    scale_color_manual(values = c("TATA fuertes" = col_misma, "Resto" = "grey55")) +
    scale_shape_manual(values = c("Interacción BH < 0.05" = 16, "Interacción n.s." = 21)) +
    facet_grid(hipotesis ~ rep, scales = "free_y", space = "free_y") +
    th +
    theme(legend.box = "vertical") +
    labs(x = "OR de actividad alta vs. baja (escala log, IC 95%)", y = NULL, color = NULL, shape = NULL)
}
save_jpg(forest_p6(filter(p6, conjunto == "Todos")), "p6_forest_alta_vs_baja", 13, 14)
save_jpg(forest_p6(filter(p6, conjunto != "Todos")), "p6_forest_alta_vs_baja_concordantes", 13, 14)

# --- p7: especificidad tisular como variable continua ----------------------------

# Gini de la expresion en las muestras de FANTOM5 (sample_specificity_gini,
# la variable de la que salen los terciles de "Alta especificidad tisular";
# NA en los no detectados) y numero de muestras de FANTOM5 con el promotor
# activo (active_fantom5_samples, 0 en los no detectados). (a) TATA fuertes
# y resto, actividad alta vs. baja: medianas, HL y Wilcoxon; interaccion
# actividad x TATA fuerte con lm sobre el percentil de la variable (dentro
# de cada replica). (b) Relacion continua: Spearman del percentil de
# actividad con la especificidad en TATA fuertes vs. resto (diferencia por
# z de Fisher) y lm(act_pct ~ gini * tata_fuerte); percentil de actividad
# medio por quintil de Gini (quintiles sobre todos los promotores de cada
# replica).

spec_vars <- c(gini = "Gini (FANTOM5)", n_muestras = "Muestras de FANTOM5 con actividad")
spec <- feat %>%
  select(seq_id, rep, act_alta, tata_fuerte) %>%
  left_join(data %>% select(seq_id, rep, act_pct, gini = sample_specificity_gini, n_muestras = active_fantom5_samples), by = c("seq_id", "rep")) %>%
  mutate(grupo = factor(if_else(tata_fuerte, "TATA fuertes", "Resto"), levels = c("TATA fuertes", "Resto")))

p7_comp <- spec %>%
  pivot_longer(c(gini, n_muestras), names_to = "variable", values_to = "valor") %>%
  filter(!is.na(valor)) %>%
  group_by(rep, variable) %>%
  mutate(valor_pct = 100 * percent_rank(valor)) %>%
  group_modify(function(d, k) {
    p_int <- coef(summary(lm(valor_pct ~ act_alta * tata_fuerte, data = d)))["act_altaTRUE:tata_fuerteTRUE", 4]
    d %>%
      group_by(grupo) %>%
      group_modify(function(g, kk) {
        w <- wilcox.test(g$valor[g$act_alta], g$valor[!g$act_alta], conf.int = TRUE, exact = FALSE)
        tibble(
          n_alta = sum(g$act_alta), n_baja = sum(!g$act_alta),
          mediana_alta = median(g$valor[g$act_alta]), mediana_baja = median(g$valor[!g$act_alta]),
          HL_alta_menos_baja = unname(w$estimate), HL_low = w$conf.int[1], HL_high = w$conf.int[2], p = w$p.value,
          dif_pct_medio = mean(g$valor_pct[g$act_alta]) - mean(g$valor_pct[!g$act_alta])
        )
      }) %>%
      ungroup() %>%
      mutate(p_interaccion = p_int)
  }) %>%
  ungroup() %>%
  mutate(variable = spec_vars[variable])
save_tsv(p7_comp, "p7_especificidad_continua_alta_vs_baja.tsv")
show(p7_comp, "p7 especificidad continua: actividad alta vs. baja, TATA fuertes y resto")

n_na_gini <- spec %>% group_by(rep, grupo, act_alta) %>% summarise(n = n(), n_sin_gini = sum(is.na(gini)), .groups = "drop")
show(n_na_gini, "p7 promotores sin Gini (no detectados en FANTOM5)")

fisher_z <- function(r1, n1, r2, n2) 2 * pnorm(-abs((atanh(r1) - atanh(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))))
p7_cor <- spec %>%
  pivot_longer(c(gini, n_muestras), names_to = "variable", values_to = "valor") %>%
  filter(!is.na(valor)) %>%
  group_by(rep, variable) %>%
  group_modify(function(d, k) {
    s <- d[d$tata_fuerte, ]
    r <- d[!d$tata_fuerte, ]
    rho_s <- cor(s$valor, s$act_pct, method = "spearman")
    rho_r <- cor(r$valor, r$act_pct, method = "spearman")
    m <- coef(summary(lm(act_pct ~ valor * tata_fuerte, data = d)))
    tibble(
      n_fuertes = nrow(s), rho_fuertes = rho_s, p_fuertes = cor.test(s$valor, s$act_pct, method = "spearman", exact = FALSE)$p.value,
      n_resto = nrow(r), rho_resto = rho_r, p_dif_rho = fisher_z(rho_s, nrow(s), rho_r, nrow(r)),
      pendiente_resto = m["valor", 1], pendiente_fuertes = m["valor", 1] + m["valor:tata_fuerteTRUE", 1], p_interaccion_lm = m["valor:tata_fuerteTRUE", 4]
    )
  }) %>%
  ungroup() %>%
  mutate(variable = spec_vars[variable])
save_tsv(p7_cor, "p7_especificidad_continua_correlacion.tsv")
show(p7_cor, "p7 actividad vs. especificidad continua (Spearman, pendiente de lm)")

# Percentil de actividad medio por quintil de Gini
p7_quintil <- spec %>%
  filter(!is.na(gini)) %>%
  group_by(rep) %>%
  mutate(quintil = ntile(gini, 5)) %>%
  group_by(rep, quintil) %>%
  mutate(gini_mediana = median(gini)) %>%
  group_by(rep, grupo, quintil, gini_mediana) %>%
  summarise(n = n(), act_media = mean(act_pct), h = qt(0.975, n() - 1) * sd(act_pct) / sqrt(n()), prop_alta = mean(act_alta), .groups = "drop") %>%
  mutate(conf.low = act_media - h, conf.high = act_media + h)
save_tsv(p7_quintil, "p7_actividad_por_quintil_gini.tsv")
show(p7_quintil %>% select(rep, grupo, quintil, gini_mediana, n, act_media, conf.low, conf.high, prop_alta), "p7 percentil de actividad por quintil de Gini")

col_grupo <- scale_color_manual(values = c("TATA fuertes" = col_misma, "Resto" = "grey55"))
panel_quintil <- ggplot(p7_quintil, aes(gini_mediana, act_media, color = grupo)) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.02)) +
  geom_line(position = position_dodge(width = 0.02)) +
  geom_point(size = 3, position = position_dodge(width = 0.02)) +
  geom_text(data = filter(p7_quintil, grupo == "TATA fuertes"), aes(y = conf.low, label = paste0("n=", n)), vjust = 1.6, size = 4, show.legend = FALSE) +
  col_grupo +
  facet_wrap(~rep) +
  th +
  labs(x = "Gini de la expresión en FANTOM5 (mediana del quintil)", y = "Percentil de actividad\n(media, IC 95%)", color = NULL)
panel_violin <- spec %>%
  filter(!is.na(gini)) %>%
  mutate(actividad = factor(if_else(act_alta, "Alta", "Baja"), levels = c("Baja", "Alta"))) %>%
  ggplot(aes(actividad, gini, fill = grupo)) +
  geom_violin(position = position_dodge(width = 0.8), alpha = 0.5, color = NA) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.15, outlier.shape = NA, show.legend = FALSE) +
  scale_fill_manual(values = c("TATA fuertes" = col_misma, "Resto" = "grey60")) +
  facet_wrap(~rep) +
  th +
  labs(x = "Actividad (percentil 75 por réplica)", y = "Gini de la expresión\nen FANTOM5", fill = NULL)
save_jpg(panel_violin / panel_quintil, "p7_especificidad_continua", 12, 11)

# --- p8: ruido de las TATA fuertes segun la especificidad tisular ---------------

# AUC del ruido de las TATA fuertes vs. sin motivo (como p4) dentro de cada
# tercil de Gini (terciles sobre todos los promotores con Gini de cada
# replica, los mismos de "Alta/Baja especificidad tisular") y de cada
# estrato de actividad; la referencia es "sin motivo" del mismo tercil y
# estrato. Los no detectados en FANTOM5 van como categoria aparte. Ademas,
# relacion continua dentro de las TATA fuertes: Spearman del rango de
# varianza con el Gini por estrato de actividad, y var_rank_sw ~ gini *
# tata_fuerte (TATA fuertes + sin motivo) por estrato.

ESPEC_LAB <- c("1" = "Baja", "2" = "Media", "3" = "Alta", "ND" = "No detectado")
ruido_spec <- data %>%
  select(seq_id, rep, var_rank_sw, actividad, clase, score_canon, gini = sample_specificity_gini) %>%
  mutate(tata_fuerte = clase == "Canónica" & score_canon >= CORTE_MEDIANA - TATA_TOL) %>%
  group_by(rep) %>%
  mutate(espec = if_else(is.na(gini), "ND", as.character(ntile(gini, 3)))) %>%
  ungroup() %>%
  mutate(espec = factor(ESPEC_LAB[espec], levels = ESPEC_LAB), actividad = factor(actividad, levels = c("Baja", "Alta")))

p8_auc <- bind_rows(
  ruido_spec %>% mutate(actividad = "Todas"),
  ruido_spec %>% mutate(actividad = as.character(actividad))
) %>%
  mutate(actividad = factor(actividad, levels = c("Todas", "Baja", "Alta"))) %>%
  group_by(rep, actividad, espec) %>%
  group_modify(~ auc_vs(.x$var_rank_sw[.x$tata_fuerte], .x$var_rank_sw[.x$clase == "Sin motivo"])) %>%
  ungroup()
save_tsv(p8_auc, "p8_auc_tata_fuertes_por_especificidad.tsv")
show(p8_auc %>% select(rep, actividad, espec, n, n_ref, AUC, conf.low, conf.high, p), "p8 AUC del ruido de las TATA fuertes vs. sin motivo, por tercil de Gini")

# Diferencia de AUC entre terciles extremos (z de DeLong, muestras independientes)
p8_dif <- p8_auc %>%
  filter(espec %in% c("Baja", "Alta")) %>%
  select(rep, actividad, espec, AUC, se) %>%
  pivot_wider(names_from = espec, values_from = c(AUC, se)) %>%
  mutate(dif_alta_menos_baja = AUC_Alta - AUC_Baja, p = auc_diff_p(AUC_Alta, se_Alta, AUC_Baja, se_Baja))
save_tsv(p8_dif, "p8_dif_auc_especificidad_alta_vs_baja.tsv")
show(p8_dif, "p8 AUC especificidad alta - baja")

p8_cont <- bind_rows(
  ruido_spec %>% mutate(actividad = "Todas"),
  ruido_spec %>% mutate(actividad = as.character(actividad))
) %>%
  filter(!is.na(gini), tata_fuerte | clase == "Sin motivo") %>%
  group_by(rep, actividad) %>%
  group_modify(function(d, k) {
    s <- d[d$tata_fuerte, ]
    m <- coef(summary(lm(var_rank_sw ~ gini * tata_fuerte, data = d)))
    ct <- cor.test(s$gini, s$var_rank_sw, method = "spearman", exact = FALSE)
    tibble(
      n_fuertes = nrow(s), rho_fuertes = unname(ct$estimate), p_rho = ct$p.value,
      pendiente_sin_motivo = m["gini", 1], pendiente_fuertes = m["gini", 1] + m["gini:tata_fuerteTRUE", 1],
      p_interaccion = m["gini:tata_fuerteTRUE", 4]
    )
  }) %>%
  ungroup()
save_tsv(p8_cont, "p8_ruido_vs_gini_continuo.tsv")
show(p8_cont, "p8 ruido vs. Gini (continuo)")

panel_p8 <- p8_auc %>%
  filter(actividad != "Todas") %>%
  ggplot(aes(espec, AUC, color = actividad)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = position_dodge(width = 0.5), linewidth = 1) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_text(aes(y = conf.high, label = paste0("n=", n)), vjust = -0.6, size = 4, position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_color_manual(values = c(Baja = col_opuesta, Alta = col_misma)) +
  facet_wrap(~rep) +
  th +
  labs(x = "Especificidad tisular (tercil de Gini en FANTOM5)", y = "AUC del ruido de las TATA fuertes\nvs. sin motivo (IC 95%)", color = "Actividad")
save_jpg(panel_p8, "p8_ruido_tata_fuertes_por_especificidad", 12, 6.5)

# --- p9: perfil posicional y controles de composicion por estrato de actividad ---

# p1 y p2 repetidos dentro de cada estrato de actividad (percentil 75 de mean
# por replica): todas las comparaciones (casos y referencia "sin motivo") se
# hacen dentro del estrato. Los graficos muestran el estrato alto; el bajo
# queda en las tablas como contraste.

por_estrato <- function(f) {
  map_dfr(c("Alta", "Baja"), ~ f(filter(data, actividad == .x)) %>% mutate(actividad = .x))
}

perfil_act <- por_estrato(function(d) profile_auc(d, hits))
save_tsv(perfil_act, "p9_auc_por_ventana_por_actividad.tsv")
show(perfil_act %>% filter(between(win, -63, 6)) %>% select(actividad, rep, hebra, win, n, AUC, conf.low, conf.high, p_adj), "p9 AUC por ventana (-63 a +6), por estrato de actividad")
show(perfil_act %>% group_by(actividad, rep, hebra) %>% summarise(n_ventanas = n(), n_BH = sum(sig), mediana_AUC_fuera_canonica = median(AUC[win != TATA_WIN_ORIGIN]), .groups = "drop"), "p9 resumen del perfil por estrato")

clases_act <- por_estrato(function(dd) {
  dd %>%
    group_by(rep) %>%
    group_modify(function(d, k) {
      ref <- sin_motivo(d)
      bind_rows(
        map_dfr(setdiff(TATA_CLASES, "Sin motivo"), ~ auc_vs(d$var_rank_sw[d$clase == .x], ref) %>% mutate(clase = .x)),
        auc_vs(d$var_rank_sw[d$seq_id %in% op_canon_solo], ref) %>% mutate(clase = "Hebra opuesta en −31 a −26, sin canónica")
      ) %>% mutate(p_adj = p.adjust(p, "BH"))
    }) %>%
    ungroup()
}) %>%
  mutate(clase = fct_rev(fct_inorder(clase)))
save_tsv(clases_act, "p9_auc_por_clase_por_actividad.tsv")
show(clases_act %>% select(actividad, rep, clase, n, n_ref, AUC, conf.low, conf.high, p_adj), "p9 AUC por clase vs. sin motivo del estrato")

perfil_perm_act <- imap_dfr(perms, function(p, nm) {
  s <- tata_scores(code, TATA_W[, p])
  h <- tata_hits(s$misma, s$opuesta, CORTE)
  por_estrato(function(d) profile_auc(d, filter(h, hebra == "Misma hebra"), ref_hits = h)) %>%
    mutate(perm = nm)
})
perfil_real_act <- filter(perfil_act, hebra == "Misma hebra")
perm_canon_act <- bind_rows(perfil_real_act %>% mutate(perm = "Matriz real"), perfil_perm_act) %>%
  filter(win == TATA_WIN_ORIGIN) %>%
  select(actividad, rep, perm, n, AUC, conf.low, conf.high, p_adj)
save_tsv(bind_rows(perfil_real_act %>% mutate(perm = "Matriz real"), perfil_perm_act), "p9_auc_por_ventana_permutadas_por_actividad.tsv")
save_tsv(perm_canon_act, "p9_permutadas_ventana_canonica_por_actividad.tsv")
show(perm_canon_act, "p9 ventana -33 a -24 (misma hebra), matriz real vs. permutadas, por estrato")

modelos_act <- data %>%
  mutate(at_prox = at_prox[seq_id], at_frag = 1 - g_c, clase = relevel(clase, ref = "Sin motivo")) %>%
  group_by(actividad, rep) %>%
  group_modify(function(d, k) {
    bind_rows(
      broom::tidy(lm(var_rank_sw ~ clase, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[1]),
      broom::tidy(lm(var_rank_sw ~ clase + at_prox, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[2]),
      broom::tidy(lm(var_rank_sw ~ clase + at_prox + at_frag, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[3])
    )
  }) %>%
  ungroup() %>%
  filter(str_detect(term, "^clase|^at_")) %>%
  mutate(modelo = factor(modelo, levels = modelos_lab))
save_tsv(modelos_act, "p9_modelo_clase_AT_por_actividad.tsv")
show(modelos_act %>% select(actividad, rep, modelo, term, estimate, conf.low, conf.high, p.value), "p9 modelo lineal var_rank_sw ~ clase + A/T, por estrato")

# Graficos (estrato alto), mismos formatos que p1 y p2
perfil_alta <- filter(perfil_act, actividad == "Alta")
panel_perfil_alta <- ggplot(perfil_alta, aes(win_center, AUC, color = hebra)) +
  franja +
  geom_hline(yintercept = 0.5, color = "grey50") +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = dodge, alpha = 0.7) +
  geom_line(aes(group = seg), position = dodge, alpha = 0.6) +
  geom_point(aes(shape = sig), position = dodge, size = 2.2, fill = "white") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "BH < 0.05", `FALSE` = "n.s."), name = NULL) +
  col_hebra +
  facet_grid(rep ~ .) +
  scale_x_continuous(breaks = seq(-240, 0, 20), labels = minus) +
  th +
  labs(x = NULL, y = "AUC del ruido vs. sin motivo\n(actividad alta, IC 95%)", color = NULL)
panel_prop_alta <- ggplot(perfil_alta, aes(win_center, 100 * prop, color = hebra)) +
  franja +
  geom_line(aes(group = seg)) +
  geom_point(size = 1.3) +
  col_hebra +
  facet_grid(. ~ rep) +
  scale_x_continuous(breaks = seq(-200, 0, 50), labels = minus) +
  th +
  theme(legend.position = "none") +
  labs(x = "Posición de la coincidencia respecto del TSS (ventanas de 10 pb)", y = "% de promotores\ncon coincidencia")
save_jpg(panel_perfil_alta / panel_prop_alta + plot_layout(heights = c(2.2, 1)), "p9_perfil_posicion_auc_actividad_alta", 14, 10)
save_jpg(forest(filter(clases_act, actividad == "Alta"), "clase", xlab = "AUC del ruido vs. sin motivo, actividad alta (IC 95%)"), "p9_clases_auc_actividad_alta", 10, 5)

panel_perm_alta <- ggplot(filter(perfil_perm_act, actividad == "Alta"), aes(win_center, AUC)) +
  franja +
  geom_hline(yintercept = 0.5, color = "grey50") +
  geom_line(aes(group = paste(perm, seg)), color = "grey65", alpha = 0.7) +
  geom_line(data = filter(perfil_real_act, actividad == "Alta"), aes(group = seg), color = col_misma, linewidth = 1.1) +
  geom_point(data = filter(perfil_real_act, actividad == "Alta"), color = col_misma, size = 2) +
  facet_wrap(~rep, ncol = 1) +
  scale_x_continuous(breaks = seq(-200, 0, 50), labels = minus) +
  th +
  labs(x = "Posición respecto del TSS (ventanas de 10 pb)", y = "AUC del ruido vs. sin motivo\n(actividad alta)")
save_jpg(panel_perm_alta, "p9_permutadas_auc_actividad_alta", 10, 8)

panel_modelo_act <- modelos_act %>%
  filter(str_detect(term, "^clase")) %>%
  mutate(
    clase = factor(str_remove(term, "^clase"), levels = rev(setdiff(TATA_CLASES, "Sin motivo"))),
    actividad = factor(paste("Actividad", tolower(actividad)), levels = c("Actividad alta", "Actividad baja"))
  ) %>%
  ggplot(aes(estimate, clase, color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.6), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  scale_color_manual(values = setNames(c("grey55", col_opuesta, col_misma), modelos_lab)) +
  facet_grid(rep ~ actividad) +
  th +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Diferencia de rango de varianza vs. sin motivo del estrato", y = NULL, color = NULL)
save_jpg(panel_modelo_act, "p9_modelo_AT_por_actividad", 14, 8)

# --- p10: clase x actividad con ajuste por A/T en un solo modelo ---------------

# var_rank_sw ~ clase * act_alta (+ A/T proximal y del fragmento), por replica,
# sobre todos los promotores. Efecto de cada clase vs. sin motivo en baja
# actividad (coeficiente de clase) y en alta (clase + interaccion, IC por
# combinacion lineal con vcov); p de la interaccion = el efecto difiere entre
# estratos. Modelos: sin A/T; con A/T de pendiente comun (el pedido); con A/T
# x actividad (sensibilidad: pendientes de A/T distintas por estrato).

modelos_int_lab <- c("Sin A/T", "+ A/T (pendiente común)", "+ A/T × actividad")
d10 <- data %>%
  mutate(at_prox = at_prox[seq_id], at_frag = 1 - g_c, clase = relevel(clase, ref = "Sin motivo"), act_alta = actividad == "Alta")
formulas_int <- list(
  var_rank_sw ~ clase * act_alta,
  var_rank_sw ~ clase * act_alta + at_prox + at_frag,
  var_rank_sw ~ clase * act_alta + (at_prox + at_frag) * act_alta
) %>% set_names(modelos_int_lab)
p10 <- d10 %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    imap_dfr(formulas_int, function(f, nm) {
      m <- lm(f, data = d)
      b <- coef(m)
      V <- vcov(m)
      map_dfr(setdiff(TATA_CLASES, "Sin motivo"), function(cl) {
        t_cl <- paste0("clase", cl)
        t_int <- paste0(t_cl, ":act_altaTRUE")
        se_alta <- sqrt(V[t_cl, t_cl] + V[t_int, t_int] + 2 * V[t_cl, t_int])
        df_res <- m$df.residual
        q <- qt(0.975, df_res)
        tibble(
          modelo = nm, clase = cl,
          efecto_baja = b[[t_cl]], baja_low = b[[t_cl]] - q * sqrt(V[t_cl, t_cl]), baja_high = b[[t_cl]] + q * sqrt(V[t_cl, t_cl]),
          p_baja = 2 * pt(-abs(b[[t_cl]] / sqrt(V[t_cl, t_cl])), df_res),
          efecto_alta = b[[t_cl]] + b[[t_int]], alta_low = efecto_alta - q * se_alta, alta_high = efecto_alta + q * se_alta,
          p_alta = 2 * pt(-abs(efecto_alta / se_alta), df_res),
          interaccion = b[[t_int]], p_interaccion = coef(summary(m))[t_int, 4]
        )
      })
    })
  }) %>%
  ungroup() %>%
  mutate(modelo = factor(modelo, levels = modelos_int_lab), clase = factor(clase, levels = rev(setdiff(TATA_CLASES, "Sin motivo"))))
save_tsv(p10, "p10_modelo_clase_x_actividad_AT.tsv")
show(p10 %>% select(rep, modelo, clase, efecto_baja, baja_low, baja_high, p_baja, efecto_alta, alta_low, alta_high, p_alta, interaccion, p_interaccion), "p10 var_rank_sw ~ clase * actividad (+ A/T)")

# Fraccion del efecto de la canonica que queda tras el ajuste, por estrato
p10_retenido <- p10 %>%
  filter(clase == "Canónica") %>%
  select(rep, modelo, efecto_baja, efecto_alta) %>%
  pivot_longer(c(efecto_baja, efecto_alta), names_to = "actividad", values_to = "efecto") %>%
  group_by(rep, actividad) %>%
  mutate(fraccion_del_sin_AT = efecto / efecto[modelo == "Sin A/T"]) %>%
  ungroup()
save_tsv(p10_retenido, "p10_fraccion_retenida_canonica.tsv")
show(p10_retenido, "p10 efecto de la canonica retenido tras el ajuste por A/T")

panel_p10 <- bind_rows(
  p10 %>% transmute(rep, modelo, clase, actividad = "Actividad alta", efecto = efecto_alta, low = alta_low, high = alta_high, p_interaccion),
  p10 %>% transmute(rep, modelo, clase, actividad = "Actividad baja", efecto = efecto_baja, low = baja_low, high = baja_high, p_interaccion)
) %>%
  mutate(actividad = factor(actividad, levels = c("Actividad alta", "Actividad baja"))) %>%
  ggplot(aes(efecto, clase, color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = low, xmax = high), position = position_dodge(width = 0.6), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  scale_color_manual(values = setNames(c("grey55", col_misma, col_opuesta), modelos_int_lab)) +
  facet_grid(rep ~ actividad) +
  th +
  guides(color = guide_legend(nrow = 1)) +
  labs(x = "Diferencia de rango de varianza vs. sin motivo (modelo clase × actividad)", y = NULL, color = NULL)
save_jpg(panel_p10, "p10_modelo_clase_x_actividad_AT", 14, 8)

message("Salidas en ", out_dir)
