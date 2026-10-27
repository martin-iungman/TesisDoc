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

out_dir <- "figures/exploratorio_tata_ruido"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

CORTE <- TATA_CORTE_EPD
CANON <- TATA_CANON
CONTROL <- TATA_CANON - 100
UMBRAL_ALTO <- -5 # puntaje "alto" de p4 si da n razonables (ver abajo)
LIMITES <- c(-Inf, -16, -14, -12, -10, CORTE, -7, -6, -5, Inf)
ETIQUETAS <- c("< −16", "−16 a −14", "−14 a −12", "−12 a −10", "−10 a −8,17", "−8,17 a −7", "−7 a −6", "−6 a −5", "≥ −5")
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
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "BH < 0,05", `FALSE` = "n.s."), name = NULL) +
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
  annotate("text", x = CORTE, y = Inf, label = "−8,17 (EPD)", hjust = 1.05, vjust = 1.4, size = 5, color = "grey30") +
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

message("Salidas en ", out_dir)
