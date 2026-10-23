# EXPLORATORIO - motivos de factores de transcripcion adicionales.
#
# Parte A (motivos dirigidos): matrices de JASPAR 2024 CORE vertebrados (no
# redundante) para NRF1, ETS (GABPA, ELK1), YY1, CREB1/ATF1, ZNF143, E-box
# (USF1, MAX, MYC) y, como controles positivos, NFYA (CCAAT) y SP1 (GC-box).
# Busqueda en ambas hebras de todo el fragmento con puntaje relativo >= 0,85
# ((puntaje - minimo) / (maximo - minimo) de la PWM log2, como el relScore de
# TFBSTools); presencia = al menos una coincidencia. Por motivo: prevalencia
# con / sin isla CpG (OR de Fisher), efecto sobre la actividad, posicion de
# las coincidencias respecto del TSS (con un control de composicion: las
# mismas busquedas con matrices de columnas permutadas) y efecto dentro vs.
# fuera de la ventana de mayor densidad. NFYA y SP1 se comparan con
# CCAAT_EPD y GCbox_EPD para calibrar el corte.
#
# Parte B (sin hipotesis): los 2080 6-meros canonicos (cada uno colapsado
# con su complemento reverso), presencia en todo el fragmento; efecto de
# cada uno en todos / con / sin isla CpG, relacion del efecto con el
# contenido de G+C y de CpG, k-meros destacados (consistentes en ambas
# replicas) y su parecido con matrices de JASPAR (universalmotif).
#
# Efectos: Wilcoxon + Hodges-Lehmann (IC 95%) de la actividad media con vs.
# sin, por replica; BH dentro de replica (parte A) o de replica x grupo
# (parte B). Isla CpG = CGI_frag. Fragmento orientado de 252 pb, -236 a +15,
# TSS de EPD = 0 en la posicion 237 (TATA_TSS_POS en R/functions/tata_scan.R).
#
# Figura de la tesis motivos_de_novo (numero vigente en docs/mapping_figuras.csv),
# escrita via fig_dir(): kmeros_efecto_vs_cpg, kmeros_volcano_motivos y
# ebox_zeb_scatter_cgi (.jpg y .pdf), en ese orden, y la figura completa
# motivos_de_novo.jpg/.pdf. kmeros_destacados es exploratorio. Todo lo demas (graficos y tablas)
# es exploratorio y va a figures/exploratorio_motivos_tf/. Los calculos
# costosos se cachean en data/processed/.
# Requiere: data/processed/activity_stats_highconf.tsv, prom_df.tsv,
# ../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa y
# data/external/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
# (ver data/external/JASPAR/README.md). Paquetes: universalmotif (TFBSTools no
# se usa: no carga en R 4.4 porque su dependencia TFMPvalue salio de CRAN).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
library(Biostrings)
library(universalmotif)
source("R/functions/tata_scan.R") # TATA_TSS_POS, TATA_FRAG_LEN, tata_hl_vs()
source("R/functions/plot_helpers.R") # thesis_clr
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename # Biostrings/S4Vectors la enmascara

out_dir <- "figures/exploratorio_motivos_tf"
dir_novo <- fig_dir("motivos_de_novo") # figura de la tesis (busqueda de novo con 6-meros)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
path_fa <- "../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa"
path_jaspar <- "data/external/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
cache_hits <- "data/processed/motivos_tf_jaspar_hits_v1.rds"
# v2: recalculado con el filtro de >1000 cuentas normalizadas por replica
cache_kmeros <- "data/processed/motivos_tf_kmeros_hl_v2.rds"

TSS_POS <- TATA_TSS_POS
REL_SCORE <- 0.85
MIN_N <- 10
N_PERM <- 10 # matrices de columnas permutadas por motivo (control de posicion)
VENTANA_ANCHO <- 30 # ancho de la ventana de mayor densidad de coincidencias
ENRIQ_MIN <- 1.5 # ventana "preferida": >= 1,5 veces la fraccion esperada por las permutadas
K <- 6
PREV_RANGO <- c(0.01, 0.99)
N_TOP <- 20
GRUPOS <- c("Todos", "Con isla CpG", "Sin isla CpG")

# Motivos de la parte A: nombre en JASPAR (altname) -> etiqueta y familia
MOTIVOS <- tribble(
  ~jaspar, ~motivo, ~familia,
  "Nrf1", "NRF1", "NRF1",
  "GABPA", "GABPA", "ETS",
  "ELK1", "ELK1", "ETS",
  "Yy1", "YY1", "YY1",
  "CREB1", "CREB1", "CRE",
  "Atf1", "ATF1", "CRE",
  "ZNF143", "ZNF143", "ZNF143",
  "USF1", "USF1", "E-box",
  "MAX", "MAX", "E-box",
  "MYC", "MYC", "E-box",
  "NFYA", "NFYA", "Control: CCAAT",
  "SP1", "SP1", "Control: GC-box"
)

col_grupo <- c("Todos" = "grey30", "Con isla CpG" = "#216869", "Sin isla CpG" = "#14AFB2")
col_hebra <- c("Misma hebra" = "#216869", "Hebra opuesta" = "#14AFB2")
col_cpg <- c("0" = "grey70", "1" = "#7FB800", "2" = "#D6741F", "3" = "#B5345A")
th <- theme_pubclean(base_size = 13) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")
minus <- function(x) sub("^-", "−", format(x, trim = TRUE))
save_tsv <- function(x, name) write_tsv(x, file.path(out_dir, name))
# paneles de la figura de la tesis: .jpg y .pdf
save_fig <- function(p, dir, name, width, height) {
  for (ext in c("jpg", "pdf")) ggsave(file.path(dir, paste0(name, ".", ext)), p, width = width, height = height, units = "in")
}
sig_lab <- function(p_adj) factor(if_else(!is.na(p_adj) & p_adj < 0.05, "BH < 0,05", "n.s."), c("BH < 0,05", "n.s."))

forest <- function(res, y, color, valores_color, x_lab, facets = . ~ rep, dodge = 0.7) {
  # group explicito: sin el, geom_point se agrupa por color x forma (sig) y
  # el dodge de los puntos no coincide con el de las barras
  ggplot(res, aes(estimate, {{ y }}, color = {{ color }}, group = {{ color }})) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = dodge), linewidth = 0.9) +
    geom_point(aes(shape = sig), size = 2.6, fill = "white", position = position_dodge(width = dodge)) +
    scale_shape_manual(values = c(`BH < 0,05` = 16, `n.s.` = 21), name = NULL) +
    scale_color_manual(values = valores_color, name = NULL) +
    facet_grid(facets) +
    th +
    labs(x = x_lab, y = NULL)
}

# --- Datos -------------------------------------------------------------------

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fa <- readDNAStringSet(path_fa)
fa <- fa[prom_df$seq_id]
stopifnot(all(width(fa) == TATA_FRAG_LEN), !anyNA(names(fa)))

data <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE) %>%
  inner_join(prom_df, by = c("seq_id", "name"))
data_g <- bind_rows(
  data %>% mutate(grupo = "Todos"),
  data %>% filter(!is.na(CGI_frag)) %>% mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG"))
) %>% mutate(grupo = factor(grupo, levels = GRUPOS))

# =============================================================================
# PARTE A - motivos dirigidos (JASPAR)
# =============================================================================

jaspar <- read_jaspar(path_jaspar)
jaspar_alt <- map_chr(jaspar, ~ .x@altname)
motivos_pcm <- map(MOTIVOS$jaspar, function(nm) {
  i <- which(jaspar_alt == nm)
  if (length(i) != 1) stop("JASPAR: se esperaba una sola matriz para ", nm, " y hay ", length(i))
  jaspar[[i]]
}) %>% set_names(MOTIVOS$motivo)
matrices <- MOTIVOS %>%
  mutate(
    jaspar_id = map_chr(motivos_pcm, ~ .x@name),
    largo = map_int(motivos_pcm, ~ ncol(.x@motif)),
    consenso = map_chr(motivos_pcm, ~ .x@consensus),
    n_sitios = map_dbl(motivos_pcm, ~ .x@nsites)
  )
print(matrices)

# Coincidencias con puntaje relativo >= REL_SCORE en ambas hebras. La
# posicion es la del centro de la coincidencia respecto del TSS; "Misma
# hebra" = hebra del promotor (las secuencias ya estan orientadas).
buscar_motivo <- function(pcm, nombre) {
  pwm <- convert_type(pcm, "PWM")@motif
  s_min <- sum(apply(pwm, 2, min))
  s_max <- sum(apply(pwm, 2, max))
  corte <- s_min + REL_SCORE * (s_max - s_min)
  res <- scan_sequences(pcm, fa,
    threshold = corte, threshold.type = "logodds.abs", RC = TRUE,
    calc.pvals = FALSE, warn.NA = FALSE, verbose = 0
  )
  as_tibble(as.data.frame(res)) %>%
    transmute(
      motivo = nombre, seq_id = as.character(sequence),
      inicio = pmin(start, stop), fin = pmax(start, stop),
      hebra = if_else(strand == "+", "Misma hebra", "Hebra opuesta"),
      score, rel_score = (score - s_min) / (s_max - s_min),
      posicion = (inicio + fin) / 2 - TSS_POS
    )
}

# matrices con las columnas permutadas (misma composicion, sin estructura)
permutar_motivo <- function(pcm, semilla) {
  set.seed(semilla)
  m <- unname(pcm@motif[, sample(ncol(pcm@motif)), drop = FALSE])
  rownames(m) <- rownames(pcm@motif)
  create_motif(m, type = "PCM", name = paste0(pcm@name, "_perm", semilla), pseudocount = pcm@pseudocount, nsites = pcm@nsites)
}

if (file.exists(cache_hits)) {
  cache <- readRDS(cache_hits)
  hits <- cache$hits
  hits_perm <- cache$hits_perm
} else {
  t0 <- Sys.time()
  hits <- imap_dfr(motivos_pcm, buscar_motivo)
  hits_perm <- imap_dfr(motivos_pcm, function(pcm, nombre) {
    map_dfr(seq_len(N_PERM), ~ buscar_motivo(permutar_motivo(pcm, 1000 + .x), nombre) %>% mutate(perm = .x))
  })
  saveRDS(list(hits = hits, hits_perm = hits_perm), cache_hits)
  message("Busqueda JASPAR: ", format(Sys.time() - t0))
}
stopifnot(all(hits$rel_score >= REL_SCORE - 1e-9))

# presencia por promotor (al menos una coincidencia en cualquier hebra)
presencia <- hits %>%
  distinct(seq_id, motivo) %>%
  mutate(presente = TRUE) %>%
  pivot_wider(names_from = motivo, values_from = presente)
prom_mot <- prom_df %>%
  select(seq_id, CGI_frag, CCAAT_EPD, GCbox_EPD) %>%
  left_join(presencia, by = "seq_id") %>%
  mutate(across(all_of(MOTIVOS$motivo), ~ replace_na(.x, FALSE)))
data_g <- data_g %>% left_join(prom_mot %>% select(seq_id, all_of(MOTIVOS$motivo)), by = "seq_id")

# --- prevalencia con / sin isla CpG (library completa) ----------------------

prevalencia <- map_dfr(MOTIVOS$motivo, function(m) {
  d <- prom_mot %>% filter(!is.na(CGI_frag))
  x <- d[[m]]
  ft <- fisher.test(table(factor(x, c(TRUE, FALSE)), factor(d$CGI_frag, c(TRUE, FALSE))))
  tibble(
    motivo = m, n = nrow(d), n_con_motivo = sum(x), pct_total = 100 * mean(x),
    pct_con_cgi = 100 * mean(x[d$CGI_frag]), pct_sin_cgi = 100 * mean(x[!d$CGI_frag]),
    odds_ratio = unname(ft$estimate), or_conf.low = ft$conf.int[1], or_conf.high = ft$conf.int[2], p_fisher = ft$p.value
  )
}) %>% mutate(p_fisher_adj = p.adjust(p_fisher, "BH"))

# --- efecto en todos / con / sin isla CpG; BH dentro de replica ------------

efecto_feats <- function(df, feats) {
  df %>%
    group_by(rep) %>%
    group_modify(function(d, k) {
      map_dfr(GRUPOS, function(g) {
        dg <- d[d$grupo == g, ]
        map_dfr(feats, function(f) {
          ind <- dg[[f]]
          ok <- !is.na(ind)
          tata_hl_vs(dg$mean[ok & ind], dg$mean[ok & !ind], MIN_N) %>%
            mutate(grupo = g, motivo = f, n = sum(ok & ind), n_ref = sum(ok & !ind), prop = mean(ind[ok]))
        })
      }) %>% mutate(p_adj = p.adjust(p, "BH"))
    }) %>%
    ungroup() %>%
    mutate(grupo = factor(grupo, GRUPOS)) %>%
    relocate(rep, grupo, motivo, n, n_ref, prop)
}
efecto <- efecto_feats(data_g, MOTIVOS$motivo)

# --- calibracion con EPD (NFYA vs. CCAAT_EPD, SP1 vs. GCbox_EPD) -------------
# La ventana "EPD" es el rango 5-95% de la posicion de las coincidencias de
# JASPAR en los promotores anotados por EPD (EPD busca en ventanas acotadas).

calib_pares <- tibble(motivo = c("NFYA", "SP1"), epd = c("CCAAT_EPD", "GCbox_EPD"))
ventanas_epd <- calib_pares %>%
  pmap_dfr(function(motivo, epd) {
    ids_epd <- prom_mot$seq_id[prom_mot[[epd]]]
    pos <- hits$posicion[hits$motivo == motivo & hits$seq_id %in% ids_epd]
    tibble(motivo = motivo, epd = epd, desde = floor(quantile(pos, 0.05)), hasta = ceiling(quantile(pos, 0.95)))
  })
solapamiento <- function(a, b) {
  tibble(
    n_ambos = sum(a & b), n_solo_jaspar = sum(a & !b), n_solo_epd = sum(!a & b), n_ninguno = sum(!a & !b),
    jaccard = sum(a & b) / sum(a | b), pct_epd_detectado = 100 * sum(a & b) / sum(b),
    pct_jaspar_en_epd = 100 * sum(a & b) / sum(a)
  )
}
calibracion <- ventanas_epd %>%
  pmap_dfr(function(motivo, epd, desde, hasta) {
    en_ventana <- hits %>%
      filter(motivo == !!motivo, between(posicion, desde, hasta)) %>%
      distinct(seq_id) %>%
      pull(seq_id)
    b <- prom_mot[[epd]]
    bind_rows(
      solapamiento(prom_mot[[motivo]], b) %>% mutate(jaspar = "Todo el fragmento"),
      solapamiento(prom_mot$seq_id %in% en_ventana, b) %>% mutate(jaspar = paste0("Ventana EPD (", minus(desde), " a ", minus(hasta), ")"))
    ) %>% mutate(motivo = motivo, epd = epd, .before = 1)
  })
# efecto de EPD y de JASPAR dentro de la ventana EPD, junto al de todo el fragmento
calib_cols <- ventanas_epd %>%
  pmap(function(motivo, epd, desde, hasta) {
    ids <- hits %>% filter(motivo == !!motivo, between(posicion, desde, hasta)) %>% distinct(seq_id) %>% pull(seq_id)
    set_names(list(ids), paste0(motivo, " (ventana EPD)"))
  }) %>%
  flatten()
data_calib <- data_g
for (nm in names(calib_cols)) data_calib[[nm]] <- data_calib$seq_id %in% calib_cols[[nm]]
efecto_calib <- efecto_feats(data_calib, c("NFYA", "CCAAT_EPD", "NFYA (ventana EPD)", "SP1", "GCbox_EPD", "SP1 (ventana EPD)"))

# --- posicion de las coincidencias ------------------------------------------
# Ventana de VENTANA_ANCHO pb con la mayor fraccion de coincidencias (centro
# de la coincidencia), comparada con la fraccion de las coincidencias de las
# matrices permutadas en esa misma ventana (control de composicion
# posicional) y con la esperada si fueran uniformes.

pos_min <- -TSS_POS + 1
pos_max <- TATA_FRAG_LEN - TSS_POS
ventana_preferida <- function(pos_real, pos_perm) {
  inicios <- seq(floor(min(pos_real)), ceiling(max(pos_real)) - VENTANA_ANCHO + 1)
  frac_real <- map_dbl(inicios, ~ mean(pos_real >= .x & pos_real < .x + VENTANA_ANCHO))
  i <- which.max(frac_real)
  desde <- inicios[i]
  hasta <- desde + VENTANA_ANCHO - 1
  frac_perm <- mean(pos_perm >= desde & pos_perm <= hasta)
  rango <- max(pos_real) - min(pos_real) + 1
  tibble(
    desde = desde, hasta = hasta, n_coinc_ventana = sum(pos_real >= desde & pos_real <= hasta), n_coinc = length(pos_real),
    frac_ventana = frac_real[i], frac_uniforme = VENTANA_ANCHO / rango, frac_permutadas = frac_perm,
    enriq_vs_uniforme = frac_real[i] / (VENTANA_ANCHO / rango), enriq_vs_permutadas = frac_real[i] / frac_perm
  )
}
posicion <- MOTIVOS$motivo %>%
  map_dfr(function(m) {
    ventana_preferida(hits$posicion[hits$motivo == m], hits_perm$posicion[hits_perm$motivo == m]) %>% mutate(motivo = m, .before = 1)
  }) %>%
  mutate(preferida = enriq_vs_permutadas >= ENRIQ_MIN & enriq_vs_uniforme >= ENRIQ_MIN & n_coinc_ventana >= 50)

# efecto dentro vs. fuera de la ventana: "Dentro" = alguna coincidencia en la
# ventana; "Solo fuera" = coincidencias pero ninguna en la ventana
efecto_ventana <- posicion %>%
  select(motivo, desde, hasta, preferida) %>%
  pmap_dfr(function(motivo, desde, hasta, preferida) {
    dentro <- hits %>% filter(motivo == !!motivo, between(posicion, desde, hasta)) %>% distinct(seq_id) %>% pull(seq_id)
    con_motivo <- hits %>% filter(motivo == !!motivo) %>% distinct(seq_id) %>% pull(seq_id)
    data_g %>%
      mutate(clase = case_when(seq_id %in% dentro ~ "Dentro", seq_id %in% con_motivo ~ "Solo fuera", TRUE ~ "Sin motivo")) %>%
      group_by(rep, grupo) %>%
      group_modify(function(d, k) {
        bind_rows(
          tata_hl_vs(d$mean[d$clase == "Dentro"], d$mean[d$clase == "Sin motivo"], MIN_N) %>% mutate(comparacion = "Dentro vs. sin motivo"),
          tata_hl_vs(d$mean[d$clase == "Solo fuera"], d$mean[d$clase == "Sin motivo"], MIN_N) %>% mutate(comparacion = "Solo fuera vs. sin motivo"),
          tata_hl_vs(d$mean[d$clase == "Dentro"], d$mean[d$clase == "Solo fuera"], MIN_N) %>% mutate(comparacion = "Dentro vs. solo fuera")
        ) %>% mutate(n_dentro = sum(d$clase == "Dentro"), n_solo_fuera = sum(d$clase == "Solo fuera"), n_sin = sum(d$clase == "Sin motivo"))
      }) %>%
      ungroup() %>%
      mutate(motivo = motivo, ventana = paste0(minus(desde), " a ", minus(hasta)), preferida = preferida)
  }) %>%
  group_by(rep) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup() %>%
  relocate(rep, grupo, motivo, ventana, preferida, comparacion)

save_tsv(
  bind_rows(
    matrices %>% mutate(analisis = "matrices"),
    prevalencia %>% mutate(analisis = "prevalencia_cgi"),
    efecto %>% mutate(analisis = "efecto_HL"),
    calibracion %>% mutate(analisis = "calibracion_epd_solapamiento"),
    efecto_calib %>% mutate(analisis = "calibracion_epd_efecto_HL"),
    posicion %>% mutate(analisis = "posicion_ventana"),
    efecto_ventana %>% mutate(analisis = "efecto_HL_ventana")
  ) %>% relocate(analisis, motivo),
  "motivos_jaspar.tsv"
)

# --- graficos parte A --------------------------------------------------------

etiq_motivo <- function(m) {
  f <- MOTIVOS$familia[match(m, MOTIVOS$motivo)]
  paste0(m, if_else(f == m, "", paste0(" (", f, ")")))
}
niveles_motivo <- rev(etiq_motivo(MOTIVOS$motivo))

p_efecto <- efecto %>%
  mutate(motivo = factor(etiq_motivo(motivo), niveles_motivo), sig = sig_lab(p_adj)) %>%
  forest(motivo, grupo, col_grupo, "Diferencia de actividad media con vs. sin el motivo\n(Hodges-Lehmann, IC 95%)")
ggsave(file.path(out_dir, "motivos_jaspar_efecto.jpg"), p_efecto, width = 11, height = 9, units = "in")

p_calib <- efecto_calib %>%
  mutate(motivo = factor(motivo, rev(unique(motivo))), sig = sig_lab(p_adj)) %>%
  forest(motivo, grupo, col_grupo, "Diferencia de actividad media con vs. sin (Hodges-Lehmann, IC 95%)\nJASPAR en todo el fragmento, anotación de EPD y JASPAR en la ventana de EPD")
ggsave(file.path(out_dir, "motivos_jaspar_calibracion_epd.jpg"), p_calib, width = 11, height = 6.5, units = "in")

ventanas_plot <- posicion %>% mutate(motivo = factor(etiq_motivo(motivo), rev(niveles_motivo)))
p_pos <- hits %>%
  mutate(motivo = factor(etiq_motivo(motivo), rev(niveles_motivo))) %>%
  ggplot(aes(posicion, fill = hebra)) +
  geom_rect(data = ventanas_plot, aes(xmin = desde - 0.5, xmax = hasta + 0.5, ymin = -Inf, ymax = Inf, alpha = preferida),
    inherit.aes = FALSE, fill = "grey60") +
  geom_histogram(binwidth = 5, boundary = 0, position = "stack") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = col_hebra, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 0.1), labels = c(`TRUE` = "Ventana preferida", `FALSE` = "Ventana de mayor densidad (no preferida)"), name = NULL) +
  scale_x_continuous(breaks = c(-200, -150, -100, -50, 0), labels = minus) +
  facet_wrap(~motivo, scales = "free_y", ncol = 3) +
  th +
  labs(
    x = "Posición del centro de la coincidencia respecto del TSS (pb, bins de 5)", y = "Coincidencias",
    caption = paste0("Sombreado: ventana de ", VENTANA_ANCHO, " pb con más coincidencias; preferida si su fracción de coincidencias es ≥ ", ENRIQ_MIN,
                     " veces la de matrices de columnas permutadas y la uniforme (y ≥ 50 coincidencias).")
  )
ggsave(file.path(out_dir, "motivos_jaspar_posicion.jpg"), p_pos, width = 13, height = 12, units = "in")

p_ventana <- efecto_ventana %>%
  mutate(motivo = factor(paste0(etiq_motivo(motivo), "\n", ventana, if_else(preferida, " (preferida)", "")),
                         levels = rev(unique(paste0(etiq_motivo(posicion$motivo), "\n", minus(posicion$desde), " a ", minus(posicion$hasta),
                                                     if_else(posicion$preferida, " (preferida)", ""))))),
         sig = sig_lab(p_adj)) %>%
  forest(motivo, comparacion, c("Dentro vs. sin motivo" = "#D6741F", "Solo fuera vs. sin motivo" = "#14AFB2", "Dentro vs. solo fuera" = "grey30"),
         "Diferencia de actividad media (Hodges-Lehmann, IC 95%)", facets = grupo ~ rep)
ggsave(file.path(out_dir, "motivos_jaspar_ventana.jpg"), p_ventana, width = 14, height = 22, units = "in")

# =============================================================================
# PARTE B - k-meros (6-meros canonicos)
# =============================================================================

kmeros_todos <- mkAllStrings(c("A", "C", "G", "T"), K)
rc_de <- as.character(reverseComplement(DNAStringSet(kmeros_todos)))
canon_de <- pmin(kmeros_todos, rc_de)
canonicos <- sort(unique(canon_de))
stopifnot(length(canonicos) == 2080)

conteo <- oligonucleotideFrequency(fa, width = K) # promotores x 4096, solo hebra del promotor
# presencia del canonico en cualquier hebra = el k-mero o su complemento en la hebra del promotor
pres_k <- sapply(canonicos, function(km) {
  cols <- unique(c(km, as.character(reverseComplement(DNAString(km)))))
  rowSums(conteo[, cols, drop = FALSE]) > 0
})
rownames(pres_k) <- names(fa)

info_k <- tibble(
  kmero = canonicos,
  rc = as.character(reverseComplement(DNAStringSet(canonicos))),
  gc = str_count(canonicos, "[GC]") / K,
  n_cpg = str_count(canonicos, "CG"),
  prevalencia = colMeans(pres_k)
) %>%
  mutate(
    analizado = between(prevalencia, PREV_RANGO[1], PREV_RANGO[2]),
    # motivos ya analizados en la tesis (busqueda en el k-mero o su complemento)
    motivo_conocido = case_when(
      str_detect(kmero, "TATA[AT]A") | str_detect(rc, "TATA[AT]A") ~ "TATA",
      str_detect(kmero, "CCAAT") | str_detect(rc, "CCAAT") ~ "CCAAT",
      str_detect(kmero, "GGGCGG|GGCGGG") | str_detect(rc, "GGGCGG|GGCGGG") ~ "GC-box",
      str_detect(kmero, "^[CGT][CGT]CA[CGT][AT]$") | str_detect(rc, "^[CGT][CGT]CA[CGT][AT]$") ~ "Iniciador (BBCABW)",
      TRUE ~ NA_character_
    )
  )
kmeros_ok <- info_k$kmero[info_k$analizado]
message(length(kmeros_ok), " de ", length(canonicos), " 6-meros canonicos con prevalencia entre 1% y 99%")

if (file.exists(cache_kmeros)) {
  hl_k <- readRDS(cache_kmeros)
} else {
  t0 <- Sys.time()
  grupos_k <- data_g %>% select(rep, grupo, seq_id, mean) %>% group_split(rep, grupo)
  bloques <- split(kmeros_ok, ceiling(seq_along(kmeros_ok) / 40))
  tarea_k <- function(bloque) {
    map_dfr(grupos_k, function(d) {
      filas <- match(d$seq_id, rownames(bloque))
      map_dfr(colnames(bloque), function(km) {
        ind <- bloque[filas, km]
        tata_hl_vs(d$mean[ind], d$mean[!ind], MIN_N) %>%
          mutate(rep = d$rep[1], grupo = d$grupo[1], kmero = km, n = sum(ind), n_ref = sum(!ind))
      })
    })
  }
  cl <- parallel::makeCluster(max(1, min(12, parallel::detectCores() - 2)))
  parallel::clusterEvalQ(cl, suppressMessages(library(tidyverse)))
  parallel::clusterExport(cl, c("grupos_k", "tata_hl_vs", "MIN_N"))
  hl_k <- bind_rows(parallel::parLapplyLB(cl, map(bloques, ~ pres_k[, .x, drop = FALSE]), tarea_k))
  parallel::stopCluster(cl)
  saveRDS(hl_k, cache_kmeros)
  message("HL de k-meros: ", format(Sys.time() - t0))
}

kmeros_efecto <- hl_k %>%
  group_by(rep, grupo) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup() %>%
  left_join(info_k, by = "kmero") %>%
  mutate(grupo = factor(grupo, GRUPOS)) %>%
  relocate(rep, grupo, kmero, rc, gc, n_cpg, prevalencia, motivo_conocido, n, n_ref)
save_tsv(kmeros_efecto, "kmeros_efecto.tsv")

# (a) relacion del efecto con la composicion del k-mero
kmeros_composicion <- kmeros_efecto %>%
  group_by(rep, grupo) %>%
  summarise(
    n_kmeros = sum(!is.na(estimate)),
    rho_gc = cor(estimate, gc, method = "spearman", use = "complete.obs"),
    p_gc = cor.test(estimate, gc, method = "spearman", exact = FALSE)$p.value,
    rho_cpg = cor(estimate, n_cpg, method = "spearman", use = "complete.obs"),
    p_cpg = cor.test(estimate, n_cpg, method = "spearman", exact = FALSE)$p.value,
    n_sig = sum(p_adj < 0.05, na.rm = TRUE),
    pct_positivos_entre_sig_con_cg = 100 * mean(estimate[p_adj < 0.05 & n_cpg > 0] > 0, na.rm = TRUE),
    pct_positivos_entre_sig_sin_cg = 100 * mean(estimate[p_adj < 0.05 & n_cpg == 0] > 0, na.rm = TRUE),
    .groups = "drop"
  )
save_tsv(kmeros_composicion, "kmeros_composicion.tsv")

# (b, c) k-meros destacados dentro de cada grupo de isla CpG: BH < 0,05 y el
# mismo signo en ambas replicas; orden por el HL medio de las dos replicas
consistentes <- kmeros_efecto %>%
  filter(grupo != "Todos") %>%
  select(grupo, kmero, rep, estimate, conf.low, conf.high, p_adj) %>%
  pivot_wider(names_from = rep, values_from = c(estimate, conf.low, conf.high, p_adj), names_glue = "{.value}_{str_remove(rep, ' ')}") %>%
  filter(p_adj_Rep1 < 0.05, p_adj_Rep2 < 0.05, sign(estimate_Rep1) == sign(estimate_Rep2)) %>%
  mutate(hl_medio = (estimate_Rep1 + estimate_Rep2) / 2, direccion = if_else(hl_medio > 0, "Positivo", "Negativo")) %>%
  left_join(info_k %>% select(kmero, rc, gc, n_cpg, prevalencia, motivo_conocido), by = "kmero")
kmeros_top <- bind_rows(
  consistentes %>% mutate(conjunto = "Todos los k-meros"),
  consistentes %>% filter(n_cpg == 0) %>% mutate(conjunto = "Sin CG")
) %>%
  group_by(conjunto, grupo, direccion) %>%
  arrange(desc(abs(hl_medio)), .by_group = TRUE) %>%
  mutate(rango = row_number()) %>%
  filter(rango <= N_TOP) %>%
  ungroup()

# parecido con JASPAR: cada k-mero como matriz de un solo sitio vs. todas las
# matrices CORE vertebrados (correlacion de Pearson por columna, ambas
# hebras); se informa la de mayor puntaje
km_destacados <- unique(kmeros_top$kmero)
km_motivos <- map(km_destacados, ~ create_motif(.x, name = .x, type = "PCM"))
comparacion <- compare_motifs(c(km_motivos, jaspar),
  compare.to = seq_along(km_motivos), method = "PCC", tryRC = TRUE,
  min.overlap = K, max.p = 1, max.e = Inf
)
mejor_jaspar <- as_tibble(as.data.frame(comparacion)) %>%
  group_by(subject) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(kmero = subject, jaspar_id = target, jaspar_factor = jaspar_alt[match(target, map_chr(jaspar, ~ .x@name))],
            jaspar_score = score, jaspar_p = Pval)
kmeros_top <- kmeros_top %>%
  left_join(mejor_jaspar, by = "kmero") %>%
  relocate(conjunto, grupo, direccion, rango, kmero, rc, hl_medio, motivo_conocido, jaspar_factor, jaspar_id, jaspar_score, jaspar_p)
save_tsv(kmeros_top, "kmeros_top.tsv")

# --- graficos parte B --------------------------------------------------------

rotulos_volcano <- kmeros_top %>%
  filter(conjunto == "Todos los k-meros", rango <= 5) %>%
  select(grupo, kmero) %>%
  inner_join(kmeros_efecto %>% mutate(grupo = as.character(grupo)), by = c("grupo", "kmero")) %>%
  mutate(grupo = factor(grupo, GRUPOS))
p_volcano <- kmeros_efecto %>%
  filter(!is.na(estimate)) %>%
  ggplot(aes(estimate, -log10(p_adj), colour = factor(n_cpg))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(size = 0.9, alpha = 0.6) +
  ggrepel::geom_text_repel(data = rotulos_volcano, aes(label = kmero), size = 3, show.legend = FALSE, min.segment.length = 0, max.overlaps = 30) +
  scale_colour_manual(values = col_cpg, name = "CpG en el 6-mero") +
  facet_grid(grupo ~ rep) +
  th +
  labs(x = "Diferencia de actividad media con vs. sin el 6-mero (HL)", y = "−log10(p ajustado, BH)",
       caption = "Rótulos: los 5 6-meros con mayor efecto positivo y negativo, consistentes en ambas réplicas, en cada grupo de isla CpG.")
ggsave(file.path(out_dir, "kmeros_volcano.jpg"), p_volcano, width = 12, height = 13, units = "in")

# Volcano para ver los motivos (el anterior colorea por CpG y sirve para la
# composicion). Solo con / sin isla CpG; color por familia de motivo
# conocido, asignada por secuencia en el 6-mero o su complemento, en este
# orden de prioridad; el resto en gris. Rotulos: los 6-meros destacados
# (BH < 0,05 y mismo signo en ambas replicas) que pertenecen a una familia,
# mas los 5 de mayor efecto positivo y negativo sin familia.
FAMILIAS_K <- c(
  "TATA" = "TATA[AT]",
  "CCAAT (NF-Y)" = "CCAAT|CAATC", # CAATC: parte del sitio de NF-Y (CCAATCAG)
  "CRE (CREB/ATF)" = "TGACG|GACGTC",
  "YY1" = "AAATGG|AATGGC|AAAATG",
  "E-box ZEB/SNAI" = "CAGGT|ACCTG",
  "GC-box (SP/KLF)" = "GGGCGG|GGCGGG"
)
col_familia <- c(
  "TATA" = "#D6741F", "CCAAT (NF-Y)" = "#216869", "CRE (CREB/ATF)" = "#B5345A", "YY1" = "#7FB800",
  "E-box ZEB/SNAI" = "#3B6FB6", "GC-box (SP/KLF)" = "#A98BD0", "Otro" = "grey75"
)
familia_kmero <- function(k, rc) {
  fam <- rep("Otro", length(k))
  for (f in rev(names(FAMILIAS_K))) {
    fam[str_detect(k, FAMILIAS_K[[f]]) | str_detect(rc, FAMILIAS_K[[f]])] <- f
  }
  factor(fam, names(col_familia))
}
volc_m <- kmeros_efecto %>%
  filter(grupo != "Todos", !is.na(estimate)) %>%
  mutate(familia = familia_kmero(kmero, rc), grupo = factor(as.character(grupo), GRUPOS[-1]))
destacados_m <- consistentes %>%
  mutate(familia = familia_kmero(kmero, rc)) %>%
  group_by(grupo, familia == "Otro", direccion) %>%
  filter(familia != "Otro" | rank(-abs(hl_medio)) <= 5) %>%
  ungroup() %>%
  select(grupo, kmero) %>%
  mutate(grupo = as.character(grupo))
rotulos_m <- volc_m %>%
  mutate(grupo_chr = as.character(grupo)) %>%
  semi_join(destacados_m, by = c("grupo_chr" = "grupo", "kmero"))
p_volc_m <- ggplot(volc_m, aes(estimate, -log10(p_adj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(data = ~ filter(.x, familia == "Otro"), colour = "grey75", size = 0.7, alpha = 0.4) +
  geom_point(data = ~ filter(.x, familia != "Otro"), aes(colour = familia), size = 1.8, alpha = 0.9) +
  ggrepel::geom_text_repel(data = rotulos_m, aes(label = kmero, colour = familia), size = 4.6, fontface = "bold",
    show.legend = FALSE, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.4, seed = 1) +
  scale_colour_manual(values = col_familia, breaks = setdiff(names(col_familia), "Otro"), drop = FALSE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  facet_grid(grupo ~ rep) +
  th +
  labs(
    x = "Diferencia de actividad media con vs. sin el 6-mero (HL)", y = "−log10(p ajustado, BH)",
    caption = paste0(
      "Color: familia asignada por secuencia (6-mero o su complemento): TATA = TATAW; CCAAT = CCAAT o CAATC; CRE = TGACG o GACGTC;\n",
      "YY1 = AAATGG, AATGGC o AAAATG; E-box ZEB/SNAI = CAGGT; GC-box = GGGCGG o GGCGGG; gris: el resto.\n",
      "Rótulos: 6-meros consistentes en ambas réplicas (BH < 0,05, mismo signo) de alguna familia, más los 5 de mayor efecto de cada signo sin familia."
    )
  )
save_fig(p_volc_m, dir_novo, "kmeros_volcano_motivos", 14, 11)

# Grafico de puntos con IC de 6-meros destacados (consistentes en ambas
# replicas con isla CpG, dos o tres por familia), con el HL de cada replica
# con y sin isla CpG. Seleccion a mano a partir de kmeros_top.tsv.
SELECCION_K <- tribble(
  ~kmero_sel, ~familia,
  "TATATA", "TATA", "TATAAA", "TATA", "CTTATA", "TATA",
  "ATTGGC", "CCAAT (NF-Y)", "CCAATC", "CCAAT (NF-Y)", "CAATCA", "CCAAT (NF-Y)", "ACCAAT", "CCAAT (NF-Y)",
  "ACGTCA", "CRE (CREB/ATF)", "GACGTC", "CRE (CREB/ATF)",
  "AATGGC", "YY1", "AAATGG", "YY1", "AAAATG", "YY1",
  "CAGGTA", "E-box ZEB/SNAI", "ACAGGT", "E-box ZEB/SNAI", "CACCTG", "E-box ZEB/SNAI",
  "CCCGCC", "GC-box (SP/KLF)", "CCGCCC", "GC-box (SP/KLF)",
  "ACACAC", "Otros negativos", "CCATGG", "Otros negativos",
  "ACGAGA", "Sin familia (no confirmado)"
) %>%
  mutate(
    kmero = pmin(kmero_sel, as.character(reverseComplement(DNAStringSet(kmero_sel)))),
    familia = factor(familia, unique(familia))
  )
destacados_k <- kmeros_efecto %>%
  filter(grupo != "Todos") %>%
  inner_join(SELECCION_K, by = "kmero") %>%
  mutate(
    etiqueta = paste0(kmero, " / ", rc),
    etiqueta = factor(etiqueta, rev(unique(paste0(SELECCION_K$kmero, " / ", info_k$rc[match(SELECCION_K$kmero, info_k$kmero)])))),
    grupo = factor(as.character(grupo), GRUPOS[-1]),
    sig = sig_lab(p_adj)
  )
stopifnot(n_distinct(destacados_k$kmero) == nrow(SELECCION_K))
write_tsv(
  destacados_k %>% select(familia, kmero, rc, rep, grupo, n, n_ref, estimate, conf.low, conf.high, p, p_adj) %>% arrange(familia, kmero, grupo, rep),
  file.path(out_dir, "kmeros_destacados.tsv")
)
p_dest <- ggplot(destacados_k, aes(estimate, etiqueta, colour = rep, group = rep)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.6), linewidth = 0.9) +
  geom_point(aes(shape = sig), size = 2.4, fill = "white", position = position_dodge(width = 0.6)) +
  scale_shape_manual(values = c(`BH < 0,05` = 16, `n.s.` = 21), name = NULL) +
  scale_colour_manual(values = c("Rep 1" = "#216869", "Rep 2" = "#14AFB2"), name = NULL) +
  facet_grid(familia ~ grupo, scales = "free_y", space = "free_y", labeller = labeller(familia = label_wrap_gen(14))) +
  th +
  theme(strip.text.y = element_text(angle = 0, hjust = 0)) +
  labs(
    x = "Diferencia de actividad media con vs. sin el 6-mero\n(Hodges-Lehmann, IC 95%)", y = NULL,
    caption = "6-meros consistentes en ambas réplicas con isla CpG (BH < 0,05, mismo signo), dos o tres por familia. Rótulo: 6-mero / complemento reverso."
  )
# exploratorio: grafico suelto, no forma parte de la figura de la tesis
ggsave(file.path(out_dir, "kmeros_destacados.jpg"), p_dest, width = 12, height = 11, units = "in")

conc <- kmeros_efecto %>%
  select(grupo, kmero, n_cpg, rep, estimate) %>%
  pivot_wider(names_from = rep, values_from = estimate) %>%
  filter(!is.na(`Rep 1`), !is.na(`Rep 2`))
rho_conc <- conc %>%
  group_by(grupo) %>%
  summarise(label = sprintf("ρ = %.2f", cor(`Rep 1`, `Rep 2`, method = "spearman")), .groups = "drop")
p_conc <- ggplot(conc, aes(`Rep 1`, `Rep 2`, colour = factor(n_cpg))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 0.9, alpha = 0.6) +
  geom_text(data = rho_conc, aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.3, inherit.aes = FALSE, size = 4) +
  scale_colour_manual(values = col_cpg, name = "CpG en el 6-mero") +
  facet_wrap(~grupo, nrow = 1) +
  coord_equal() +
  th +
  labs(x = "HL, Rep 1", y = "HL, Rep 2")
ggsave(file.path(out_dir, "kmeros_concordancia.jpg"), p_conc, width = 14, height = 5.5, units = "in")

# efecto de cada 6-mero segun su numero de CpG, por grupo de isla CpG: en
# todos los promotores el efecto sube con los CpG (la composicion domina);
# dentro de cada grupo la relacion se aplana. HL medio de las dos replicas;
# R2 = variancia del efecto explicada por los CpG (0, 1, >= 2) + G+C, por
# replica. Los 6-meros con 2 y 3 CpG van juntos (hay un solo 6-mero canonico
# con 3: CGCGCG).
CPG_NIVELES <- c("0", "1", "≥ 2")
cpg_cat <- function(n) factor(if_else(n >= 2, "≥ 2", as.character(n)), CPG_NIVELES)
r2_comp <- kmeros_efecto %>%
  filter(!is.na(estimate)) %>%
  mutate(cpg = cpg_cat(n_cpg)) %>%
  group_by(rep, grupo) %>%
  summarise(r2 = summary(lm(estimate ~ cpg + gc))$r.squared, .groups = "drop") %>%
  group_by(grupo) %>%
  summarise(label = paste0("R² composición = ", sprintf("%.2f", r2[rep == "Rep 1"]), " / ", sprintf("%.2f", r2[rep == "Rep 2"]),
                           "\n(Rep 1 / Rep 2)") %>% str_replace_all("\\.", ","), .groups = "drop")
hl_medio_k <- kmeros_efecto %>%
  filter(!is.na(estimate)) %>%
  group_by(grupo, kmero, n_cpg) %>%
  summarise(hl_medio = mean(estimate), n_rep = n(), .groups = "drop") %>%
  filter(n_rep == 2) %>%
  mutate(cpg = cpg_cat(n_cpg))
n_por_cpg <- hl_medio_k %>% count(grupo, cpg)
p_cpg <- ggplot(hl_medio_k, aes(cpg, hl_medio)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_jitter(aes(colour = cpg), position = position_jitter(width = 0.25, seed = 1), size = 0.8, alpha = 0.3, show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5, linewidth = 0.6) +
  geom_text(data = n_por_cpg, aes(y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 3.5, colour = "grey30") +
  geom_text(data = r2_comp, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2, size = 3.8, inherit.aes = FALSE) +
  scale_colour_manual(values = set_names(col_cpg[c("0", "1", "2")], CPG_NIVELES)) +
  facet_wrap(~grupo, nrow = 1) +
  th +
  labs(
    x = "CpG en el 6-mero", y = "Diferencia de actividad media con vs. sin el 6-mero\n(HL, media de las dos réplicas)",
    caption = "Un punto por 6-mero canónico. R²: variancia del efecto explicada por el número de CpG y el contenido de G+C del 6-mero."
  )
save_fig(p_cpg, dir_novo, "kmeros_efecto_vs_cpg", 13, 5.5)

# =============================================================================
# Graficos con el formato de R2.2 (tata_actividad): YY1 y E-box tipo ZEB/SNAI
# =============================================================================
# motif_scatter() y violin_by_motif() son copia de las de
# R/14_motivos_cgi_actividad/motivos_cgi_actividad.R (mismo formato); los bins
# de 100 promotores por actividad media se arman con add_mean_sw_bins()
# (plot_helpers.R), por separado en cada replica y en cada universo.
# YY1 = coincidencia de MA0095.4 y CREB1 = de MA0018.5 (parte A, puntaje
# relativo >= 0,85).
# E-box tipo ZEB/SNAI = CAGGTR en cualquier hebra (agrupa los 6-meros
# CAGGTG/CACCTG y CAGGTA/TACCTG, negativos en la parte B).

presence_labels <- c("FALSE" = "Ausencia", "TRUE" = "Presencia")
motif_scatter <- function(data, motif, y_label, title) {
  data %>%
    group_by(mean_sw, rep) %>%
    summarise(prop = sum(.data[[motif]]) / n(), .groups = "drop") %>%
    ggplot(aes(as.numeric(mean_sw), prop)) +
    geom_point(col = "#14AFB2") +
    geom_smooth(col = "#216869") +
    theme_pubclean() +
    theme(text = element_text(size = 30), panel.spacing = unit(3, "lines")) +
    facet_wrap(~rep) +
    ylim(0, 1) +
    scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = c(0.03, 0.1))) +
    labs(x = "Bins de actividad media", y = y_label) +
    ggtitle(title)
}
violin_by_motif <- function(data, rep_id, motif, fill_label) {
  p <- ggviolin(data %>% filter(rep == rep_id),
    x = motif, y = "mean", fill = motif, draw_quantiles = 0.5,
    add = "median_q1q3", palette = c("#14AFB2", "#216869"), alpha = 0.7
  )
  p + stat_compare_means() +
    theme_pubclean() +
    scale_x_discrete(labels = presence_labels) +
    labs(fill = fill_label, x = NULL, y = "Actividad media") +
    theme(legend.position = "none", axis.ticks.x = element_blank(), text = element_text(size = 40)) +
    ylim(c(0.5, 6.2))
}

seq_chr <- set_names(as.character(fa), names(fa))
# CREB1 solo en su ventana preferida (la de mayor densidad de coincidencias,
# ver "posicion"; hoy -71 a -42): al menos una coincidencia centrada ahi
vent_creb1 <- posicion %>% filter(motivo == "CREB1")
ids_creb1_vent <- hits %>%
  filter(motivo == "CREB1", between(posicion, vent_creb1$desde, vent_creb1$hasta)) %>%
  distinct(seq_id) %>%
  pull(seq_id)
lab_vent_creb1 <- paste0(minus(vent_creb1$desde), " a ", minus(vent_creb1$hasta))
data_r22 <- data %>%
  left_join(prom_mot %>% select(seq_id, YY1, CREB1), by = "seq_id") %>%
  mutate(
    EBOX_ZEB = str_detect(seq_chr[seq_id], "CAGGT[AG]|[CT]ACCTG"),
    CREB1_VENT = seq_id %in% ids_creb1_vent
  )
universos_r22 <- list(
  todos = list(datos = data_r22, sufijo = ""),
  cgi = list(datos = filter(data_r22, CGI_frag %in% TRUE), sufijo = " (con isla CpG)")
)
motivos_r22 <- tribble(
  ~col, ~archivo, ~titulo, ~y_lab,
  "YY1", "yy1", "YY1", "Proporción de promotores\ncon YY1",
  "CREB1", "creb1", "CREB1", "Proporción de promotores\ncon CREB1",
  "CREB1_VENT", "creb1_ventana", paste0("CREB1 ", lab_vent_creb1), paste0("Proporción de promotores\ncon CREB1 en ", lab_vent_creb1),
  "EBOX_ZEB", "ebox_zeb", "E-box ZEB/SNAI", "Proporción de promotores\ncon E-box CAGGTR"
)
for (u in names(universos_r22)) {
  d <- add_mean_sw_bins(universos_r22[[u]]$datos)
  for (i in seq_len(nrow(motivos_r22))) {
    m <- motivos_r22[i, ]
    # el scatter de la E-box con isla CpG es el panel C de motivos_de_novo
    p_sc <- motif_scatter(d, m$col, m$y_lab, paste0(m$titulo, universos_r22[[u]]$sufijo))
    if (m$archivo == "ebox_zeb" && u == "cgi") {
      p_ebox_cgi <- p_sc
      save_fig(p_sc, dir_novo, paste0(m$archivo, "_scatter_", u), 12, 6.75)
    } else {
      ggsave(file.path(out_dir, paste0(m$archivo, "_scatter_", u, ".jpg")), p_sc, width = 12, height = 6.75, units = "in")
    }
    for (r in c("Rep 1", "Rep 2")) {
      ggsave(file.path(out_dir, paste0(m$archivo, "_violin_", str_remove(tolower(r), " "), "_", u, ".jpg")),
        violin_by_motif(d, r, m$col, m$titulo), width = 9, height = 6.75, units = "in")
    }
  }
}
# prevalencia de las dos definiciones (para el texto)
save_tsv(
  bind_rows(map(names(universos_r22), function(u) {
    universos_r22[[u]]$datos %>%
      group_by(rep) %>%
      summarise(n = n(), pct_yy1 = 100 * mean(YY1), pct_creb1 = 100 * mean(CREB1), pct_creb1_ventana = 100 * mean(CREB1_VENT), pct_ebox_zeb = 100 * mean(EBOX_ZEB), .groups = "drop") %>%
      mutate(universo = u, .before = 1)
  })),
  "yy1_ebox_zeb_prevalencia.tsv"
)

# =============================================================================
# Figura completa motivos_de_novo
# =============================================================================
# A: efecto segun CpG del 6-mero (ancho completo); B: volcano con / sin isla
# CpG, a todo el ancho y alto para que se lean los rotulos; abajo C: E-box
# ZEB/SNAI por bin de actividad (con isla CpG). Los 6-meros destacados
# (kmeros_destacados) quedan como grafico exploratorio aparte. El
# scatter de C se rearma con el mismo tamano de texto que el resto (en su
# archivo suelto queda con el formato de tata_actividad).
th_fig <- theme(text = element_text(size = 15), plot.caption = element_text(size = 11))
fig_novo <- wrap_plots(
  A = p_cpg + th_fig,
  B = p_volc_m + th_fig + theme(legend.text = element_text(size = 15)),
  C = p_ebox_cgi + theme(text = element_text(size = 15), panel.spacing = unit(1.5, "lines"), plot.title = element_text(size = 15)),
  design = "AAAA
BBBB
CC##", heights = c(5.5, 12, 6.5)
) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 24, face = "bold"))
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(dir_novo, paste0("motivos_de_novo.", ext)), fig_novo, width = 16, height = 24.5, units = "in")
}
message("Figura motivos_de_novo guardada en ", dir_novo)
