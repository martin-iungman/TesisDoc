# Busqueda de la TATA-box en los fragmentos de la library con la regla con
# la que EPDnew anoto la TATA-box, y analisis posicional de su efecto sobre
# la actividad. Compartido por R/30_tata_posicion/tata_posicion.R (figura
# tata_posicion) y R/30_tata_posicion/exploratorio_tata_posicion.qmd.
#
# Regla (reproduce TATA_EPD exactamente: 1653/1653, 0 falsos positivos):
#   - Puntaje: suma de los pesos W de la matriz de Bucher (1990) publicados
#     en https://epd.expasy.org/epd/promoter_elements/tata_old.php (la mejor
#     base de cada columna vale 0; maximo posible = 0).
#   - Posicion: la de la columna 4 de la matriz (posicion 0 de EPD, la
#     segunda T de TATA), en la orientacion del promotor, en ambas hebras.
#   - Canonica: misma hebra, columna 4 entre -31 y -26 del TSS (la pagina da
#     -36 a -20 como "region preferida", pero con esa ventana sobran 171
#     promotores que EPD no anota).
#   - Corte: -8.17 (EPD publica -8.16, pero los pesos publicados estan
#     redondeados a 2 decimales y con -8.16 se pierden 5 promotores anotados).
#
# Geometria del fragmento: 252 pb de -236 a +15 respecto del TSS de EPD
# (TSS = 0, en la posicion 237 de la secuencia orientada segun la hebra;
# verificado para los 20851 promotores FP). Solo promotores FP (EPDnew).
#
# Requiere tidyverse cargado; tata_load_library() requiere ademas
# rtracklayer y BSgenome.Hsapiens.UCSC.hg38.

TATA_TSS_POS <- 237L
TATA_FRAG_LEN <- 252L
TATA_ANCHOR_COL <- 4L
TATA_CANON <- c(-31, -26)
TATA_CORTE_EPD <- -8.17
TATA_WIN_ORIGIN <- -33 # ventanas de 10 pb alineadas para que -33 a -24 sea una de ellas
TATA_TOL <- 1e-9 # tolerancia numerica al comparar puntajes con un corte
TATA_CLASES <- c("Canónica", "Solo no canónica", "Solo hebra opuesta", "Sin motivo")

# Conteos de Bucher (1990, tabla 3; identicos a JASPAR MA0108.2 y TRANSFAC
# V$TATA_01 / M00252), solo para el logo.
TATA_COUNTS <- rbind(
  A = c(61, 16, 352, 3, 354, 268, 360, 222, 155, 56, 83, 82, 82, 68, 77),
  C = c(145, 46, 0, 10, 0, 0, 3, 2, 44, 135, 147, 127, 118, 107, 101),
  G = c(152, 18, 2, 2, 5, 0, 20, 44, 157, 150, 128, 128, 128, 139, 140),
  T = c(31, 309, 35, 374, 30, 121, 6, 121, 33, 48, 31, 52, 61, 75, 71)
)

# Pesos W de EPD, columnas = posiciones -3..11 de EPD
TATA_W <- rbind(
  A = c(-1.02, -3.05, 0.00, -4.61, 0.00, 0.00, 0.00, 0.00, -0.01, -0.94, -0.54, -0.48, -0.48, -0.74, -0.62),
  C = c(-0.28, -2.06, -5.22, -3.49, -5.17, -4.63, -4.12, -3.74, -1.13, -0.05, 0.00, -0.05, -0.11, -0.28, -0.40),
  G = c(0.00, -2.74, -4.28, -4.61, -3.77, -4.73, -2.65, -1.50, 0.00, 0.00, -0.09, 0.00, 0.00, 0.00, 0.00),
  T = c(-1.68, 0.00, -2.28, 0.00, -2.34, -0.52, -3.65, -0.37, -1.40, -0.97, -1.40, -0.82, -0.66, -0.54, -0.61)
)

# Fragmentos FP de la library orientados segun la hebra (getSeq respeta la
# hebra). `code`: matriz entera promotores x posicion (A=1, C=2, G=3, T=4,
# N = NA), con rownames = seq_id; `has_n`: fragmentos con alguna N.
tata_load_library <- function(bed = "data/raw/library.bed") {
  lib <- rtracklayer::import.bed(bed)
  lib <- lib[grepl("^FP", lib$name)]
  dna <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, lib)
  names(dna) <- lib$name
  seq_mat <- as.matrix(dna)
  code <- matrix(match(seq_mat, c("A", "C", "G", "T")), nrow = nrow(seq_mat), dimnames = list(rownames(seq_mat), NULL))
  list(lib = lib, dna = dna, seq_mat = seq_mat, code = code, has_n = rowSums(is.na(code)) > 0)
}

# Puntaje de cada posicion de inicio (filas = promotores, columnas = inicio
# 1..ncol(code) - ncol(w) + 1). `w` tiene filas A/C/G/T.
tata_score_matrix <- function(w, code) {
  k <- ncol(w)
  n_start <- ncol(code) - k + 1
  s <- matrix(0, nrow(code), n_start, dimnames = list(rownames(code), NULL))
  for (j in seq_len(k)) {
    cj <- code[, j:(j + n_start - 1), drop = FALSE]
    s <- s + matrix(w[cbind(as.vector(cj), j)], nrow(code))
  }
  s
}

tata_rev_comp <- function(w) {
  rc <- w[c("T", "G", "C", "A"), ncol(w):1]
  rownames(rc) <- c("A", "C", "G", "T")
  rc
}

# Puntajes en ambas hebras: list(misma, opuesta)
tata_scores <- function(code, w = TATA_W) {
  list(misma = tata_score_matrix(w, code), opuesta = tata_score_matrix(tata_rev_comp(w), code))
}

# Coincidencias >= min_score en ambas hebras. `pos` = posicion de la
# columna de referencia respecto del TSS: en la misma hebra esta
# anchor_col - 1 pb a la derecha del inicio; en la opuesta (matriz
# invertida) ncol(w) - anchor_col pb a la derecha del extremo izquierdo.
# `win` = inicio de la ventana de 10 pb que la contiene.
tata_hits <- function(s_misma, s_opuesta, min_score, exclude = character(),
                      motif_len = ncol(TATA_W), anchor_col = TATA_ANCHOR_COL,
                      tss_pos = TATA_TSS_POS, win_origin = TATA_WIN_ORIGIN) {
  to_tbl <- function(s, hebra, offset) {
    idx <- which(s >= min_score - TATA_TOL, arr.ind = TRUE)
    tibble(seq_id = rownames(s)[idx[, 1]], hebra = hebra, start_idx = idx[, 2], pos = idx[, 2] + offset - tss_pos, score = s[idx])
  }
  bind_rows(
    to_tbl(s_misma, "Misma hebra", anchor_col - 1L),
    to_tbl(s_opuesta, "Hebra opuesta", motif_len - anchor_col)
  ) %>%
    filter(!seq_id %in% exclude) %>%
    mutate(
      win = floor((pos - win_origin) / 10) * 10 + win_origin,
      hebra = factor(hebra, levels = c("Misma hebra", "Hebra opuesta"))
    )
}

# Una clase por promotor (prioridad: canonica > solo no canonica > solo
# hebra opuesta > sin motivo), a partir de las coincidencias que superan
# el corte (`hits_u`).
tata_classify <- function(hits_u, ids, canon = TATA_CANON) {
  canon_ids <- hits_u %>% filter(hebra == "Misma hebra", between(pos, canon[1], canon[2])) %>% pull(seq_id)
  misma <- hits_u %>% filter(hebra == "Misma hebra") %>% pull(seq_id)
  opuesta <- hits_u %>% filter(hebra == "Hebra opuesta") %>% pull(seq_id)
  tibble(seq_id = ids) %>%
    mutate(clase = case_when(
      seq_id %in% canon_ids ~ "Canónica",
      seq_id %in% misma ~ "Solo no canónica",
      seq_id %in% opuesta ~ "Solo hebra opuesta",
      TRUE ~ "Sin motivo"
    ) %>% factor(levels = TATA_CLASES))
}

# Mejor puntaje (misma hebra) con la columna de referencia entre from y to
tata_best_in_window <- function(s, from, to, anchor_col = TATA_ANCHOR_COL, tss_pos = TATA_TSS_POS) {
  cols <- (from:to) - (anchor_col - 1L) + tss_pos
  apply(s[, cols, drop = FALSE], 1, max)
}

# Fraccion de A/T de cada fragmento entre from y to (respecto del TSS)
tata_at_prox <- function(code, from = -50, to = -1, tss_pos = TATA_TSS_POS) {
  reg <- code[, (tss_pos + from):(tss_pos + to)]
  setNames(rowMeans(matrix(reg %in% c(1L, 4L), nrow = nrow(reg))), rownames(code))
}

# Hodges-Lehmann (diferencia de localizacion x - y) + IC 95% y p de Wilcoxon
tata_hl_vs <- function(x, y, min_n = 10) {
  if (length(x) < min_n || length(y) < min_n) {
    return(tibble(estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p = NA_real_))
  }
  wt <- wilcox.test(x, y, conf.int = TRUE, exact = FALSE)
  tibble(estimate = unname(wt$estimate), conf.low = wt$conf.int[1], conf.high = wt$conf.int[2], p = wt$p.value)
}

# Compara cada clase contra "Sin motivo" dentro de cada grupo (df con
# columnas mean y clase); BH dentro de cada grupo (familia = las 3
# comparaciones de un grupo).
tata_compare_classes <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(function(d, k) {
      ref <- d$mean[d$clase == "Sin motivo"]
      map_dfr(setdiff(TATA_CLASES, "Sin motivo"), function(cl) {
        tata_hl_vs(d$mean[d$clase == cl], ref) %>% mutate(clase = cl, n = sum(d$clase == cl), n_ref = length(ref))
      }) %>% mutate(p_adj = p.adjust(p, "BH"))
    }) %>%
    ungroup() %>%
    mutate(
      clase = factor(clase, levels = rev(setdiff(TATA_CLASES, "Sin motivo"))),
      efecto = case_when(
        p_adj < 0.05 & estimate > 0 ~ "Mayor actividad",
        p_adj < 0.05 & estimate < 0 ~ "Menor actividad",
        TRUE ~ "n.s."
      ) %>% factor(levels = c("Mayor actividad", "Menor actividad", "n.s."))
    )
}

# Perfil del efecto por posicion: para cada ventana de 10 pb y hebra, HL de
# los promotores con alguna coincidencia en esa ventana vs. los promotores
# sin ninguna coincidencia en `ref_hits` (por defecto, las mismas
# coincidencias: sin ninguna coincidencia en todo el fragmento, en ninguna
# hebra). BH dentro de cada grupo. `seg` agrupa ventanas contiguas para que
# las lineas no unan ventanas salteadas (solo ventanas con >= min_n).
tata_profile_effect <- function(act_df, hits_u, group_vars = "rep", min_n = 20, ref_hits = hits_u) {
  hit_ids <- unique(ref_hits$seq_id)
  wins <- hits_u %>% distinct(seq_id, hebra, win)
  act_df %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(function(d, k) {
      ref <- d$mean[!d$seq_id %in% hit_ids]
      n_tot <- nrow(d)
      d %>%
        inner_join(wins, by = "seq_id", relationship = "many-to-many") %>%
        group_by(hebra, win) %>%
        filter(n() >= min_n) %>%
        summarise(tata_hl_vs(mean, ref), n = n(), n_ref = length(ref), prop = n() / n_tot, .groups = "drop") %>%
        mutate(p_adj = p.adjust(p, "BH")) %>%
        arrange(hebra, win) %>%
        group_by(hebra) %>%
        mutate(seg = paste(hebra, cumsum(c(1, diff(win) != 10)))) %>%
        ungroup()
    }) %>%
    ungroup() %>%
    mutate(win_center = win + 4.5, sig = p_adj < 0.05)
}

# Permutaciones de las columnas de la matriz (control de composicion),
# distintas entre si y de la identidad. La columna de referencia sigue
# siendo la 4a de cada matriz permutada.
tata_permutations <- function(n = 10, seed = 2024, motif_len = ncol(TATA_W)) {
  set.seed(seed)
  perms <- list()
  while (length(perms) < n) {
    p <- sample(motif_len)
    if (!identical(p, seq_len(motif_len)) && !any(map_lgl(perms, ~ identical(.x, p)))) perms[[length(perms) + 1]] <- p
  }
  set_names(perms, paste0("perm", seq_along(perms)))
}
