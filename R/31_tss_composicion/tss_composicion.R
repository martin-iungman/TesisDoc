# EXPLORATORIO - composicion de bases en el TSS y actividad.
#
# Punto 1: ¿el efecto de las caracteristicas del TSS (dinucleotido -1/0
# CG/CA/TG/TA, TSS no canonico, Inr fuerte YCASW, TCT; columnas *_TSS de
# prom_df) depende de la isla CpG (CGI_frag)? Efecto = Hodges-Lehmann (IC
# 95%) de la actividad media con vs. sin la caracteristica, en todos los
# promotores y por separado con / sin isla CpG; prevalencia con / sin isla
# CpG (OR y p de Fisher); y modelo lineal del percentil de actividad
# ajustando por CGI_frag, TATA_EPD, CCAAT_EPD y GCbox_EPD.
#
# Punto 2: perfil posicional de -10 a +10 (bases y dinucleotidos (p, p+1)):
# HL de tener esa base / dinucleotido en esa posicion vs. no tenerlo, BH
# dentro de replica x grupo. Control de composicion: el mismo perfil en una
# ventana alejada (-130 a -110) y un modelo por posicion
# percentil ~ base + GC_local (G+C de las 10 bases de +-5 pb, sin contar la
# propia posicion, para que GC_local no incluya a la base evaluada).
#
# Punto 3: iniciador con el consenso humano (INR_humano, BBCA+1BW de Vo
# ngoc et al. 2017) vs. YCASW: prevalencia con / sin isla CpG, solapamiento,
# efecto por grupo y modelo crudo / ajustado.
# Punto 4: interaccion TATA-box x iniciador: efecto de TATA_EPD con y sin
# iniciador, las 4 combinaciones vs. "ninguno" (con la suma solo TATA +
# solo iniciador como referencia de aditividad) y el termino de
# interaccion del modelo lineal; en todos y por grupo de isla CpG.
# Punto 5: TCT y genes de proteinas ribosomales citoplasmaticas (RP): cuanto
# del efecto del TCT se explica por los promotores RP.
#
# Geometria: fragmento de 252 pb de -236 a +15, TSS de EPD (= 0, sin +1
# separado) en la posicion 237 de la secuencia orientada (TATA_TSS_POS en
# R/functions/tata_scan.R). Se verifica contra TSS_seq y las columnas
# *_TSS de prom_df antes de seguir.
#
# Figuras de la tesis (numero vigente en docs/mapping_figuras.csv), escritas
# via fig_dir():
#   - perfil_efectos_promotor: tss_perfil_efectos_promotor_hl.jpg (HL, sin ajuste por G+C) y su tabla.
#   - tss_composicion: tss_perfil_YR_todos.jpg y tss_modelo.jpg (en ese orden),
#     con tss_perfil_bases.tsv, tss_perfil_dinucleotidos.tsv y tss_modelo.tsv.
# Todo lo demas (graficos y tablas) es exploratorio y va a
# figures/exploratorio_tss_composicion/.
# Requiere: data/processed/activity_stats_highconf.tsv,
# data/processed/prom_df.tsv y la secuencia de la library
# (../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
library(Biostrings)
source("R/functions/tata_scan.R") # TATA_TSS_POS, TATA_FRAG_LEN, tata_hl_vs()
source("R/functions/plot_helpers.R") # thesis_clr
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename # Biostrings/S4Vectors la enmascara

out_dir <- "figures/exploratorio_tss_composicion"
dir_perfil_efectos <- fig_dir("perfil_efectos_promotor") # figura de la tesis
dir_tss <- fig_dir("tss_composicion") # figura de la tesis
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
path_fa <- "../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa"

TSS_POS <- TATA_TSS_POS
POS_TSS <- -10:10
POS_CONTROL <- -130:-110
GC_FLANCO <- 5
MIN_N <- 10
BASES <- c("A", "C", "G", "T")
SETS_BASE <- list(A = "A", C = "C", G = "G", T = "T")
SETS_YR <- list(Y = c("C", "T"), R = c("A", "G"))
DINUCS <- as.vector(outer(BASES, BASES, paste0))
SETS_DINUC <- set_names(as.list(DINUCS), DINUCS)
SET_YR_DINUC <- list(YR = c("CA", "CG", "TA", "TG"))

FEATS <- c(
  CG_TSS = "CG en −1/0", CA_TSS = "CA en −1/0", TG_TSS = "TG en −1/0", TA_TSS = "TA en −1/0",
  other_TSS = "No canónico", INR_strong_TSS = "Inr fuerte (YCASW)", TCT_TSS = "TCT"
)
GRUPOS <- c("Todos", "Con isla CpG", "Sin isla CpG")

col_alto <- "#D6741F"
col_bajo <- "#7FB800"
col_grupo <- c("Todos" = "grey30", "Con isla CpG" = "#216869", "Sin isla CpG" = "#14AFB2")
# mismos colores por base que perfil_nucleotidos_library
# (R/32_perfil_nucleotidos/perfil_nucleotidos.R), para que sean comparables
col_bases <- c(A = "#7FB800", C = thesis_clr, G = "#D6741F", T = "#B5345A")
col_efecto <- scale_color_manual(values = c("Mayor actividad" = col_alto, "Menor actividad" = col_bajo, "n.s." = "grey60"), drop = FALSE)
th <- theme_pubclean(base_size = 13) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")
minus <- function(x) sub("^-", "−", format(x, trim = TRUE))
breaks_pares <- function(l) seq(2 * ceiling(l[1] / 2), floor(l[2]), by = 2)
save_tsv <- function(x, name) write_tsv(x, file.path(out_dir, name))
efecto_lab <- function(estimate, p_adj) {
  case_when(
    p_adj < 0.05 & estimate > 0 ~ "Mayor actividad",
    p_adj < 0.05 & estimate < 0 ~ "Menor actividad",
    TRUE ~ "n.s."
  ) %>% factor(levels = c("Mayor actividad", "Menor actividad", "n.s."))
}

# --- Datos -------------------------------------------------------------------

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fa <- readDNAStringSet(path_fa)
fa <- fa[names(fa) %in% prom_df$seq_id]
stopifnot(all(width(fa) == TATA_FRAG_LEN), all(prom_df$seq_id %in% names(fa)))
seq_mat <- as.matrix(fa) # promotores x posicion, rownames = seq_id

# Verificacion de la geometria: el 5-mero -2..+2 debe ser TSS_seq y el
# dinucleotido -1/0 debe reproducir CA/CG/TA/TG_TSS.
chk <- prom_df %>%
  mutate(
    seq5 = apply(seq_mat[seq_id, TSS_POS + (-2:2)], 1, paste, collapse = ""),
    din = paste0(seq_mat[seq_id, TSS_POS - 1], seq_mat[seq_id, TSS_POS])
  )
verif <- tibble(
  chequeo = c("TSS_seq (-2..+2)", "CA_TSS", "CG_TSS", "TA_TSS", "TG_TSS"),
  coincidencias = c(
    sum(chk$seq5 == chk$TSS_seq, na.rm = TRUE),
    sum((chk$din == "CA") == chk$CA_TSS, na.rm = TRUE), sum((chk$din == "CG") == chk$CG_TSS, na.rm = TRUE),
    sum((chk$din == "TA") == chk$TA_TSS, na.rm = TRUE), sum((chk$din == "TG") == chk$TG_TSS, na.rm = TRUE)
  ),
  n = nrow(chk)
)
print(verif)
stopifnot(all(verif$coincidencias == verif$n))

# Iniciador humano (Vo ngoc et al. 2017, Genes Dev 31:6): BBCA+1BW, con la
# A del consenso en su +1 = nuestro 0 (no tienen posicion 0). En nuestras
# coordenadas: B(-3) B(-2) C(-1) A(0) B(+1) W(+2); B = C/G/T, W = A/T.
INR_HUMANO_POS <- -3:2
prom_df <- prom_df %>%
  mutate(
    inr_seq = apply(seq_mat[seq_id, TSS_POS + INR_HUMANO_POS], 1, paste, collapse = ""),
    INR_humano = str_detect(inr_seq, "^[CGT][CGT]CA[CGT][AT]$")
  )

data <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE) %>%
  inner_join(prom_df, by = c("seq_id", "name")) %>%
  group_by(rep) %>%
  mutate(act_pct = 100 * percent_rank(mean)) %>%
  ungroup()

# todos + con / sin isla CpG (la referencia de cada comparacion es la del mismo grupo)
data_g <- bind_rows(
  data %>% mutate(grupo = "Todos"),
  data %>% filter(!is.na(CGI_frag)) %>% mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG"))
) %>% mutate(grupo = factor(grupo, levels = GRUPOS))

# =============================================================================
# PUNTO 1 - caracteristicas del TSS x isla CpG
# =============================================================================

# HL de la actividad media con vs. sin cada caracteristica (columnas
# logicas `feats`), por replica x `group_vars`; BH dentro de cada grupo.
efecto_feats <- function(df, feats, group_vars = c("rep", "grupo")) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(function(d, k) {
      map_dfr(feats, function(f) {
        ind <- d[[f]]
        ok <- !is.na(ind)
        tata_hl_vs(d$mean[ok & ind], d$mean[ok & !ind], MIN_N) %>%
          mutate(caracteristica = f, n = sum(ok & ind), n_ref = sum(ok & !ind), prop = mean(ind[ok]))
      }) %>% mutate(p_adj = p.adjust(p, "BH"))
    }) %>%
    ungroup() %>%
    mutate(efecto = efecto_lab(estimate, p_adj)) %>%
    relocate(all_of(group_vars), caracteristica, n, n_ref, prop)
}

# prevalencia con / sin isla CpG y OR de Fisher, por universo
prevalencia_feats <- function(universos, feats) {
  universos %>%
    group_by(universo) %>%
    group_modify(function(d, k) {
      map_dfr(feats, function(f) {
        d2 <- d[!is.na(d[[f]]), ]
        ft <- fisher.test(table(factor(d2[[f]], c(TRUE, FALSE)), factor(d2$CGI_frag, c(TRUE, FALSE))))
        tibble(
          caracteristica = f,
          n_con_cgi = sum(d2$CGI_frag), n_feat_con_cgi = sum(d2[[f]] & d2$CGI_frag),
          n_sin_cgi = sum(!d2$CGI_frag), n_feat_sin_cgi = sum(d2[[f]] & !d2$CGI_frag),
          odds_ratio = unname(ft$estimate), or_conf.low = ft$conf.int[1], or_conf.high = ft$conf.int[2], p_fisher = ft$p.value
        ) %>% mutate(pct_con_cgi = 100 * n_feat_con_cgi / n_con_cgi, pct_sin_cgi = 100 * n_feat_sin_cgi / n_sin_cgi, .after = n_feat_sin_cgi)
      })
    }) %>%
    ungroup() %>%
    mutate(p_fisher_adj = p.adjust(p_fisher, "BH"))
}

# modelos lineales del percentil de actividad, por replica y caracteristica
modelos_lab <- c("Crudo", "Ajustado (CGI + TATA + CCAAT + GC-box)", "Ajustado + interacción con CGI")
modelos_feats <- function(df, feats, group_vars = "rep", interaccion = TRUE) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(function(d, k) {
      map_dfr(feats, function(f) {
        d2 <- d %>% mutate(feat = .data[[f]])
        bind_rows(
          broom::tidy(lm(act_pct ~ feat, data = d2), conf.int = TRUE) %>% mutate(modelo = modelos_lab[1]),
          broom::tidy(lm(act_pct ~ feat + CGI_frag + TATA_EPD + CCAAT_EPD + GCbox_EPD, data = d2), conf.int = TRUE) %>% mutate(modelo = modelos_lab[2]),
          if (interaccion) broom::tidy(lm(act_pct ~ feat * CGI_frag + TATA_EPD + CCAAT_EPD + GCbox_EPD, data = d2), conf.int = TRUE) %>% mutate(modelo = modelos_lab[3])
        ) %>% mutate(caracteristica = f)
      })
    }) %>%
    ungroup() %>%
    mutate(term = str_replace_all(term, "TRUE", "")) %>%
    relocate(all_of(group_vars), caracteristica, modelo, term)
}

efecto <- efecto_feats(data_g, names(FEATS))
save_tsv(efecto, "tss_efecto_cgi.tsv")

# prevalencia con / sin isla CpG: library completa y promotores con actividad de cada replica
prev_universos <- bind_rows(
  prom_df %>% mutate(universo = "Library (promotores FP)"),
  data %>% mutate(universo = paste0("Con actividad (", rep, ")"))
) %>% filter(!is.na(CGI_frag))
prevalencia <- prevalencia_feats(prev_universos, names(FEATS))
save_tsv(prevalencia, "tss_prevalencia_cgi.tsv")

modelo <- modelos_feats(data, names(FEATS))
write_tsv(modelo, file.path(dir_tss, "tss_modelo.tsv"))

# grafico de puntos con IC (formato de las clases de TATA)
feat_levels <- rev(unname(FEATS))
panel_efecto <- efecto %>%
  mutate(caracteristica = factor(FEATS[caracteristica], levels = feat_levels)) %>%
  ggplot(aes(estimate, caracteristica, color = grupo, shape = p_adj < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.7), linewidth = 0.9) +
  geom_point(size = 2.6, fill = "white", position = position_dodge(width = 0.7)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "BH < 0,05", `FALSE` = "n.s."), name = NULL) +
  scale_color_manual(values = col_grupo, name = NULL) +
  guides(color = guide_legend(reverse = FALSE)) +
  facet_grid(. ~ rep) +
  th +
  labs(x = "Diferencia de actividad media con vs. sin la característica\n(Hodges-Lehmann, IC 95%)", y = NULL)
ggsave(file.path(out_dir, "tss_efecto_cgi.jpg"), panel_efecto, width = 11, height = 7, units = "in")

panel_modelo <- modelo %>%
  filter(term == "feat", modelo != modelos_lab[3], caracteristica != "TCT_TSS") %>% # el TCT se ve en tss_tct_rp.jpg
  mutate(caracteristica = factor(FEATS[caracteristica], levels = setdiff(feat_levels, FEATS[["TCT_TSS"]]))) %>%
  ggplot(aes(estimate, caracteristica, color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.6), linewidth = 0.9) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("grey55", "#216869"), name = NULL) +
  facet_grid(. ~ rep) +
  th +
  labs(x = "Coeficiente de la característica sobre el percentil de actividad\n(modelo lineal, IC 95%)", y = NULL)
ggsave(file.path(dir_tss, "tss_modelo.jpg"), panel_modelo, width = 11, height = 6, units = "in")

# =============================================================================
# PUNTO 2 - perfil posicional
# =============================================================================

# Efecto de un indicador: HL sobre la actividad media + modelos del
# percentil (crudo y ajustado por GC_local).
efecto_ind <- function(mean, pct, ind, gc = NULL) {
  out <- tata_hl_vs(mean[ind], mean[!ind], MIN_N) %>% mutate(n = sum(ind), n_ref = sum(!ind), prop = mean(ind))
  if (is.null(gc)) return(out)
  if (sum(ind) < MIN_N || sum(!ind) < MIN_N) {
    return(out %>% mutate(lm_crudo_est = NA_real_, lm_crudo_p = NA_real_, lm_gc_est = NA_real_,
                          lm_gc_low = NA_real_, lm_gc_high = NA_real_, lm_gc_p = NA_real_, gc_coef = NA_real_))
  }
  m0 <- coef(summary(lm(pct ~ ind)))
  m1 <- lm(pct ~ ind + gc)
  s1 <- coef(summary(m1))
  ci <- confint(m1)["indTRUE", ]
  out %>% mutate(
    lm_crudo_est = m0["indTRUE", 1], lm_crudo_p = m0["indTRUE", 4],
    lm_gc_est = s1["indTRUE", 1], lm_gc_low = ci[1], lm_gc_high = ci[2], lm_gc_p = s1["indTRUE", 4],
    gc_coef = s1["gc", 1]
  )
}

# bases (A/C/G/T) y Y/R en cada posicion; GC_local = fraccion G+C de +-5 pb sin la propia posicion
perfil_bases <- function(d, posiciones, ventana) {
  M <- seq_mat[d$seq_id, , drop = FALSE]
  map_dfr(posiciones, function(pos) {
    col <- TSS_POS + pos
    b <- M[, col]
    S <- M[, setdiff(col + (-GC_FLANCO:GC_FLANCO), col), drop = FALSE]
    gc <- rowMeans(S == "G" | S == "C")
    ok <- b %in% BASES
    bind_rows(
      imap_dfr(SETS_BASE, ~ efecto_ind(d$mean[ok], d$act_pct[ok], b[ok] %in% .x, gc[ok]) %>% mutate(base = .y, conjunto = "Base")),
      imap_dfr(SETS_YR, ~ efecto_ind(d$mean[ok], d$act_pct[ok], b[ok] %in% .x, gc[ok]) %>% mutate(base = .y, conjunto = "Y/R"))
    ) %>% mutate(posicion = pos, ventana = ventana)
  })
}

# dinucleotidos en (p, p+1)
perfil_dinuc <- function(d, posiciones, sets, ventana, conjunto) {
  M <- seq_mat[d$seq_id, , drop = FALSE]
  map_dfr(posiciones, function(pos) {
    col <- TSS_POS + pos
    ok <- M[, col] %in% BASES & M[, col + 1] %in% BASES
    dn <- paste0(M[ok, col], M[ok, col + 1])
    imap_dfr(sets, ~ efecto_ind(d$mean[ok], d$act_pct[ok], dn %in% .x) %>% mutate(dinucleotido = .y)) %>%
      mutate(posicion = pos, conjunto = conjunto, ventana = ventana)
  })
}

# BH dentro de replica x grupo x ventana x conjunto
ajustar_bh <- function(res) {
  res %>%
    group_by(rep, grupo, ventana, conjunto) %>%
    mutate(p_adj = p.adjust(p, "BH"), across(any_of("lm_gc_p"), ~ p.adjust(.x, "BH"), .names = "lm_gc_p_adj")) %>%
    ungroup() %>%
    mutate(efecto = efecto_lab(estimate, p_adj))
}

ventana_tss <- "TSS (−10 a +10)"
ventana_control <- "Control (−130 a −110)"

# Los ~3500 Wilcoxon con IC tardan ~45 min en serie: se reparten por
# replica x grupo x posicion en un cluster local y el resultado crudo (sin
# BH) se guarda en data/processed/ para retocar graficos sin recalcular.
# Borrar el .rds para forzar el recalculo.
cache_perfiles <- "data/processed/tss_composicion_perfiles_v1.rds"
if (file.exists(cache_perfiles)) {
  crudo <- readRDS(cache_perfiles)
} else {
  t0 <- Sys.time()
  d_split <- data_g %>% select(rep, grupo, seq_id, mean, act_pct) %>% group_split(rep, grupo)
  tareas <- bind_rows(
    expand_grid(g = seq_along(d_split), posicion = POS_TSS, tipo = "bases_tss"),
    expand_grid(g = seq_along(d_split), posicion = POS_CONTROL, tipo = "bases_control"),
    expand_grid(g = seq_along(d_split), posicion = head(POS_TSS, -1), tipo = "dinuc_tss"),
    expand_grid(g = seq_along(d_split), posicion = head(POS_CONTROL, -1), tipo = "yr_control")
  )
  correr_tarea <- function(i) {
    t <- tareas[i, ]
    d <- d_split[[t$g]]
    res <- switch(t$tipo,
      bases_tss = perfil_bases(d, t$posicion, ventana_tss),
      bases_control = perfil_bases(d, t$posicion, ventana_control),
      dinuc_tss = bind_rows(
        perfil_dinuc(d, t$posicion, SETS_DINUC, ventana_tss, "Dinucleótido"),
        perfil_dinuc(d, t$posicion, SET_YR_DINUC, ventana_tss, "YR (CA/CG/TA/TG)")
      ),
      yr_control = perfil_dinuc(d, t$posicion, SET_YR_DINUC, ventana_control, "YR (CA/CG/TA/TG)")
    )
    res %>% mutate(rep = d$rep[1], grupo = d$grupo[1], tipo = if_else(str_starts(t$tipo, "bases"), "bases", "dinuc"))
  }
  cl <- parallel::makeCluster(max(1, min(12, parallel::detectCores() - 2)))
  parallel::clusterEvalQ(cl, suppressMessages(library(tidyverse)))
  parallel::clusterExport(cl, c(
    "tareas", "d_split", "seq_mat", "perfil_bases", "perfil_dinuc", "efecto_ind", "tata_hl_vs",
    "TSS_POS", "GC_FLANCO", "MIN_N", "BASES", "SETS_BASE", "SETS_YR", "SETS_DINUC", "SET_YR_DINUC",
    "ventana_tss", "ventana_control"
  ))
  crudo <- bind_rows(parallel::parLapplyLB(cl, seq_len(nrow(tareas)), correr_tarea))
  parallel::stopCluster(cl)
  saveRDS(crudo, cache_perfiles)
  message("Perfiles: ", format(Sys.time() - t0))
}

perfil <- crudo %>%
  filter(tipo == "bases") %>%
  select(where(~ !all(is.na(.x)))) %>%
  ajustar_bh() %>%
  relocate(rep, grupo, ventana, conjunto, posicion, base)
dinuc <- crudo %>%
  filter(tipo == "dinuc") %>%
  select(where(~ !all(is.na(.x)))) %>%
  ajustar_bh() %>%
  relocate(rep, grupo, ventana, conjunto, posicion, dinucleotido)

write_tsv(perfil %>% filter(ventana == ventana_tss), file.path(dir_tss, "tss_perfil_bases.tsv"))
write_tsv(dinuc %>% filter(ventana == ventana_tss), file.path(dir_tss, "tss_perfil_dinucleotidos.tsv"))
save_tsv(
  bind_rows(perfil %>% filter(ventana == ventana_control), dinuc %>% filter(ventana == ventana_control) %>% rename(base = dinucleotido)),
  "tss_perfil_control.tsv"
)

# --- (a) heatmaps posicion x base -------------------------------------------

lim_hl <- max(abs(perfil$estimate[perfil$conjunto == "Base"]), na.rm = TRUE)
escala_hl <- scale_fill_gradient2(low = col_bajo, mid = "white", high = col_alto, midpoint = 0,
                                  limits = c(-lim_hl, lim_hl), oob = scales::squish, name = "HL")
heatmap_pos <- function(res, y, x_lab, lims = escala_hl) {
  ggplot(res, aes(posicion, {{ y }}, fill = estimate)) +
    geom_tile(color = "white") +
    geom_text(aes(label = if_else(p_adj < 0.05, "*", "")), size = 4, vjust = 0.75) +
    lims +
    scale_x_continuous(breaks = breaks_pares, labels = minus, expand = c(0, 0)) +
    facet_grid(grupo ~ rep) +
    th +
    theme(legend.position = "right", panel.grid = element_blank()) +
    labs(x = x_lab, y = NULL, caption = "* BH < 0,05 (dentro de réplica × grupo × ventana). Color: HL de la actividad media, tener vs. no tener.")
}
p_heat_tss <- perfil %>%
  filter(ventana == ventana_tss, conjunto == "Base") %>%
  mutate(base = factor(base, rev(BASES))) %>%
  heatmap_pos(base, "Posición respecto del TSS (TSS = 0)") +
  geom_vline(xintercept = c(-1.5, 0.5), linewidth = 0.3, color = "grey30")
ggsave(file.path(out_dir, "tss_heatmap_bases.jpg"), p_heat_tss, width = 12, height = 7.5, units = "in")

p_heat_ctl <- perfil %>%
  filter(ventana == ventana_control, conjunto == "Base") %>%
  mutate(base = factor(base, rev(BASES))) %>%
  heatmap_pos(base, "Posición respecto del TSS (ventana control, misma escala de color que el TSS)")
ggsave(file.path(out_dir, "tss_heatmap_control.jpg"), p_heat_ctl, width = 12, height = 7.5, units = "in")

lim_dn <- max(abs(dinuc$estimate[dinuc$conjunto == "Dinucleótido"]), na.rm = TRUE)
p_heat_dn <- dinuc %>%
  filter(conjunto == "Dinucleótido") %>%
  mutate(dinucleotido = factor(dinucleotido, rev(DINUCS))) %>%
  heatmap_pos(dinucleotido, "Posición p del dinucleótido (p, p+1) respecto del TSS",
              scale_fill_gradient2(low = col_bajo, mid = "white", high = col_alto, midpoint = 0,
                                   limits = c(-lim_dn, lim_dn), oob = scales::squish, name = "HL")) +
  geom_vline(xintercept = c(-1.5, -0.5), linewidth = 0.3, color = "grey30")
ggsave(file.path(out_dir, "tss_heatmap_dinucleotidos.jpg"), p_heat_dn, width = 12, height = 13, units = "in")

# --- (b, c) perfiles de lineas: Y en p, R en p, YR en (p, p+1) ---------------

perfil_yr <- bind_rows(
  perfil %>% filter(conjunto == "Y/R") %>% mutate(serie = paste(base, "en p")),
  dinuc %>% filter(conjunto == "YR (CA/CG/TA/TG)") %>% mutate(serie = "YR en (p, p+1)")
) %>% mutate(serie = factor(serie, c("Y en p", "R en p", "YR en (p, p+1)")))
lim_yr <- range(c(perfil_yr$conf.low, perfil_yr$conf.high), na.rm = TRUE)
# Y violeta claro, R oliva/dorado (promedio RGB de la A y la G de col_bases), YR negro
col_serie <- scale_color_manual(values = c("Y en p" = "#A98BD0", "R en p" = "#AA960F", "YR en (p, p+1)" = "black"), name = NULL)
lineas_pos <- function(res, x_lab, marcas = NULL, facets = grupo ~ rep, ylim = lim_yr) {
  ggplot(res, aes(posicion, estimate, color = serie)) +
    {if (!is.null(marcas)) geom_vline(xintercept = marcas, linetype = "dotted", color = "grey40")} +
    geom_hline(yintercept = 0, color = "grey50") +
    geom_linerange(aes(ymin = conf.low, ymax = conf.high), alpha = 0.6) +
    geom_line(alpha = 0.7) +
    geom_point(aes(shape = p_adj < 0.05), size = 1.8, fill = "white") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "BH < 0,05", `FALSE` = "n.s."), name = NULL) +
    col_serie +
    scale_x_continuous(breaks = breaks_pares, labels = minus) +
    coord_cartesian(ylim = ylim) +
    facet_grid(facets) +
    th +
    labs(x = x_lab, y = "Diferencia de actividad media\n(HL, IC 95%)")
}
p_yr_tss <- perfil_yr %>%
  filter(ventana == ventana_tss) %>%
  lineas_pos("Posición p respecto del TSS (TSS = 0; punteadas en −1 y 0)", marcas = c(-1, 0))
ggsave(file.path(out_dir, "tss_perfil_YR.jpg"), p_yr_tss, width = 12, height = 9, units = "in")
# version solo con todos los promotores (sin el estrato por isla CpG)
p_yr_tss_todos <- perfil_yr %>%
  filter(ventana == ventana_tss, grupo == "Todos") %>%
  lineas_pos("Posición p respecto del TSS (TSS = 0; punteadas en −1 y 0)", marcas = c(-1, 0), facets = . ~ rep, ylim = NULL)
ggsave(file.path(dir_tss, "tss_perfil_YR_todos.jpg"), p_yr_tss_todos, width = 12, height = 4.5, units = "in")
p_yr_ctl <- perfil_yr %>%
  filter(ventana == ventana_control) %>%
  lineas_pos("Posición p respecto del TSS (ventana control; misma escala que el TSS)")
ggsave(file.path(out_dir, "tss_perfil_YR_control.jpg"), p_yr_ctl, width = 12, height = 9, units = "in")

# --- efecto de cada base por posicion (todos los promotores) ----------------
p_efecto_bases <- perfil %>%
  filter(ventana == ventana_tss, conjunto == "Base", grupo == "Todos") %>%
  ggplot(aes(posicion, estimate, colour = base)) +
  geom_vline(xintercept = c(-1, 0), linetype = "dotted", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5) +
  geom_line(alpha = 0.8) +
  geom_point(aes(shape = p_adj < 0.05), size = 2, fill = "white") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  scale_colour_manual(values = col_bases, name = "Base") +
  scale_x_continuous(breaks = breaks_pares, labels = minus) +
  facet_grid(. ~ rep) +
  th +
  labs(
    x = "Posición respecto del TSS (TSS = 0; punteadas en −1 y 0)",
    y = "Diferencia de actividad media\ncon vs. sin la base (HL, IC 95%)",
    caption = "Todos los promotores. Punto lleno: BH < 0,05 (dentro de réplica); vacío: n.s."
  )
ggsave(file.path(out_dir, "tss_efecto_bases.jpg"), p_efecto_bases, width = 12, height = 4.5, units = "in")

# --- ajuste por G+C local: coeficiente crudo vs. ajustado (todos los promotores)
p_gc <- perfil %>%
  filter(grupo == "Todos", conjunto == "Base") %>%
  select(rep, ventana, posicion, base, Crudo = lm_crudo_est, `Ajustado por G+C local` = lm_gc_est, lm_gc_p_adj) %>%
  pivot_longer(c(Crudo, `Ajustado por G+C local`), names_to = "modelo", values_to = "coef") %>%
  mutate(ventana = factor(ventana, c(ventana_tss, ventana_control)), modelo = factor(modelo, c("Crudo", "Ajustado por G+C local"))) %>%
  ggplot(aes(posicion, coef, color = base, linetype = modelo)) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_line() +
  geom_point(data = ~ filter(.x, modelo != "Crudo", lm_gc_p_adj < 0.05), size = 1.6) +
  scale_color_manual(values = col_bases, name = "Base") +
  scale_x_continuous(breaks = breaks_pares, labels = minus) +
  facet_grid(rep ~ ventana, scales = "free_x") +
  th +
  labs(x = "Posición respecto del TSS", y = "Coeficiente sobre el percentil de actividad\n(base vs. no base)", linetype = NULL,
       caption = "Todos los promotores. Puntos: ajustado por G+C local con BH < 0,05.")
ggsave(file.path(out_dir, "tss_ajuste_gc.jpg"), p_gc, width = 13, height = 8, units = "in")

# =============================================================================
# Puntos 3-5: grafico de puntos con IC comun (formato de tss_efecto_cgi.jpg)
# =============================================================================

forest <- function(res, y, color, valores_color, x_lab, facets = . ~ rep, dodge = 0.7) {
  # group explicito: sin el, geom_point se agrupa por color x forma (sig) y
  # el dodge de los puntos no coincide con el de las barras
  ggplot(res, aes(estimate, {{ y }}, color = {{ color }}, group = {{ color }})) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = dodge), linewidth = 0.9) +
    geom_point(aes(shape = sig), size = 2.6, fill = "white", position = position_dodge(width = dodge)) +
    scale_shape_manual(values = c(`BH < 0,05` = 16, `n.s.` = 21, `Suma (aditiva)` = 23), name = NULL) +
    scale_color_manual(values = valores_color, name = NULL) +
    facet_grid(facets) +
    th +
    labs(x = x_lab, y = NULL)
}
SIG_NIVELES <- c("BH < 0,05", "n.s.", "Suma (aditiva)")
sig_lab <- function(p_adj) factor(if_else(!is.na(p_adj) & p_adj < 0.05, "BH < 0,05", "n.s."), SIG_NIVELES)
x_lab_hl <- "Diferencia de actividad media con vs. sin la característica\n(Hodges-Lehmann, IC 95%)"
col_modelo <- c("grey55", "#216869") %>% set_names(modelos_lab[1:2])

# =============================================================================
# PUNTO 3 - iniciador humano (BBCA+1BW) vs. YCASW
# =============================================================================

INR_DEFS <- c(INR_humano = "Inr humano (BBCA+1BW)", INR_strong_TSS = "Inr fuerte (YCASW)")
col_inr <- c("#D6741F", "#216869") %>% set_names(INR_DEFS)

inr_prev <- prevalencia_feats(prev_universos, names(INR_DEFS))

# solapamiento de las dos definiciones, por universo y grupo de isla CpG
inr_solap <- prev_universos %>%
  mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")) %>%
  bind_rows(mutate(., grupo = "Todos")) %>%
  group_by(universo, grupo) %>%
  summarise(
    n = n(),
    n_ambos = sum(INR_humano & INR_strong_TSS), n_solo_humano = sum(INR_humano & !INR_strong_TSS),
    n_solo_ycasw = sum(!INR_humano & INR_strong_TSS), n_ninguno = sum(!INR_humano & !INR_strong_TSS),
    .groups = "drop"
  ) %>%
  mutate(
    pct_humano = 100 * (n_ambos + n_solo_humano) / n, pct_ycasw = 100 * (n_ambos + n_solo_ycasw) / n,
    jaccard = n_ambos / (n_ambos + n_solo_humano + n_solo_ycasw),
    pct_humano_que_es_ycasw = 100 * n_ambos / (n_ambos + n_solo_humano),
    pct_ycasw_que_es_humano = 100 * n_ambos / (n_ambos + n_solo_ycasw)
  )

inr_efecto <- efecto_feats(data_g, names(INR_DEFS))
inr_modelo <- modelos_feats(data, names(INR_DEFS), interaccion = FALSE)

save_tsv(
  bind_rows(
    inr_prev %>% mutate(analisis = "prevalencia_cgi"),
    inr_solap %>% mutate(analisis = "solapamiento"),
    inr_efecto %>% mutate(analisis = "efecto_HL"),
    inr_modelo %>% filter(term == "feat") %>% mutate(analisis = "modelo_lm")
  ) %>% relocate(analisis, caracteristica),
  "tss_inr_humano.tsv"
)

p_inr_hl <- inr_efecto %>%
  mutate(caracteristica = factor(INR_DEFS[caracteristica], rev(INR_DEFS)), sig = sig_lab(p_adj)) %>%
  forest(caracteristica, grupo, col_grupo, x_lab_hl)
p_inr_lm <- inr_modelo %>%
  filter(term == "feat") %>%
  mutate(caracteristica = factor(INR_DEFS[caracteristica], rev(INR_DEFS)), sig = sig_lab(p.value)) %>%
  forest(caracteristica, modelo, col_modelo, "Coeficiente sobre el percentil de actividad\n(modelo lineal, IC 95%; p sin corregir)", dodge = 0.6)
ggsave(file.path(out_dir, "tss_inr_humano.jpg"), p_inr_hl / p_inr_lm, width = 11, height = 8, units = "in")

# =============================================================================
# PUNTO 4 - interaccion TATA-box x iniciador
# =============================================================================

COMBOS <- c("Solo TATA", "Solo iniciador", "Ambos")
ADITIVA <- "Suma solo TATA + solo iniciador"
MIN_CELDA <- 50

tata_inr <- data_g %>%
  group_by(rep, grupo) %>%
  group_modify(function(d, k) {
    map_dfr(names(INR_DEFS), function(f) {
      tata <- d$TATA_EPD
      inr <- d[[f]]
      celda <- case_when(!tata & !inr ~ "Ninguno", tata & !inr ~ "Solo TATA", !tata & inr ~ "Solo iniciador", TRUE ~ "Ambos")
      ref <- d$mean[celda == "Ninguno"]
      combos <- map_dfr(COMBOS, ~ tata_hl_vs(d$mean[celda == .x], ref, MIN_N) %>%
        mutate(comparacion = .x, n = sum(celda == .x), n_ref = length(ref)))
      aditiva <- tibble(comparacion = ADITIVA, estimate = sum(combos$estimate[combos$comparacion != "Ambos"]))
      tata_cond <- bind_rows(
        tata_hl_vs(d$mean[tata & inr], d$mean[!tata & inr], MIN_N) %>%
          mutate(comparacion = "TATA con iniciador", n = sum(tata & inr), n_ref = sum(!tata & inr)),
        tata_hl_vs(d$mean[tata & !inr], d$mean[!tata & !inr], MIN_N) %>%
          mutate(comparacion = "TATA sin iniciador", n = sum(tata & !inr), n_ref = sum(!tata & !inr))
      )
      bind_rows(combos, aditiva, tata_cond) %>%
        mutate(
          p_adj = p.adjust(p, "BH"), iniciador = f,
          n_ninguno = length(ref), n_solo_tata = sum(celda == "Solo TATA"),
          n_solo_inr = sum(celda == "Solo iniciador"), n_ambos = sum(celda == "Ambos")
        )
    })
  }) %>%
  ungroup()

# interaccion en el modelo lineal (sin CGI_frag como covariable dentro de un grupo de isla CpG)
tata_inr_lm <- data_g %>%
  group_by(rep, grupo) %>%
  group_modify(function(d, k) {
    cov <- if (k$grupo == "Todos") "+ CGI_frag + CCAAT_EPD + GCbox_EPD" else "+ CCAAT_EPD + GCbox_EPD"
    map_dfr(names(INR_DEFS), function(f) {
      d2 <- d %>% mutate(inr = .data[[f]])
      broom::tidy(lm(as.formula(paste("act_pct ~ TATA_EPD * inr", cov)), data = d2), conf.int = TRUE) %>% mutate(iniciador = f)
    })
  }) %>%
  ungroup() %>%
  mutate(term = str_replace_all(term, "TRUE", ""))

celdas_chicas <- tata_inr %>%
  distinct(rep, grupo, iniciador, n_ninguno, n_solo_tata, n_solo_inr, n_ambos) %>%
  pivot_longer(starts_with("n_"), names_to = "celda", values_to = "n") %>%
  filter(n < MIN_CELDA)
if (nrow(celdas_chicas) > 0) {
  message("Celdas TATA x iniciador con menos de ", MIN_CELDA, " promotores:")
  print(celdas_chicas, n = Inf)
}

save_tsv(
  bind_rows(
    tata_inr %>% mutate(analisis = "HL_vs_ninguno_y_TATA_condicional"),
    tata_inr_lm %>% mutate(analisis = "modelo_lm")
  ) %>% relocate(analisis, rep, grupo, iniciador),
  "tss_tata_inr.tsv"
)

niveles_comp <- rev(c(COMBOS, ADITIVA, "TATA con iniciador", "TATA sin iniciador"))
p_tata_inr <- tata_inr %>%
  mutate(
    comparacion = factor(comparacion, niveles_comp),
    iniciador = factor(INR_DEFS[iniciador], INR_DEFS),
    sig = if_else(comparacion == ADITIVA, factor("Suma (aditiva)", SIG_NIVELES), sig_lab(p_adj))
  ) %>%
  forest(comparacion, iniciador, col_inr,
    "Diferencia de actividad media (Hodges-Lehmann, IC 95%): combinaciones vs. «ninguno»;\nTATA con / sin iniciador vs. sin TATA con el mismo estado de iniciador",
    facets = grupo ~ rep
  )
p_tata_inr_lm <- tata_inr_lm %>%
  filter(term == "TATA_EPD:inr") %>%
  mutate(iniciador = factor(INR_DEFS[iniciador], INR_DEFS), grupo = factor(grupo, rev(GRUPOS)), sig = sig_lab(p.value)) %>%
  forest(grupo, iniciador, col_inr,
    "Interacción TATA-box × iniciador sobre el percentil de actividad\n(modelo lineal, IC 95%; p sin corregir)",
    dodge = 0.6
  )
ggsave(file.path(out_dir, "tss_tata_inr.jpg"), p_tata_inr / p_tata_inr_lm + plot_layout(heights = c(3, 1)),
  width = 12, height = 13, units = "in"
)

# =============================================================================
# PUNTO 5 - TCT y genes de proteinas ribosomales citoplasmaticas
# =============================================================================
# RP = RPL<n>[A/X], RPS<n>[A/X], RPLP0-2 y RPSA. Quedan afuera los
# mitocondriales (MRPL/MRPS), los paralogos tipo "-L" (RPL3L, RPL39L,
# RPS27L, RPL22L1, RPL26L1, RPL7L1, RPL36AL) y otros genes con prefijo RP
# que no son proteinas ribosomales (RPS6K*, RPS19BP1).

RP_REGEX <- "^(RP[LS]\\d+[AX]?|RPLP[0-2]|RPSA)$"
RP_CANDIDATOS <- "^(M?RP[LS]|RPLP|RPSA)"
prom_df <- prom_df %>% mutate(es_rp = str_detect(gene_sym, RP_REGEX))
data <- data %>% left_join(prom_df %>% select(seq_id, es_rp), by = "seq_id")

rp_lista <- prom_df %>%
  filter(str_detect(gene_sym, RP_CANDIDATOS)) %>%
  mutate(
    incluido = es_rp,
    motivo_exclusion = case_when(
      es_rp ~ NA_character_,
      str_detect(gene_sym, "^MRP") ~ "mitocondrial",
      str_detect(gene_sym, "^RPS6K|BP1$") ~ "no es proteína ribosomal",
      TRUE ~ "parálogo tipo -L"
    ),
    activo_rep1 = seq_id %in% data$seq_id[data$rep == "Rep 1"],
    activo_rep2 = seq_id %in% data$seq_id[data$rep == "Rep 2"]
  ) %>%
  select(incluido, motivo_exclusion, gene_sym, seq_id, TSS_seq, TCT_TSS, activo_rep1, activo_rep2) %>%
  arrange(desc(incluido), motivo_exclusion, gene_sym, seq_id)
save_tsv(rp_lista, "tss_tct_rp_promotores.tsv")

rp_conteos <- bind_rows(
  prom_df %>% mutate(universo = "Library (promotores FP)"),
  data %>% mutate(universo = paste0("Con actividad (", rep, ")"))
) %>%
  group_by(universo) %>%
  summarise(
    n = n(), n_tct = sum(TCT_TSS), n_tct_rp = sum(TCT_TSS & es_rp),
    n_rp = sum(es_rp), n_genes_rp = n_distinct(gene_sym[es_rp]), n_rp_tct = sum(es_rp & TCT_TSS),
    n_genes_rp_con_tct = n_distinct(gene_sym[es_rp & TCT_TSS]),
    .groups = "drop"
  ) %>%
  mutate(pct_tct_que_son_rp = 100 * n_tct_rp / n_tct, pct_rp_con_tct = 100 * n_rp_tct / n_rp)

subconj_tct <- bind_rows(
  data %>% mutate(subconjunto = "Todos"),
  data %>% filter(!es_rp) %>% mutate(subconjunto = "Sin genes RP"),
  data %>% filter(es_rp) %>% mutate(subconjunto = "Solo genes RP")
)
tct_efecto <- efecto_feats(subconj_tct, "TCT_TSS", c("rep", "subconjunto"))
tct_modelo <- modelos_feats(subconj_tct, "TCT_TSS", c("rep", "subconjunto"), interaccion = FALSE)

rp_contexto <- data %>%
  group_by(rep) %>%
  group_modify(~ tata_hl_vs(.x$mean[.x$es_rp], .x$mean[!.x$es_rp], MIN_N) %>%
    mutate(
      n_rp = sum(.x$es_rp), n_resto = sum(!.x$es_rp),
      pct_medio_rp = mean(.x$act_pct[.x$es_rp]), pct_medio_resto = mean(.x$act_pct[!.x$es_rp]),
      pct_mediana_rp = median(.x$act_pct[.x$es_rp]), pct_mediana_resto = median(.x$act_pct[!.x$es_rp])
    )) %>%
  ungroup()

save_tsv(
  bind_rows(
    rp_conteos %>% mutate(analisis = "conteos"),
    tct_efecto %>% mutate(analisis = "efecto_HL_TCT"),
    tct_modelo %>% filter(term == "feat") %>% mutate(analisis = "modelo_lm_TCT"),
    rp_contexto %>% mutate(analisis = "actividad_RP_vs_resto")
  ) %>% relocate(analisis),
  "tss_tct_rp.tsv"
)

niveles_sub <- rev(c("Todos", "Sin genes RP", "Solo genes RP"))
col_sub <- c("Todos" = "grey30", "Sin genes RP" = "#14AFB2", "Solo genes RP" = "#D6741F")
p_tct_hl <- tct_efecto %>%
  mutate(subconjunto = factor(subconjunto, niveles_sub), sig = sig_lab(p_adj)) %>%
  forest(subconjunto, subconjunto, col_sub, "Diferencia de actividad media con vs. sin TCT\n(Hodges-Lehmann, IC 95%)") +
  geom_text(aes(label = paste0("n = ", n, " / ", n_ref)), vjust = -1.1, size = 3, color = "grey30") +
  guides(color = "none")
p_tct_lm <- tct_modelo %>%
  filter(term == "feat") %>%
  mutate(subconjunto = factor(subconjunto, niveles_sub), sig = sig_lab(p.value)) %>%
  forest(subconjunto, modelo, col_modelo, "Coeficiente del TCT sobre el percentil de actividad\n(modelo lineal, IC 95%; p sin corregir)", dodge = 0.6)
p_rp <- rp_contexto %>%
  mutate(y = "Genes RP vs. resto", sig = sig_lab(p)) %>%
  forest(y, y, c("Genes RP vs. resto" = "#D6741F"), "Diferencia de actividad media, promotores RP vs. resto\n(Hodges-Lehmann, IC 95%)") +
  geom_text(aes(label = sprintf("percentil medio %.0f vs. %.0f", pct_medio_rp, pct_medio_resto)), vjust = -1.3, size = 3.2, color = "grey30") +
  guides(color = "none")
ggsave(file.path(out_dir, "tss_tct_rp.jpg"), p_tct_hl / p_tct_lm / p_rp + plot_layout(heights = c(1.3, 1.3, 0.6)),
  width = 11, height = 11, units = "in"
)

# =============================================================================
# PUNTO 6 - perfil de efectos en todo el promotor y frecuencia vs. efecto
# =============================================================================
# Mismo modelo que tss_ajuste_gc (percentil ~ base y percentil ~ base +
# GC_local, GC_local = G+C de +-5 pb sin la propia posicion; en los bordes
# del fragmento la ventana queda truncada) en las 252 posiciones (-236 a
# +15), por replica x grupo de isla CpG. BH dentro de replica x grupo (las
# 252 x 4 pruebas del modelo ajustado). Por velocidad se usa lm.fit.
#
# Frecuencia vs. efecto: enriquecimiento posicional = log2(frecuencia de la
# base en la posicion / frecuencia media de esa base en las 252 posiciones),
# calculado sobre los mismos promotores con actividad de cada replica x
# grupo (no sobre la library completa de perfil_nucleotidos.tsv). Efecto =
# coeficiente ajustado por G+C local. Spearman y regresion ponderada por
# 1/SE^2, en todas las posiciones y en subconjuntos.

POS_TODAS <- seq_len(TATA_FRAG_LEN) - TSS_POS # -236..+15
ZOOM <- c(-50, 15)
TATA_BANDA <- TATA_CANON
th_perfil <- theme_pubclean(base_size = 14) +
  theme(strip.text = element_text(face = "bold"), legend.position = "right")
marca_tata <- annotate("rect", xmin = TATA_BANDA[1] - 0.5, xmax = TATA_BANDA[2] + 0.5, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.4)
marca_tss <- geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")
marca_cero <- geom_hline(yintercept = 0, colour = "grey60")
x_completo <- scale_x_continuous(breaks = c(-200, -150, -100, -50, 0), labels = minus, expand = expansion(0.01))
x_zoom <- scale_x_continuous(breaks = seq(-50, 10, 10), labels = minus, expand = expansion(0.02))

# coeficiente de la columna 2 de X (la base) con su SE y p
coef_base <- function(y, X) {
  fit <- lm.fit(X, y)
  if (fit$rank < ncol(X)) return(c(est = NA_real_, se = NA_real_, p = NA_real_))
  s2 <- sum(fit$residuals^2) / fit$df.residual
  se <- sqrt(diag(chol2inv(qr.R(fit$qr)))[2] * s2)
  est <- unname(fit$coefficients[2])
  c(est = est, se = se, p = 2 * pt(-abs(est / se), fit$df.residual))
}

efectos_promotor <- function(d) {
  M <- seq_mat[d$seq_id, , drop = FALSE]
  es_gc <- (M == "G" | M == "C") * 1
  y <- d$act_pct
  map_dfr(seq_len(ncol(M)), function(col) {
    vecinas <- setdiff(max(1, col - GC_FLANCO):min(ncol(M), col + GC_FLANCO), col)
    gc_local <- rowMeans(es_gc[, vecinas, drop = FALSE])
    b <- M[, col]
    ok <- b %in% BASES
    map_dfr(BASES, function(base) {
      ind <- (b[ok] == base) * 1
      out <- tibble(posicion = col - TSS_POS, base = base, n = sum(ind), n_ref = sum(1 - ind), frecuencia = mean(ind))
      if (sum(ind) < MIN_N || sum(1 - ind) < MIN_N) return(out)
      cr <- coef_base(y[ok], cbind(1, ind))
      aj <- coef_base(y[ok], cbind(1, ind, gc_local[ok]))
      out %>% mutate(
        crudo_est = cr[["est"]], crudo_se = cr[["se"]], crudo_p = cr[["p"]],
        estimate = aj[["est"]], se = aj[["se"]], p = aj[["p"]]
      )
    })
  })
}

efectos_prom <- data_g %>%
  group_by(rep, grupo) %>%
  group_modify(~ efectos_promotor(.x)) %>%
  ungroup() %>%
  group_by(rep, grupo) %>%
  mutate(
    conf.low = estimate - qnorm(0.975) * se, conf.high = estimate + qnorm(0.975) * se,
    p_adj = p.adjust(p, "BH"),
    enriquecimiento = log2(frecuencia / ave(frecuencia, base, FUN = mean))
  ) %>%
  ungroup() %>%
  relocate(rep, grupo, posicion, base)
save_tsv(efectos_prom, "tss_perfil_efectos_promotor.tsv")

# --- grafico: fila superior fragmento completo, inferior zoom (formato M12)
perfil_efectos_plot <- function(res, facets, y_lab = "Efecto ajustado por G+C local\n(Δ percentil de actividad)") {
  y_lims <- range(res$estimate, na.rm = TRUE)
  completo <- ggplot(res, aes(posicion, estimate, colour = base)) +
    marca_tata + marca_cero + marca_tss +
    geom_line(linewidth = 0.5) +
    facet_grid(facets) +
    scale_colour_manual(values = col_bases, name = "Base") +
    x_completo +
    coord_cartesian(ylim = y_lims) +
    labs(x = "Posición respecto del TSS (pb)", y = y_lab) +
    th_perfil
  zoom <- res %>%
    filter(between(posicion, ZOOM[1], ZOOM[2])) %>%
    ggplot(aes(posicion, estimate, colour = base)) +
    marca_tata + marca_cero + marca_tss +
    geom_line(linewidth = 0.6) +
    geom_point(aes(shape = p_adj < 0.05), size = 1.4, fill = "white", show.legend = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21)) +
    facet_grid(facets) +
    scale_colour_manual(values = col_bases, name = "Base") +
    x_zoom +
    coord_cartesian(ylim = y_lims) +
    labs(x = paste0("Posición respecto del TSS (pb), ", minus(ZOOM[1]), " a +", ZOOM[2]),
         y = y_lab,
         caption = "Banda gris: −31 a −26 (TATA canónica). Abajo, punto lleno: BH < 0,05 (dentro de réplica × grupo).") +
    th_perfil
  completo / zoom + plot_layout(guides = "collect")
}
ggsave(file.path(out_dir, "tss_perfil_efectos_promotor.jpg"),
  perfil_efectos_plot(filter(efectos_prom, grupo == "Todos"), . ~ rep), width = 14, height = 8, units = "in")
ggsave(file.path(out_dir, "tss_perfil_efectos_promotor_cgi.jpg"),
  perfil_efectos_plot(efectos_prom, grupo ~ rep), width = 14, height = 16, units = "in")

# --- version HL (sin ajuste por G+C) del perfil de efectos -------------------
# HL (IC 95%) de la actividad media con vs. sin la base, en cada posicion de
# -236 a +15, por replica x grupo de isla CpG; BH dentro de replica x grupo.
# Son ~6000 Wilcoxon con IC: se reparten en un cluster local y el resultado
# crudo (sin BH) se cachea en data/processed/ (borrar el .rds para recalcular).
cache_hl <- "data/processed/tss_composicion_perfil_hl_v1.rds"
if (file.exists(cache_hl)) {
  efectos_hl <- readRDS(cache_hl)
} else {
  t0 <- Sys.time()
  d_split_hl <- data_g %>% select(rep, grupo, seq_id, mean) %>% group_split(rep, grupo)
  tareas_hl <- expand_grid(g = seq_along(d_split_hl), col = seq_len(TATA_FRAG_LEN))
  hl_tarea <- function(i) {
    d <- d_split_hl[[tareas_hl$g[i]]]
    col <- tareas_hl$col[i]
    b <- seq_mat[d$seq_id, col]
    ok <- b %in% BASES
    map_dfr(BASES, function(base) {
      ind <- b[ok] == base
      tata_hl_vs(d$mean[ok][ind], d$mean[ok][!ind], MIN_N) %>%
        mutate(rep = d$rep[1], grupo = d$grupo[1], posicion = col - TSS_POS, base = base,
               n = sum(ind), n_ref = sum(!ind), frecuencia = mean(ind))
    })
  }
  cl <- parallel::makeCluster(max(1, min(12, parallel::detectCores() - 2)))
  parallel::clusterEvalQ(cl, suppressMessages(library(tidyverse)))
  parallel::clusterExport(cl, c("d_split_hl", "tareas_hl", "seq_mat", "tata_hl_vs", "TSS_POS", "BASES", "MIN_N"))
  efectos_hl <- bind_rows(parallel::parLapplyLB(cl, seq_len(nrow(tareas_hl)), hl_tarea))
  parallel::stopCluster(cl)
  saveRDS(efectos_hl, cache_hl)
  message("Perfil HL: ", format(Sys.time() - t0))
}
efectos_hl <- efectos_hl %>%
  group_by(rep, grupo) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup() %>%
  relocate(rep, grupo, posicion, base, n, n_ref, frecuencia)
write_tsv(efectos_hl, file.path(dir_perfil_efectos, "tss_perfil_efectos_promotor_hl.tsv"))

y_lab_hl <- "Diferencia de actividad media\ncon vs. sin la base (HL)"
ggsave(file.path(dir_perfil_efectos, "tss_perfil_efectos_promotor_hl.jpg"),
  perfil_efectos_plot(filter(efectos_hl, grupo == "Todos"), . ~ rep, y_lab_hl), width = 14, height = 8, units = "in")
ggsave(file.path(out_dir, "tss_perfil_efectos_promotor_hl_cgi.jpg"),
  perfil_efectos_plot(efectos_hl, grupo ~ rep, y_lab_hl), width = 14, height = 16, units = "in")

# --- frecuencia vs. efecto -------------------------------------------------

SUBCONJ <- c(
  "a) Todas (−236 a +15)", "b) −50 a +15", "c) −236 a −51", "d) Sin −1, 0 ni −31 a −26"
)
en_subconj <- function(pos, s) {
  switch(s,
    "a) Todas (−236 a +15)" = rep(TRUE, length(pos)),
    "b) −50 a +15" = between(pos, -50, 15),
    "c) −236 a −51" = pos <= -51,
    "d) Sin −1, 0 ni −31 a −26" = !(pos %in% c(-1, 0) | between(pos, TATA_BANDA[1], TATA_BANDA[2]))
  )
}
cor_frec_efecto <- function(d) {
  d <- filter(d, !is.na(estimate), is.finite(enriquecimiento))
  if (nrow(d) < 10) return(NULL)
  sp <- cor.test(d$enriquecimiento, d$estimate, method = "spearman", exact = FALSE)
  w <- broom::tidy(lm(estimate ~ enriquecimiento, data = d, weights = 1 / se^2), conf.int = TRUE) %>%
    filter(term == "enriquecimiento")
  tibble(
    n_puntos = nrow(d), rho = unname(sp$estimate), p_spearman = sp$p.value,
    pendiente_pond = w$estimate, pendiente_conf.low = w$conf.low, pendiente_conf.high = w$conf.high, p_pond = w$p.value
  )
}
frec_efecto <- expand_grid(subconjunto = SUBCONJ, bases = c("Todas", BASES)) %>%
  pmap_dfr(function(subconjunto, bases) {
    efectos_prom %>%
      filter(en_subconj(posicion, subconjunto), bases == "Todas" | base == bases) %>%
      group_by(rep, grupo) %>%
      group_modify(~ cor_frec_efecto(.x)) %>%
      ungroup() %>%
      mutate(subconjunto = subconjunto, bases = bases)
  }) %>%
  relocate(rep, grupo, subconjunto, bases)
save_tsv(frec_efecto, "tss_frecuencia_vs_efecto.tsv")

frec_efecto_plot <- function(res, facets) {
  res <- res %>%
    filter(!is.na(estimate)) %>%
    mutate(
      zona = case_when(
        posicion %in% c(-1, 0) ~ "−1 / 0",
        between(posicion, TATA_BANDA[1], TATA_BANDA[2]) ~ "TATA (−31 a −26)",
        TRUE ~ "Resto"
      ),
      etiqueta = if_else(zona == "−1 / 0", paste0(base, minus(posicion)), NA_character_)
    )
  rotulos <- frec_efecto %>%
    filter(bases == "Todas", subconjunto %in% SUBCONJ[c(1, 4)]) %>%
    semi_join(distinct(res, rep, grupo), by = c("rep", "grupo")) %>%
    group_by(rep, grupo) %>%
    summarise(label = paste0("ρ todas = ", sprintf("%.2f", rho[subconjunto == SUBCONJ[1]]),
                             "\nρ sin −1/0/TATA = ", sprintf("%.2f", rho[subconjunto == SUBCONJ[4]])), .groups = "drop")
  ggplot(res, aes(enriquecimiento, estimate)) +
    marca_cero +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(data = ~ filter(.x, zona == "Resto"), aes(colour = base), size = 1.1, alpha = 0.45) +
    geom_point(data = ~ filter(.x, zona != "Resto"), aes(fill = base, shape = zona), size = 3, colour = "black", stroke = 0.6) +
    ggrepel::geom_text_repel(aes(label = etiqueta, colour = base), size = 3.5, na.rm = TRUE, show.legend = FALSE, min.segment.length = 0) +
    geom_text(data = rotulos, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2, size = 3.5, inherit.aes = FALSE) +
    scale_colour_manual(values = col_bases, name = "Base") +
    scale_fill_manual(values = col_bases, guide = "none") +
    scale_shape_manual(values = c("−1 / 0" = 21, "TATA (−31 a −26)" = 24), name = NULL) +
    facet_grid(facets) +
    th_perfil +
    labs(
      x = "Enriquecimiento posicional\nlog2(frecuencia en la posición / frecuencia media de la base en el fragmento)",
      y = "Efecto ajustado por G+C local\n(Δ percentil de actividad)",
      caption = "Un punto por posición × base (−236 a +15). ρ de Spearman con todas las bases juntas."
    )
}
ggsave(file.path(out_dir, "tss_frecuencia_vs_efecto.jpg"),
  frec_efecto_plot(filter(efectos_prom, grupo == "Todos"), . ~ rep), width = 13, height = 6, units = "in")
ggsave(file.path(out_dir, "tss_frecuencia_vs_efecto_cgi.jpg"),
  frec_efecto_plot(efectos_prom, grupo ~ rep), width = 13, height = 15, units = "in")
