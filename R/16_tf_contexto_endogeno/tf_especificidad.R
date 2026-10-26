# EXPLORATORIO - especificidad de NF-Y y Sp frente al resto de los
# factores de ReMap2022. Rehace las comparaciones motivo vs. TF de
# motivo_vs_tf.R y motivo_tf_estratificado.R / sp2_estratificado_ccaat.R
# (que usan la union de NFYA en cualquier linea celular, ~64% de los
# promotores) con definiciones mas estrictas:
#   - NFYA: interseccion HeLa-S3 / Hep-G2 / K-562 (excluye GM12878, igual
#     que NFYA_scatter_intersect en tf_contexto_endogeno.R);
#   - SP1/SP2: union solo en muestras HEK293 (str_detect "HEK293":
#     HEK293, HEK293T, HEK293-FT, igual que tf_contexto_endogeno.R).
# y las ubica en la distribucion de efectos de todos los TFs de ReMap:
#   (1) motivo vs. TF y estratificados por motivo  -> *_intersect / *_hek
#   (2a) percentil del efecto de NFYA/SP1/SP2 entre todos los TFs con
#        >100 promotores unidos (cualquier linea: reutiliza el ranking de
#        remap_tf_activity_wilcoxon.tsv de R/17_remap_activity_gsea;
#        HEK293: calculado aca)
#   (2b) control de cobertura: TFs con n de promotores unidos +-25%, y
#        residuo de loess efecto ~ log10(n unidos)
#   (2c) especificidad motivo x factor: efecto(con motivo) - efecto(sin
#        motivo) para cada TF; NFYA x CCAAT-box, SP1/SP2 x GC-box
#   (2d) efecto de SP1/SP2 (HEK) y NFYA (interseccion) ajustado por
#        actividad endogena en HEK293 (hek_tpm) y DNase (mean_dnase),
#        lm sobre el percentil de actividad
#
# Efecto = estimador de Hodges-Lehmann (unido - no unido) de la
# actividad media del reportero, por replica. En (1) y para NFYA
# interseccion en (2a) se usa coin::wilcox_test como en motivo_vs_tf.R /
# R/17; para los ~1000 TFs x estratos de (2a HEK)-(2c) se usa hl_test()
# (abajo), mismo estimador calculado por biseccion sobre las diferencias
# ordenadas, con IC por estadisticos de orden (aproximacion normal) - se
# valida contra la tabla de R/17 al correr.
#
# Salidas (todas en figures/exploratorio_tf_contexto_endogeno/, no pisa las
# tablas viejas): motivo_vs_tf_efecto_intersect.tsv/.jpg,
# nfya_estratificado_ccaat_intersect.tsv/.jpg,
# sp1_estratificado_gcbox_hek.tsv/.jpg, sp2_estratificado_ccaat_hek.tsv/.jpg,
# tf_efecto_todos_intersect.tsv (2a, todos los TFs),
# tf_percentil_intersect.tsv (2a resumen), tf_cobertura_intersect.tsv (2b),
# tf_motivo_especificidad_intersect.tsv (2c, todos los TFs),
# tf_motivo_especificidad_percentil_intersect.tsv (2c resumen),
# tf_ajuste_endogeno_hek.tsv (2d), tf_efecto_distribucion_intersect.jpg,
# tf_motivo_especificidad_intersect.jpg.
#
# Requiere: remap_tf_hits.tsv, activity_stats_highconf.tsv, prom_df.tsv y
# la tabla remap_tf_activity_wilcoxon.tsv de R/17_remap_activity_gsea.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(coin)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

out_dir <- "figures/exploratorio_tf_contexto_endogeno" # exploratorio (no es figura de la tesis)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))
reps <- sort(unique(data$rep))

min_bound <- 100 # TFs con >100 promotores unidos en ambas replicas (como R/17)
min_strat <- 20 # minimo de unidos y de no unidos por estrato del motivo (2c)
cov_window <- 0.25 # control de cobertura (2b): n unidos +-25%

hits <- remap_hits %>% filter(seq_id %in% data$seq_id)
hits_any <- hits %>% distinct(TF, seq_id)
hits_hek <- hits %>% filter(str_detect(sample, "HEK293")) %>% distinct(TF, seq_id)

# NFYA: interseccion de lineas celulares (mismo criterio que
# tf_contexto_endogeno.R, GM12878 excluida)
nfya_lines <- remap_hits %>%
  filter(TF == "NFYA") %>%
  separate_rows(sample, sep = ",") %>%
  filter(sample != "GM12878") %>%
  group_by(sample) %>%
  summarise(fps = list(unique(seq_id)), .groups = "drop") %>%
  deframe()
stopifnot(setequal(names(nfya_lines), c("HeLa-S3", "Hep-G2", "K-562")))
nfya_intersect <- Reduce(intersect, nfya_lines)
message("NFYA interseccion: ", length(nfya_intersect), " promotores (", paste(names(nfya_lines), collapse = ", "), ")")

sets_any <- split(hits_any$seq_id, hits_any$TF)
sets_hek <- split(hits_hek$seq_id, hits_hek$TF)
sets_any[["NFYA_intersect"]] <- nfya_intersect

data <- data %>%
  mutate(
    NFYA_intersect = seq_id %in% nfya_intersect,
    SP1_hek = seq_id %in% sets_hek[["SP1"]],
    SP2_hek = seq_id %in% sets_hek[["SP2"]]
  )

# --- Hodges-Lehmann rapido -----------------------------------------------
# k-esima diferencia mas chica y_i - x_j (xs ordenado) por biseccion:
# #(y_i - x_j <= d) = sum_i #(x_j >= y_i - d).
kth_diff <- function(y, xs, k) {
  nx <- length(xs)
  lo <- min(y) - max(xs)
  hi <- max(y) - min(xs)
  for (i in 1:60) {
    mid <- (lo + hi) / 2
    if (sum(nx - findInterval(y - mid, xs, left.open = TRUE)) >= k) hi <- mid else lo <- mid
  }
  hi
}

# efecto (HL, unido - no unido), IC 95% y p de Wilcoxon (aprox. normal con
# correccion por empates)
hl_test <- function(v, g, conf = 0.95) {
  y <- v[g]
  xs <- sort(v[!g])
  ny <- length(y)
  nx <- length(xs)
  N <- ny * nx
  est <- if (N %% 2 == 1) kth_diff(y, xs, (N + 1) / 2) else (kth_diff(y, xs, N / 2) + kth_diff(y, xs, N / 2 + 1)) / 2
  n <- ny + nx
  r <- rank(v)
  t <- table(r)
  sigma <- sqrt(nx * ny / 12 * ((n + 1) - sum(t^3 - t) / (n * (n - 1))))
  W <- sum(r[g]) - ny * (ny + 1) / 2
  ca <- max(1, floor(N / 2 - qnorm(1 - (1 - conf) / 2) * sigma))
  tibble(
    estimate = est, P2.5 = kth_diff(y, xs, ca), P97.5 = kth_diff(y, xs, N - ca + 1),
    p = 2 * pnorm(-abs((W - N / 2) / sigma))
  )
}

# efectos de una lista de sets de seq_id sobre la actividad, por replica
# (y opcionalmente por estrato de un motivo)
tf_effects <- function(sets, tfs, motif = NULL) {
  d <- data %>% select(seq_id, rep, mean, any_of(motif))
  if (is.null(motif)) d$estrato <- "todos" else d$estrato <- ifelse(d[[motif]], "con", "sin")
  d %>%
    group_by(rep, estrato) %>%
    group_modify(~ {
      dd <- .x
      map(tfs, function(tf) {
        g <- dd$seq_id %in% sets[[tf]]
        if (sum(g) < min_strat || sum(!g) < min_strat) return(NULL)
        hl_test(dd$mean, g) %>% mutate(TF = tf, n_unidos = sum(g), n = length(g), .before = 1)
      }) %>% list_rbind()
    }) %>%
    ungroup()
}

n_bound <- function(sets, tfs) {
  map(reps, function(r) {
    ids <- data$seq_id[data$rep == r]
    tibble(TF = tfs, rep = r, n_unidos = map_int(tfs, ~ sum(ids %in% sets[[.x]])))
  }) %>% list_rbind()
}

eligible <- function(sets) {
  n_bound(sets, names(sets)) %>%
    group_by(TF) %>%
    filter(all(n_unidos > min_bound), n() == length(reps)) %>%
    ungroup() %>%
    filter(!str_detect(TF, "^H\\d")) %>% # marcas de histona, como R/17
    distinct(TF) %>%
    pull(TF)
}

# percentil (% de TFs de referencia con efecto <= x) y rango (1 = mayor)
pct_of <- function(x, ref) 100 * mean(ref <= x)
rank_of <- function(x, ref) sum(ref > x) + 1L

# =========================================================================
# (1) Motivo vs. TF con las definiciones estrictas
# =========================================================================

features <- tibble(
  rep = data$rep, mean = data$mean,
  `CCAAT-box` = data$CCAAT_EPD,
  `NFYA (ReMap, intersección)` = data$NFYA_intersect,
  `GC-box` = data$GCbox_EPD,
  `SP1 (ReMap, HEK293)` = data$SP1_hek,
  `SP2 (ReMap, HEK293)` = data$SP2_hek
)
feature_levels <- setdiff(names(features), c("rep", "mean"))

n_feat <- features %>%
  group_by(rep) %>%
  summarise(n = n(), across(all_of(feature_levels), sum)) %>%
  pivot_longer(all_of(feature_levels), names_to = "feature", values_to = "n_presente") %>%
  mutate(pct_presente = 100 * n_presente / n)
print(n_feat)

tidy_data <- features %>%
  mutate(rep = as.factor(rep), across(all_of(feature_levels), ~ factor(.x, levels = c("TRUE", "FALSE"))))
motif_tf_df <- wilcox_effect_summary(tidy_data) %>%
  mutate(feature = as.character(feature), rep = as.character(rep)) %>%
  left_join(n_feat %>% mutate(rep = as.character(rep)), by = c("feature", "rep"))
write_tsv(motif_tf_df, file.path(out_dir, "motivo_vs_tf_efecto_intersect.tsv"))

p <- motif_tf_df %>%
  select(-pval, -n, -n_presente, -pct_presente) %>%
  pivot_wider(names_from = val, values_from = estimate) %>%
  mutate(
    feature = factor(as.character(feature), levels = rev(feature_levels)),
    grupo = ifelse(str_detect(feature, "CCAAT|NFYA"), "CCAAT / NFYA", "GC-box / SP1 / SP2")
  ) %>%
  ggplot(aes(x = estimate, y = feature, group = fct_rev(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = grupo, alpha = rep)) +
  geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(0.9), height = 0.15, col = "#777777", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.6)) +
  scale_fill_manual(values = c("CCAAT / NFYA" = "#D6741F", "GC-box / SP1 / SP2" = "#216869")) +
  labs(x = "Efecto sobre la actividad (Wilcoxon)", y = NULL, fill = NULL, alpha = NULL) +
  theme_pubr(base_size = 18) +
  theme(legend.position = "top")
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("motivo_vs_tf_efecto_intersect.", ext)), p, width = 10, height = 5, units = "in")
}

# estratificados (mismo formato que motivo_tf_estratificado.R, + IC y %)
stratified_tf <- function(tf_col, tf_label, motif_col, motif_label, file_stem) {
  d <- data %>%
    mutate(
      TF = factor(ifelse(.data[[tf_col]], "Presencia", "Ausencia"), levels = c("Ausencia", "Presencia")),
      motivo = factor(ifelse(.data[[motif_col]], paste("Con", motif_label), paste("Sin", motif_label)),
        levels = paste(c("Sin", "Con"), motif_label)
      )
    )
  effect <- d %>%
    group_by(rep, motivo) %>%
    group_modify(~ {
      w <- wilcox.test(mean ~ TF, data = .x, conf.int = TRUE)
      tibble(
        n_tf = sum(.x$TF == "Presencia"), n = nrow(.x), pct_unidos = 100 * n_tf / n,
        estimate = -unname(w$estimate), P2.5 = -w$conf.int[2], P97.5 = -w$conf.int[1], p = w$p.value
      )
    }) %>%
    ungroup()
  print(effect)
  write_tsv(effect, file.path(out_dir, paste0(file_stem, ".tsv")))

  lab <- effect %>% mutate(label = sprintf("efecto = %.2f [%.2f; %.2f]\nn(%s) = %d / %d (%.0f%%)", estimate, P2.5, P97.5, tf_label, n_tf, n, pct_unidos))
  p <- ggplot(d, aes(TF, mean, fill = TF)) +
    geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
    stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), geom = "pointrange", size = 0.4) +
    geom_text(data = lab, aes(x = 1.5, y = 6.6, label = label), inherit.aes = FALSE, size = 3.6) +
    facet_grid(rep ~ motivo) +
    scale_fill_manual(values = c("Ausencia" = "#14AFB2", "Presencia" = "#216869")) +
    coord_cartesian(ylim = c(0.5, 7.2)) +
    labs(x = paste0("Unión de ", tf_label, " (ReMap)"), y = "Actividad media") +
    theme_pubclean(base_size = 16) +
    theme(legend.position = "none")
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(file_stem, ".", ext)), p, width = 9, height = 7, units = "in")
  }
  effect
}

strat_1 <- list(
  stratified_tf("NFYA_intersect", "NFYA ∩", "CCAAT_EPD", "CCAAT-box", "nfya_estratificado_ccaat_intersect"),
  stratified_tf("SP1_hek", "SP1 HEK293", "GCbox_EPD", "GC-box", "sp1_estratificado_gcbox_hek"),
  stratified_tf("SP2_hek", "SP2 HEK293", "CCAAT_EPD", "CCAAT-box", "sp2_estratificado_ccaat_hek")
)

# =========================================================================
# (2a) Efecto de todos los TFs: percentil de NFYA / SP1 / SP2
# =========================================================================

# Cualquier linea: reutiliza el ranking de R/17 (coin, mismo estimador)
r17 <- read_tsv(file.path(fig_dir("remap_tf_activity_gsea", create = FALSE), "remap_tf_activity_wilcoxon.tsv"), show_col_types = FALSE) %>%
  select(TF = feature, rep, val, estimate) %>%
  pivot_wider(names_from = val, values_from = estimate)
nfya_int_coin <- motif_tf_df %>%
  filter(feature == "NFYA (ReMap, intersección)") %>%
  select(rep, val, estimate) %>%
  pivot_wider(names_from = val, values_from = estimate) %>%
  mutate(TF = "NFYA_intersect")
tfs_any <- unique(r17$TF)

eff_any <- bind_rows(r17, nfya_int_coin) %>%
  left_join(n_bound(sets_any, c(tfs_any, "NFYA_intersect")), by = c("TF", "rep")) %>%
  mutate(universo = "Cualquier línea", fuente = ifelse(TF == "NFYA_intersect", "coin (este script)", "R/17 remap_tf_activity_wilcoxon.tsv"))

# validacion de hl_test() contra la tabla de R/17 (NFYA/SP1/SP2 + 20 al azar)
set.seed(1)
check_tfs <- c("NFYA", "SP1", "SP2", sample(setdiff(tfs_any, c("NFYA", "SP1", "SP2")), 20))
check <- tf_effects(sets_any, check_tfs) %>%
  select(TF, rep, est_fast = estimate, lo_fast = P2.5, hi_fast = P97.5) %>%
  inner_join(r17, by = c("TF", "rep"))
message(sprintf(
  "Validacion hl_test vs. R/17 (%d TF x rep): max |dif| estimate = %.2g, IC = %.2g",
  nrow(check), max(abs(check$est_fast - check$estimate)),
  max(abs(c(check$lo_fast - check$P2.5, check$hi_fast - check$P97.5)))
))

# HEK293: TFs con datos en HEK y >100 promotores unidos (solo HEK) en ambas replicas
tfs_hek <- eligible(sets_hek)
message("TFs elegibles: cualquier linea = ", length(tfs_any), ", HEK293 = ", length(tfs_hek))
eff_hek <- tf_effects(sets_hek, tfs_hek) %>%
  select(TF, rep, n_unidos, estimate, P2.5, P97.5, p) %>%
  mutate(universo = "HEK293", fuente = "hl_test (este script)")

eff_all <- bind_rows(eff_any, eff_hek) %>%
  mutate(
    TF_focal = case_when(
      universo == "Cualquier línea" & TF == "NFYA_intersect" ~ "NFYA ∩",
      universo == "Cualquier línea" & TF %in% c("NFYA", "SP1", "SP2") ~ paste(TF, "(cualquier línea)"),
      universo == "HEK293" & TF %in% c("NFYA", "SP1", "SP2") ~ paste(TF, "HEK293"),
      TRUE ~ NA_character_
    )
  ) %>%
  select(universo, rep, TF, TF_focal, n_unidos, estimate, P2.5, P97.5, p, fuente)
write_tsv(eff_all, file.path(out_dir, "tf_efecto_todos_intersect.tsv"))

# referencia = todos los TFs del universo salvo el propio focal (y, para
# NFYA ∩, sin la union de NFYA: es el mismo factor)
same_factor <- function(tf) if (tf %in% c("NFYA", "NFYA_intersect")) c("NFYA", "NFYA_intersect") else tf

pct_tbl <- eff_all %>%
  filter(!is.na(TF_focal)) %>%
  pmap(function(universo, rep, TF, TF_focal, n_unidos, estimate, ...) {
    ref <- eff_all %>% filter(universo == !!universo, rep == !!rep, !TF %in% same_factor(!!TF))
    tibble(
      universo, rep, TF, TF_focal, n_unidos, estimate,
      n_ref = nrow(ref), percentil = pct_of(estimate, ref$estimate), rango = rank_of(estimate, ref$estimate),
      mediana_ref = median(ref$estimate), p90_ref = quantile(ref$estimate, 0.9, names = FALSE)
    )
  }) %>%
  list_rbind()
print(pct_tbl)
write_tsv(pct_tbl, file.path(out_dir, "tf_percentil_intersect.tsv"))

# =========================================================================
# (2b) Control de cobertura
# =========================================================================

cov_tbl <- pct_tbl %>%
  pmap(function(universo, rep, TF, TF_focal, n_unidos, estimate, ...) {
    ref <- eff_all %>% filter(universo == !!universo, rep == !!rep, !TF %in% same_factor(!!TF))
    win <- ref %>% filter(abs(n_unidos - !!n_unidos) <= cov_window * !!n_unidos)
    fit <- loess(estimate ~ log10(n_unidos), data = ref, span = 0.5)
    pred <- predict(fit, tibble(n_unidos = n_unidos))
    if (is.na(pred)) pred <- predict(lm(estimate ~ log10(n_unidos), data = ref), tibble(n_unidos = n_unidos))
    res_ref <- residuals(fit)
    tibble(
      universo, rep, TF, TF_focal, n_unidos, estimate,
      n_ventana = nrow(win), rango_n_ventana = paste(min(win$n_unidos), max(win$n_unidos), sep = "-"),
      percentil_ventana = if (nrow(win) > 0) pct_of(estimate, win$estimate) else NA_real_,
      mediana_ventana = if (nrow(win) > 0) median(win$estimate) else NA_real_,
      max_ventana = if (nrow(win) > 0) max(win$estimate) else NA_real_,
      TF_ventana_top5 = win %>% slice_max(estimate, n = 5) %>% pull(TF) %>% paste(collapse = ","),
      esperado_loess = unname(pred), residuo_loess = estimate - unname(pred),
      percentil_residuo = pct_of(estimate - unname(pred), res_ref)
    )
  }) %>%
  list_rbind()
print(cov_tbl %>% select(-TF_ventana_top5))
write_tsv(cov_tbl, file.path(out_dir, "tf_cobertura_intersect.tsv"))

# =========================================================================
# (2c) Especificidad motivo x factor
# =========================================================================

spec <- list(
  list(universo = "Cualquier línea", sets = sets_any, tfs = c(tfs_any, "NFYA_intersect")),
  list(universo = "HEK293", sets = sets_hek, tfs = tfs_hek)
) %>%
  map(function(u) {
    map(c("CCAAT_EPD", "GCbox_EPD"), function(m) {
      tf_effects(u$sets, u$tfs, motif = m) %>%
        select(rep, estrato, TF, n_unidos, n, estimate) %>%
        pivot_wider(names_from = estrato, values_from = c(n_unidos, n, estimate)) %>%
        filter(!is.na(estimate_con), !is.na(estimate_sin)) %>%
        mutate(universo = u$universo, motivo = m, dif_con_sin = estimate_con - estimate_sin, .before = 1)
    }) %>% list_rbind()
  }) %>%
  list_rbind() %>%
  mutate(motivo = recode(motivo, CCAAT_EPD = "CCAAT-box", GCbox_EPD = "GC-box"))
write_tsv(spec, file.path(out_dir, "tf_motivo_especificidad_intersect.tsv"))

spec_focal <- tribble(
  ~universo, ~motivo, ~TF, ~comparacion,
  "Cualquier línea", "CCAAT-box", "NFYA_intersect", "principal",
  "Cualquier línea", "CCAAT-box", "NFYA", "union laxa (referencia)",
  "HEK293", "GC-box", "SP1", "principal",
  "HEK293", "GC-box", "SP2", "principal",
  "HEK293", "CCAAT-box", "SP2", "control (motivo no propio)",
  "HEK293", "CCAAT-box", "SP1", "control (motivo no propio)",
  "Cualquier línea", "GC-box", "NFYA_intersect", "control (motivo no propio)"
)
spec_pct <- spec_focal %>%
  inner_join(spec, by = c("universo", "motivo", "TF")) %>%
  pmap(function(universo, motivo, TF, comparacion, rep, dif_con_sin, estimate_con, estimate_sin, n_unidos_con, n_unidos_sin, ...) {
    ref <- spec %>% filter(universo == !!universo, motivo == !!motivo, rep == !!rep, !TF %in% same_factor(!!TF))
    tibble(
      universo, motivo, TF, comparacion, rep, n_unidos_con, n_unidos_sin, estimate_con, estimate_sin, dif_con_sin,
      n_ref = nrow(ref), percentil = pct_of(dif_con_sin, ref$dif_con_sin), rango = rank_of(dif_con_sin, ref$dif_con_sin),
      mediana_ref = median(ref$dif_con_sin), p90_ref = quantile(ref$dif_con_sin, 0.9, names = FALSE)
    )
  }) %>%
  list_rbind()
print(spec_pct)
write_tsv(spec_pct, file.path(out_dir, "tf_motivo_especificidad_percentil_intersect.tsv"))

# =========================================================================
# (2d) Ajuste por actividad endogena (HEK293) y DNase
# =========================================================================

adj <- data %>%
  group_by(rep) %>%
  mutate(pct_act = 100 * percent_rank(mean)) %>%
  ungroup() %>%
  mutate(log_hek_tpm = log10(hek_tpm + 1)) %>%
  filter(!is.na(log_hek_tpm), !is.na(mean_dnase))
message("2d: ", n_distinct(adj$seq_id), " promotores con hek_tpm y mean_dnase")

covars <- list("sin ajuste" = NULL, "hek_tpm" = "log_hek_tpm", "mean_dnase" = "mean_dnase", "hek_tpm + mean_dnase" = c("log_hek_tpm", "mean_dnase"))
adj_tbl <- expand_grid(tf = c("SP1_hek", "SP2_hek", "NFYA_intersect"), rep = reps, ajuste = names(covars)) %>%
  pmap(function(tf, rep, ajuste) {
    d <- adj %>% filter(rep == !!rep)
    fit <- lm(reformulate(c(tf, covars[[ajuste]]), "pct_act"), data = d)
    term <- paste0(tf, "TRUE")
    ci <- confint(fit)[term, ]
    tibble(
      TF = tf, rep, ajuste, n = nobs(fit),
      coef = coef(fit)[[term]], ci_low = ci[[1]], ci_high = ci[[2]],
      p = summary(fit)$coefficients[term, 4], r2 = summary(fit)$r.squared
    )
  }) %>%
  list_rbind()
print(adj_tbl)
write_tsv(adj_tbl, file.path(out_dir, "tf_ajuste_endogeno_hek.tsv"))

# =========================================================================
# Figuras
# =========================================================================

focal_clr <- c(
  "NFYA ∩" = "#D6741F", "NFYA (cualquier línea)" = "#F2B680", "NFYA (laxa)" = "#F2B680",
  "SP1 (cualquier línea)" = "#8FC1C2", "SP2 (cualquier línea)" = "#8FC1C2",
  "SP1 HEK293" = "#216869", "SP2 HEK293" = "#14AFB2", "NFYA HEK293" = "#F2B680"
)

# etiquetas en una linea, escalonadas verticalmente dentro de cada panel
lab_2a <- pct_tbl %>%
  arrange(universo, rep, desc(estimate)) %>%
  group_by(universo, rep) %>%
  mutate(label = sprintf("%s (p%.0f)", TF_focal, percentil), vj = 1.5 + 1.5 * (row_number() - 1)) %>%
  ungroup()
p2a <- eff_all %>%
  filter(is.na(TF_focal)) %>%
  ggplot(aes(estimate)) +
  geom_histogram(bins = 60, fill = "grey75", color = NA) +
  geom_vline(data = lab_2a, aes(xintercept = estimate, color = TF_focal), linewidth = 1) +
  geom_text(data = lab_2a, aes(x = estimate, y = Inf, label = label, color = TF_focal, vjust = vj), hjust = 1.03, size = 3.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_grid(universo ~ rep, scales = "free_y") +
  scale_color_manual(values = focal_clr) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = "Efecto de la unión sobre la actividad (Hodges-Lehmann)", y = "N° de factores (ReMap)", color = NULL) +
  theme_pubr(base_size = 15) +
  theme(legend.position = "none")
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("tf_efecto_distribucion_intersect.", ext)), p2a, width = 11, height = 7.5, units = "in")
}

spec_panels <- tribble(
  ~universo, ~motivo, ~panel,
  "Cualquier línea", "CCAAT-box", "CCAAT-box (cualquier línea; NFYA ∩)",
  "HEK293", "GC-box", "GC-box (HEK293; SP1/SP2)"
)
spec_lab <- spec_pct %>%
  filter(comparacion %in% c("principal", "union laxa (referencia)")) %>%
  inner_join(spec_panels, by = c("universo", "motivo")) %>%
  mutate(
    TF_focal = case_when(TF == "NFYA_intersect" ~ "NFYA ∩", TF == "NFYA" ~ "NFYA (laxa)", TRUE ~ paste(TF, "HEK293")),
    label = sprintf("%s (p%.0f)", TF_focal, percentil)
  ) %>%
  arrange(panel, rep, desc(dif_con_sin)) %>%
  group_by(panel, rep) %>%
  mutate(vj = 1.5 + 1.5 * (row_number() - 1)) %>%
  ungroup()
p2c <- spec %>%
  inner_join(spec_panels, by = c("universo", "motivo")) %>%
  filter(!TF %in% c("NFYA", "NFYA_intersect") | universo == "HEK293") %>%
  filter(!(universo == "HEK293" & TF %in% c("SP1", "SP2"))) %>%
  ggplot(aes(dif_con_sin)) +
  geom_histogram(bins = 60, fill = "grey75", color = NA) +
  geom_vline(data = spec_lab, aes(xintercept = dif_con_sin, color = TF_focal), linewidth = 1) +
  geom_text(data = spec_lab, aes(x = dif_con_sin, y = Inf, label = label, color = TF_focal, vjust = vj), hjust = 1.03, size = 3.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  facet_grid(rep ~ panel, scales = "free") +
  scale_color_manual(values = focal_clr) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = "Efecto con motivo − efecto sin motivo (Hodges-Lehmann)", y = "N° de factores (ReMap)", color = NULL) +
  theme_pubr(base_size = 15) +
  theme(legend.position = "none")
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("tf_motivo_especificidad_intersect.", ext)), p2c, width = 11, height = 7.5, units = "in")
}

message("Tablas y figuras guardadas en ", out_dir)
