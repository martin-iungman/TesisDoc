# Exploratorio: panorama de caracteristicas asociadas al ruido
# (ruido_summary) estratificado por isla CpG (CGI_frag, la definicion por
# composicion del fragmento que usa todo R3).
#
# Motivo: en figures/exploratorio_ausencia_senal/handoff.md el gradiente de
# ocupacion ReMap y el conteo de ausencias viven casi enteramente en los
# promotores sin isla. Aca se pregunta lo mismo para cada caracteristica
# del panorama.
#
# p1 AUC (pROC, igual que ruido_summary: P(ruido | con > ruido | sin)),
#    IC 95% y p de Wald con la varianza de DeLong, dentro de cada estrato
#    (con / sin isla) y replica. BH dentro de cada estrato (sobre
#    caracteristicas x replicas), igual que el BH global de ruido_summary
#    es sobre caracteristicas x replicas.
# p2 n con / sin la caracteristica en cada celda. Minimo: N_MIN promotores
#    en cada grupo dentro del estrato; las celdas que no llegan quedan sin
#    evaluar (n reportado igual).
# p3 Tabla comparativa global (ruido_summary_auc.tsv, se recalcula y se
#    verifica que coincida) vs. cada estrato, Rep 1 / Rep 2; diferencia
#    con isla - sin isla (z con la varianza de DeLong de cada estrato,
#    muestras independientes), BH sobre caracteristicas x replicas.
# p4 Criterio de consistencia (pval_corr < 0.05 y mismo sentido en ambas
#    replicas) global y por estrato; cambios de estado; compatibilidad del
#    IC del estrato con el AUC global.
# p5 Sensibilidad: var_rank_sw principal = rango global (add_noise_rank
#    sobre el set completo, como R3) restringido a cada estrato; version
#    alternativa = add_noise_rank recalculado dentro de cada estrato y
#    replica (bins de 100 dentro del estrato).
# p6 Caracteristicas ligadas a la definicion de isla: prevalencia dentro de
#    cada estrato.
# p7 Figura: panorama con global / con isla / sin isla por caracteristica,
#    seq y endo separados como en ruido_summary, con n visible.
#
# Las caracteristicas se definen sobre el set completo (cortes de tercil
# globales de build_tidy_features()) y despues se estratifica, asi que
# "Alta especificidad tisular" etc. significan lo mismo que en el global.
#
# Salida en figures/exploratorio_ruido_summary_cgi/ (prefijo por punto) y
# handoff.md. Requiere data/processed/activity_stats_highconf.tsv y
# prom_df.tsv. Corre en segundos. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(pROC)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter

out_dir <- "figures/exploratorio_ruido_summary_cgi"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

N_MIN <- 50
estratos <- c("Global", "Con isla CpG", "Sin isla CpG")

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d, big.mark = "")

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>% add_noise_rank()

# Rango recalculado dentro de cada estrato (sensibilidad, p5). add_mean_sw_
# bins descarta el resto de la division por 100 dentro de cada estrato, asi
# que hasta 99 promotores por estrato y replica quedan sin este rango.
rank_est <- data %>%
  select(seq_id, rep, mean, var, CGI_frag) %>%
  group_split(CGI_frag) %>%
  map(add_noise_rank) %>%
  list_rbind() %>%
  select(seq_id, rep, var_rank_est = var_rank_sw)
data <- left_join(data, rank_est, by = c("seq_id", "rep"))
stopifnot(!anyNA(data$CGI_frag))

tidy_data <- build_tidy_features(data, keep = c("var_rank_sw", "var_rank_est")) %>%
  mutate(cgi = `Islas CpG` == "TRUE")
vbles <- setdiff(names(tidy_data), c("rep", "var_rank_sw", "var_rank_est", "cgi"))
repname <- levels(tidy_data$rep)

# --- AUC por celda ---------------------------------------------------------
auc_celda <- function(df, f, score) {
  resp <- df[[f]]
  ok <- !is.na(resp) & !is.na(df[[score]])
  n_con <- sum(resp[ok] == "TRUE")
  n_sin <- sum(resp[ok] == "FALSE")
  base <- tibble(feature = f, n_con = n_con, n_sin = n_sin)
  if (min(n_con, n_sin) < N_MIN) {
    return(base %>% mutate(evaluada = FALSE, AUC = NA_real_, ci2.5 = NA_real_, ci97.5 = NA_real_, se = NA_real_, pval = NA_real_))
  }
  roc_obj <- pROC::roc(resp[ok], df[[score]][ok], ci = TRUE, direction = ">", quiet = TRUE)
  auc_val <- as.numeric(roc_obj$ci[2])
  se_delong <- sqrt(pROC::var(roc_obj, method = "delong"))
  base %>% mutate(
    evaluada = TRUE, AUC = auc_val, ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]),
    se = se_delong, pval = 2 * pnorm(-abs((auc_val - 0.5) / se_delong))
  )
}

correr <- function(score) {
  celdas <- list(
    "Global" = tidy_data,
    "Con isla CpG" = filter(tidy_data, cgi),
    "Sin isla CpG" = filter(tidy_data, !cgi)
  )
  imap(celdas, function(df_e, est) {
    map(repname, function(r) {
      df <- filter(df_e, rep == r)
      map(vbles, ~ auc_celda(df, .x, score)) %>% list_rbind() %>% mutate(rep = r, n_estrato = nrow(df))
    }) %>% list_rbind() %>% mutate(estrato = est)
  }) %>%
    list_rbind() %>%
    group_by(estrato) %>%
    mutate(pval_corr = p.adjust(pval, method = "BH")) %>%
    ungroup() %>%
    mutate(estrato = factor(estrato, levels = estratos)) %>%
    select(estrato, feature, rep, n_estrato, n_con, n_sin, evaluada, AUC, ci2.5, ci97.5, se, pval, pval_corr)
}

# Caracteristicas excluidas a priori del analisis estratificado (p6): la
# isla misma no tiene variacion dentro del estrato por definicion.
excl_definicion <- "Islas CpG"

auc_main <- correr("var_rank_sw") %>%
  mutate(evaluada = evaluada & !(estrato != "Global" & feature %in% excl_definicion)) %>%
  mutate(across(c(AUC, ci2.5, ci97.5, se, pval), ~ if_else(evaluada, .x, NA_real_))) %>%
  group_by(estrato) %>%
  mutate(pval_corr = p.adjust(pval, method = "BH")) %>%
  ungroup()

# Verificacion contra el global de ruido_summary
auc_ref <- read_tsv(file.path(fig_dir("ruido_summary"), "ruido_summary_auc.tsv"), show_col_types = FALSE)
chk <- auc_main %>% filter(estrato == "Global") %>% inner_join(auc_ref, by = c("feature", "rep"), suffix = c("", "_ref"))
stopifnot(nrow(chk) == nrow(auc_ref), max(abs(chk$AUC - chk$AUC_ref)) < 1e-9, max(abs(chk$pval_corr - chk$pval_corr_ref)) < 1e-9)
message("Global recalculado coincide con ruido_summary_auc.tsv")

# --- p6: prevalencia de las caracteristicas ligadas a la isla --------------
p6 <- auc_main %>%
  mutate(prevalencia = n_con / (n_con + n_sin)) %>%
  select(feature, estrato, rep, n_con, n_sin, prevalencia, evaluada) %>%
  arrange(feature, estrato, rep)
write_tsv(p6, file.path(out_dir, "p6_prevalencia_por_estrato.tsv"))

# --- p1/p2: tabla larga por celda ------------------------------------------
consist <- function(df) {
  df %>%
    group_by(estrato, feature) %>%
    summarise(
      evaluada_ambas = all(evaluada),
      sig_ambas = evaluada_ambas && all(pval_corr < 0.05, na.rm = FALSE),
      mismo_sentido = evaluada_ambas && n_distinct(sign(AUC - 0.5)) == 1,
      consistente = sig_ambas & mismo_sentido,
      sentido = if_else(consistente, if_else(mean(AUC) > 0.5, "Ruido alto", "Ruido bajo"), NA_character_),
      .groups = "drop"
    ) %>%
    mutate(estado = case_when(
      !evaluada_ambas ~ "sin evaluar",
      consistente ~ sentido,
      TRUE ~ "no consistente"
    ))
}
cons_main <- consist(auc_main)

p1 <- auc_main %>% left_join(cons_main %>% select(estrato, feature, consistente, estado), by = c("estrato", "feature"))
write_tsv(p1, file.path(out_dir, "p1_p2_auc_por_estrato.tsv"))
p2_sin_evaluar <- auc_main %>%
  filter(!evaluada, estrato != "Global") %>%
  mutate(motivo = if_else(feature %in% excl_definicion, "define el estrato", paste0("min(n con, n sin) < ", N_MIN))) %>%
  select(feature, estrato, rep, n_con, n_sin, motivo) %>%
  arrange(estrato, feature, rep)
write_tsv(p2_sin_evaluar, file.path(out_dir, "p2_sin_evaluar.tsv"))

# --- p3: comparativa global vs. estratos -----------------------------------
glob <- auc_main %>% filter(estrato == "Global") %>% select(feature, rep, AUC_glob = AUC)
ancho <- auc_main %>%
  filter(estrato != "Global") %>%
  left_join(glob, by = c("feature", "rep")) %>%
  mutate(compatible_global = ci2.5 <= AUC_glob & ci97.5 >= AUC_glob)

dif <- ancho %>%
  select(feature, rep, estrato, AUC, se) %>%
  pivot_wider(names_from = estrato, values_from = c(AUC, se)) %>%
  mutate(
    dif_con_menos_sin = `AUC_Con isla CpG` - `AUC_Sin isla CpG`,
    z = dif_con_menos_sin / sqrt(`se_Con isla CpG`^2 + `se_Sin isla CpG`^2),
    p_dif = 2 * pnorm(-abs(z))
  ) %>%
  mutate(p_dif_corr = p.adjust(p_dif, method = "BH")) %>%
  select(feature, rep, dif_con_menos_sin, p_dif, p_dif_corr)

# AUC combinado dentro de estratos ("ajustado por isla"): promedio de los
# AUC de los dos estratos ponderado por varianza inversa (DeLong), solo
# cuando la caracteristica se evaluo en ambos. Tiene mas potencia que cada
# estrato por separado: sirve para distinguir "no significativo en el
# estrato sin isla por falta de n" de "el efecto desaparece al comparar
# dentro de la misma clase de isla".
comb <- ancho %>%
  group_by(feature, rep) %>%
  filter(all(evaluada), n() == 2) %>%
  summarise(w = list(1 / se^2), AUC = sum(AUC / se^2) / sum(1 / se^2), se = sqrt(1 / sum(1 / se^2)), .groups = "drop") %>%
  select(-w) %>%
  mutate(ci2.5 = AUC - 1.96 * se, ci97.5 = AUC + 1.96 * se, pval = 2 * pnorm(-abs((AUC - 0.5) / se)), pval_corr = p.adjust(pval, method = "BH"))
cons_comb <- comb %>%
  group_by(feature) %>%
  summarise(estado_comb = if (n() == 2 && all(pval_corr < 0.05) && n_distinct(sign(AUC - 0.5)) == 1) if_else(mean(AUC) > 0.5, "Ruido alto", "Ruido bajo") else "no consistente", .groups = "drop")
write_tsv(comb %>% left_join(cons_comb, by = "feature") %>% left_join(glob, by = c("feature", "rep")), file.path(out_dir, "p3_auc_combinado_dentro_de_estratos.tsv"))

celda_txt <- function(auc, lo, hi, p) if_else(is.na(auc), "—", paste0(fmt(auc), " [", fmt(lo), ", ", fmt(hi), "]", if_else(p < 0.05, "*", "")))
p3 <- auc_main %>%
  mutate(txt = celda_txt(AUC, ci2.5, ci97.5, pval_corr), n_txt = paste0(n_con, "/", n_sin)) %>%
  select(feature, estrato, rep, txt, n_txt) %>%
  pivot_wider(names_from = c(estrato, rep), values_from = c(txt, n_txt), names_glue = "{estrato} {rep} {.value}") %>%
  left_join(cons_main %>% select(estrato, feature, estado) %>% pivot_wider(names_from = estrato, values_from = estado, names_glue = "estado {estrato}"), by = "feature") %>%
  left_join(ancho %>% select(feature, rep, estrato, compatible_global) %>% pivot_wider(names_from = c(estrato, rep), values_from = compatible_global, names_glue = "compatible_global {estrato} {rep}"), by = "feature") %>%
  left_join(dif %>% pivot_wider(names_from = rep, values_from = c(dif_con_menos_sin, p_dif, p_dif_corr), names_glue = "{.value} {rep}"), by = "feature") %>%
  left_join(comb %>% mutate(txt = celda_txt(AUC, ci2.5, ci97.5, pval_corr)) %>% select(feature, rep, txt) %>% pivot_wider(names_from = rep, values_from = txt, names_glue = "Combinado {rep} txt"), by = "feature") %>%
  left_join(cons_comb, by = "feature") %>%
  left_join(feature_groups, by = "feature") %>%
  relocate(group, .after = feature) %>%
  arrange(group, feature)
write_tsv(p3, file.path(out_dir, "p3_comparativa_global_estratos.tsv"))

# --- p4: consistencia y cambios de estado ----------------------------------
p4_resumen <- cons_main %>%
  group_by(estrato) %>%
  summarise(
    evaluadas = sum(evaluada_ambas), consistentes = sum(consistente),
    ruido_alto = sum(estado == "Ruido alto"), ruido_bajo = sum(estado == "Ruido bajo"), .groups = "drop"
  )
write_tsv(p4_resumen, file.path(out_dir, "p4_consistencia_resumen.tsv"))
p4_estado <- cons_main %>%
  select(feature, estrato, estado) %>%
  pivot_wider(names_from = estrato, values_from = estado) %>%
  left_join(feature_groups, by = "feature") %>%
  mutate(cambia = `Con isla CpG` != Global | `Sin isla CpG` != Global) %>%
  arrange(group, feature)
write_tsv(p4_estado, file.path(out_dir, "p4_consistencia_estado.tsv"))

# --- p5: sensibilidad con rango recalculado dentro del estrato -------------
auc_est <- correr("var_rank_est") %>%
  filter(estrato != "Global") %>%
  mutate(evaluada = evaluada & !(feature %in% excl_definicion)) %>%
  mutate(across(c(AUC, ci2.5, ci97.5, se, pval), ~ if_else(evaluada, .x, NA_real_))) %>%
  group_by(estrato) %>%
  mutate(pval_corr = p.adjust(pval, method = "BH")) %>%
  ungroup()
cons_est <- consist(auc_est)
p5 <- auc_main %>%
  filter(estrato != "Global") %>%
  select(estrato, feature, rep, n_con, n_sin, AUC, ci2.5, ci97.5, pval_corr) %>%
  left_join(auc_est %>% select(estrato, feature, rep, n_con_est = n_con, n_sin_est = n_sin, AUC_est = AUC, ci2.5_est = ci2.5, ci97.5_est = ci97.5, pval_corr_est = pval_corr), by = c("estrato", "feature", "rep")) %>%
  left_join(cons_main %>% select(estrato, feature, estado), by = c("estrato", "feature")) %>%
  left_join(cons_est %>% select(estrato, feature, estado_est = estado), by = c("estrato", "feature")) %>%
  mutate(dif_AUC = AUC_est - AUC)
write_tsv(p5, file.path(out_dir, "p5_sensibilidad_rango_por_estrato.tsv"))
p5_resumen <- p5 %>%
  filter(!is.na(AUC), !is.na(AUC_est)) %>%
  group_by(estrato) %>%
  summarise(
    n_celdas = n(), cor_pearson = cor(AUC, AUC_est), max_abs_dif = max(abs(dif_AUC)), mediana_abs_dif = median(abs(dif_AUC)),
    cambios_sig = sum((pval_corr < 0.05) != (pval_corr_est < 0.05)), .groups = "drop"
  ) %>%
  left_join(
    p5 %>% distinct(estrato, feature, estado, estado_est) %>% group_by(estrato) %>% summarise(cambios_estado = sum(estado != estado_est), .groups = "drop"),
    by = "estrato"
  )
write_tsv(p5_resumen, file.path(out_dir, "p5_sensibilidad_resumen.tsv"))

# --- p7: figura -------------------------------------------------------------
# Filas: caracteristicas consistentes en al menos uno de los tres conjuntos.
# Columnas: global / con isla / sin isla. Barras AUC - 0.5 con IC 95%,
# color por sentido si la celda (caracteristica x estrato) cumple el
# criterio de consistencia, gris si no; alpha por replica como en
# ruido_summary. A la derecha de cada columna, n con la caracteristica
# (Rep 1 / Rep 2); "n < 50" si no se evaluo.
filas <- cons_main %>% filter(consistente) %>% distinct(feature) %>% pull(feature)
plot_df <- auc_main %>%
  filter(feature %in% filas) %>%
  left_join(cons_main %>% select(estrato, feature, estado), by = c("estrato", "feature")) %>%
  mutate(fill = case_when(estado %in% c("Ruido alto", "Ruido bajo") ~ estado, TRUE ~ "n.s. / no consistente"))
orden <- plot_df %>% filter(estrato == "Global") %>% group_by(feature) %>% summarise(m = mean(AUC)) %>% arrange(m) %>% pull(feature)
x_lim <- range(c(plot_df$ci2.5, plot_df$ci97.5, 0.5), na.rm = TRUE) - 0.5
x_txt <- x_lim[2] + 0.02
n_lab <- plot_df %>%
  group_by(estrato, feature) %>%
  summarise(lab = if (all(evaluada)) paste(n_con, collapse = " / ") else if (feature[1] %in% excl_definicion) "—" else paste0("n < ", N_MIN, " (", paste(n_con, collapse = " / "), ")"), .groups = "drop")

panel_cgi <- function(grupo, y_label) {
  fs <- intersect(orden, feature_groups$feature[feature_groups$group == grupo])
  d <- plot_df %>% filter(feature %in% fs) %>% mutate(feature = factor(feature, levels = fs), rep = factor(rep, levels = c("Rep 2", "Rep 1")))
  nl <- n_lab %>% filter(feature %in% fs) %>% mutate(feature = factor(feature, levels = fs))
  ggplot(d, aes(x = AUC - 0.5, y = feature, group = rep)) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    geom_col(orientation = "y", position = position_dodge(0.9), width = 0.9, aes(fill = fill, alpha = rep), na.rm = TRUE) +
    geom_errorbarh(aes(xmin = ci2.5 - 0.5, xmax = ci97.5 - 0.5), position = position_dodge(0.9), height = 0.05, col = "#777777", linewidth = 1, na.rm = TRUE) +
    geom_text(data = nl, aes(x = x_txt, y = feature, label = lab), inherit.aes = FALSE, hjust = 0, size = 4.2, col = "grey30") +
    facet_wrap(~estrato, nrow = 1) +
    scale_x_continuous(labels = function(x) fmt(x + 0.5, 2), breaks = seq(-0.2, 0.2, 0.1)) +
    coord_cartesian(xlim = c(x_lim[1], x_lim[2] + 0.13), clip = "off") +
    scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7), breaks = c("Rep 1", "Rep 2")) +
    scale_fill_manual(values = c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800", "n.s. / no consistente" = "grey75"), drop = FALSE) +
    theme_pubr(base_size = 18) +
    theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(size = 18, face = "bold"), panel.spacing = unit(1.5, "lines")) +
    labs(fill = "Efecto", x = "AUC (efecto sobre el ruido)", y = y_label, alpha = "")
}
panel_seq <- panel_cgi("seq", "Elementos del promotor")
panel_endo <- panel_cgi("endo", "Características del promotor")
n_feat <- map_int(list(panel_seq, panel_endo), ~ n_distinct(.x$data$feature))
figura <- (panel_seq / panel_endo) +
  plot_layout(heights = n_feat + 3, guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    caption = paste0("n = promotores con la característica (Rep 1 / Rep 2). Color: significativo (BH < 0.05) y mismo sentido en ambas réplicas dentro del estrato. Celdas con menos de ", N_MIN, " promotores en algún grupo, sin evaluar.")
  ) &
  theme(legend.position = "top", plot.tag = element_text(size = 26, face = "bold"), plot.caption = element_text(size = 12, hjust = 0))
save_jpg_pdf(figura, "p7_summary_noise_cgi", 19, 2.5 + 0.5 * sum(n_feat + 3))
write_tsv(auc_main %>% left_join(cons_main %>% select(estrato, feature, estado), by = c("estrato", "feature")) %>% left_join(feature_groups, by = "feature"),
  file.path(out_dir, "p7_summary_noise_cgi.tsv"))

# --- Resumen en consola -----------------------------------------------------
print(count(tidy_data, rep, cgi))
print(p4_resumen)
print(p5_resumen)
print(p4_estado, n = Inf)
print(p2_sin_evaluar, n = Inf)
message("Salidas en ", out_dir)
