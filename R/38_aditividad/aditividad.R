# EXPLORATORIO - ¿los elementos del promotor basal se suman?
#
# Elementos con efecto positivo sobre la actividad (prom_df.tsv): isla CpG
# (CGI_frag), TATA_EPD, CCAAT_EPD, YY1_motif y CRE_motif. GCbox_EPD se
# excluye: esta muy asociada a la isla CpG (casi redundante con ella).
# INR_strong_TSS tambien queda afuera: sumado como sexto elemento casi no
# cambia el R2 (+0.001) y baja el Spearman del conteo; su efecto aparece
# solo dentro de isla CpG.
#
# Parte 1 - numero de elementos (0 a 5) por promotor. Las categorias con
#   menos de MIN_CAT promotores (en alguna replica) se agrupan con la vecina
#   ("k o mas" en la cola alta, "k o menos" en la baja). Violines de la
#   actividad media por categoria, una faceta por replica. Spearman entre
#   el numero de elementos (sin agrupar) y la actividad media, y HL (IC 95%)
#   + p de Wilcoxon de cada categoria vs. 0 elementos (BH dentro de replica).
# Parte 2 - igual que la 1 pero contando solo los cuatro motivos (sin la
#   isla), por separado en promotores con y sin isla CpG; la referencia es
#   "0 motivos" del mismo grupo (BH dentro de replica x grupo).
# Parte 3 - modelo lineal por replica: mean ~ cinco efectos principales vs.
#   mean ~ (cinco)^2 (+ las 10 interacciones de a pares). R2 de cada uno,
#   F de comparacion (anova) e interacciones con p ajustada por BH (familia
#   = las 10 de una replica). Coeficiente > 0 = sinergia (la combinacion
#   rinde mas que la suma), < 0 = subaditivo / redundante. Resumen visual:
#   para cada combinacion de los cinco elementos con >= MIN_CAT promotores,
#   actividad media observada vs. la predicha por el modelo aditivo.
#   Ademas, modelo aditivo con los cinco positivos + los tres negativos
#   (Parte 4) como efectos principales; F vs. el aditivo de solo positivos.
#   Los tres modelos se repiten con el percentil de actividad dentro de
#   cada replica como respuesta (rango promedio en empates), para ver si
#   las interacciones sobreviven a un cambio de escala o son un efecto de
#   piso de los bins. Techo de varianza explicable: Pearson (y r2) entre
#   Rep 1 y Rep 2 en los promotores medidos en ambas (con un predictor
#   perfecto del valor "verdadero", el R2 esperado contra una replica es
#   ~ r entre replicas, no r2).
# Parte 4 - caracteristicas asociadas con menor actividad: TSS no canonico
#   (other_TSS), E-box CAGGTR (Ebox_CAGGTR) y tracto poli(dA:dT) proximal
#   (polyAT_prox). Los elementos transponibles (TE_superclass no NA) quedan
#   afuera: su efecto crudo negativo desaparece al ajustar por la isla CpG
#   en el modelo lineal (los promotores con TE tienen mucha menos isla), asi
#   que contarlos repetiria la falta de isla. Conteo de
#   negativos (0 a 3, agrupado como en la Parte 1); grilla positivos x
#   negativos con la actividad mediana y n por celda (heatmap; se omiten
#   celdas con < MIN_CELDA promotores); puntaje neto = positivos -
#   negativos: violines, Spearman (puntaje sin agrupar) y HL de cada valor
#   vs. puntaje 0 (BH dentro de replica).
#
# Cada replica por separado. La figura de la tesis (slug
# aditividad_elementos_promotor en docs/mapping_figuras.csv, carpeta via
# fig_dir()) es aditividad_elementos_promotor.jpg/.pdf: A =
# violines_n_elementos (Parte 1), B = grilla_pos_x_neg (Parte 4; un heatmap
# por replica, cada uno con su propia escala de color); los dos paneles
# tambien se guardan sueltos en .jpg y .pdf. El resto es exploratorio y va a
# figures/exploratorio_aditividad/.
# Requiere: data/processed/prom_df.tsv (con YY1_motif, CRE_motif,
# Ebox_CAGGTR y polyAT_prox,
# R/00_prom_features/build_prom_features.R) y
# data/processed/activity_stats_highconf.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/tata_scan.R") # tata_hl_vs()
source("R/functions/fig_paths.R")

out_dir <- "figures/exploratorio_aditividad"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_out <- fig_dir("aditividad_elementos_promotor")

MIN_CAT <- 100
ELEMENTOS <- c(CGI_frag = "Isla CpG", TATA_EPD = "TATA", CCAAT_EPD = "CCAAT", YY1_motif = "YY1", CRE_motif = "CRE")
MOTIVOS <- ELEMENTOS[-1]
NEGATIVOS <- c(other_TSS = "TSS no canónico", Ebox_CAGGTR = "E-box CAGGTR", polyAT_prox = "Poli(dA:dT) proximal")
MIN_CELDA <- 50

th <- theme_pubclean(base_size = 16) +
  theme(strip.text = element_text(face = "bold"), legend.position = "none")
save_tsv <- function(x, name) write_tsv(x, file.path(out_dir, name))
save_fig <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(fig_out, paste0(name, ".", ext)), plot, width = width, height = height, units = "in",
      device = if (ext == "pdf") cairo_pdf else NULL)
  }
}

# --- Datos -------------------------------------------------------------------

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>%
  filter(type == "promoter") %>%
  select(seq_id, all_of(c(names(ELEMENTOS), names(NEGATIVOS))))
stopifnot(!anyNA(prom_df))

act <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE) %>%
  select(seq_id, rep, mean) %>%
  inner_join(prom_df, by = "seq_id") %>%
  mutate(
    n_elem = rowSums(across(all_of(names(ELEMENTOS)))),
    n_motivos = rowSums(across(all_of(names(MOTIVOS)))),
    n_neg = rowSums(across(all_of(names(NEGATIVOS)))),
    puntaje = n_elem - n_neg,
    grupo_cgi = factor(if_else(CGI_frag, "Con isla CpG", "Sin isla CpG"), levels = c("Con isla CpG", "Sin isla CpG"))
  ) %>%
  mutate(percentil = 100 * (rank(mean, ties.method = "average") - 0.5) / n(), .by = rep)

# Agrupa las categorias de un conteo entero con < min_n en alguna replica
# con su vecina hacia el centro. Devuelve los cortes como vector de
# etiquetas por valor (nombres = valor del conteo).
agrupar_conteo <- function(df, var, min_n = MIN_CAT) {
  tab <- df %>%
    count(rep, k = .data[[var]]) %>%
    complete(rep, k = full_seq(c(0, k), 1), fill = list(n = 0)) %>%
    group_by(k) %>%
    summarise(n = min(n), .groups = "drop") %>%
    arrange(k)
  ks <- tab$k
  lo <- 1
  hi <- length(ks)
  # cola alta: fusionar hacia abajo hasta que el bin agregado tenga >= min_n
  while (hi > lo && sum(tab$n[hi:length(ks)]) < min_n) hi <- hi - 1
  # cola baja (0 es la referencia; solo se fusiona si es necesario)
  while (lo < hi && sum(tab$n[1:lo]) < min_n) lo <- lo + 1
  lab <- as.character(ks)
  if (hi < length(ks)) lab[hi:length(ks)] <- paste(ks[hi], "o más")
  if (lo > 1) lab[1:lo] <- paste(ks[lo], "o menos")
  set_names(lab, ks)
}

aplicar_cat <- function(x, cortes) factor(cortes[as.character(x)], levels = unique(cortes))

# Conteo n / % por categoria, HL vs. la categoria ref (por defecto la
# primera) y Spearman
resumen_conteo <- function(df, var_cat, var_num, by = "rep", ref_cat = NULL) {
  conteo <- df %>%
    count(across(all_of(c(by, var_cat)))) %>%
    group_by(across(all_of(by))) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  hl <- df %>%
    group_by(across(all_of(by))) %>%
    group_modify(function(d, k) {
      lv <- levels(droplevels(d[[var_cat]]))
      ref_lv <- if (is.null(ref_cat)) lv[1] else ref_cat
      ref <- d$mean[d[[var_cat]] == ref_lv]
      n_ref <- length(ref)
      map_dfr(setdiff(lv, ref_lv), ~ tata_hl_vs(d$mean[d[[var_cat]] == .x], ref) %>%
        mutate(categoria = .x, ref = ref_lv, n = sum(d[[var_cat]] == .x), n_ref = n_ref)) %>%
        mutate(p_adj = p.adjust(p, "BH"))
    }) %>%
    ungroup() %>%
    relocate(categoria, ref, n, n_ref, .after = all_of(by))
  sp <- df %>%
    group_by(across(all_of(by))) %>%
    summarise(
      n = n(),
      rho = cor(.data[[var_num]], mean, method = "spearman"),
      p = cor.test(.data[[var_num]], mean, method = "spearman", exact = FALSE)$p.value,
      .groups = "drop"
    )
  list(conteo = conteo, hl = hl, spearman = sp)
}

violines <- function(df, var_cat, facet, xlab, fill_col) {
  lab_n <- df %>%
    group_by(across(all_of(c(setdiff(all.vars(facet), "."), var_cat)))) %>%
    summarise(n = n(), med = median(mean), .groups = "drop")
  ymax <- max(df$mean)
  ggplot(df, aes(.data[[var_cat]], mean)) +
    geom_violin(fill = fill_col, alpha = 0.7, draw_quantiles = 0.5, scale = "width") +
    geom_point(data = lab_n, aes(y = med), size = 2) +
    geom_text(data = lab_n, aes(y = ymax + 0.35, label = paste0("n=", n)), size = 3.6) +
    facet_grid(facet) +
    coord_cartesian(ylim = c(min(df$mean), ymax + 0.5)) +
    labs(x = xlab, y = "Actividad media") +
    th
}

# =============================================================================
# PARTE 1 - numero de elementos
# =============================================================================

cortes_elem <- agrupar_conteo(act, "n_elem")
act <- act %>% mutate(cat_elem = aplicar_cat(n_elem, cortes_elem))
r1 <- resumen_conteo(act, "cat_elem", "n_elem")
save_tsv(r1$conteo %>% rename(categoria = cat_elem), "p1_conteo_n_elementos.tsv")
save_tsv(r1$hl, "p1_hl_vs_0_elementos.tsv")
save_tsv(r1$spearman, "p1_spearman_n_elementos.tsv")

p1 <- violines(act, "cat_elem", . ~ rep, "Número de elementos (isla CpG, TATA, CCAAT, YY1, CRE)", "#216869")
save_fig(p1, "violines_n_elementos", 12, 6)

# =============================================================================
# PARTE 2 - motivos dentro de cada grupo de isla CpG
# =============================================================================

act2 <- act %>%
  group_split(grupo_cgi) %>%
  map_dfr(function(d) {
    cortes <- agrupar_conteo(d, "n_motivos")
    d %>% mutate(cat_motivos = cortes[as.character(n_motivos)])
  }) %>%
  mutate(cat_motivos = factor(cat_motivos, levels = unique(cat_motivos[order(n_motivos)])))
r2 <- resumen_conteo(act2, "cat_motivos", "n_motivos", by = c("rep", "grupo_cgi"))
save_tsv(r2$conteo %>% rename(categoria = cat_motivos), "p2_conteo_n_motivos_por_cgi.tsv")
save_tsv(r2$hl, "p2_hl_vs_0_motivos_por_cgi.tsv")
save_tsv(r2$spearman, "p2_spearman_n_motivos_por_cgi.tsv")

p2 <- violines(act2, "cat_motivos", grupo_cgi ~ rep, "Número de motivos (TATA, CCAAT, YY1, CRE)", "#14AFB2") +
  scale_x_discrete(drop = TRUE)
p2 <- p2 + facet_grid(rep ~ grupo_cgi, scales = "free_x", space = "free_x")
ggsave(file.path(out_dir, "p2_violines_n_motivos_por_cgi.jpg"), p2, width = 13, height = 9, units = "in")

# =============================================================================
# PARTE 3 - ¿aditivo o con interaccion?
# =============================================================================

# Ajusta los tres modelos por replica con la respuesta resp
ajustar_modelos <- function(df, resp) {
  f_adit <- reformulate(names(ELEMENTOS), resp)
  f_inter <- reformulate(paste0("(", paste(names(ELEMENTOS), collapse = " + "), ")^2"), resp)
  f_pos_neg <- reformulate(c(names(ELEMENTOS), names(NEGATIVOS)), resp)
  df %>%
    group_by(rep) %>%
    group_map(function(d, k) {
      m_a <- lm(f_adit, data = d)
      m_i <- lm(f_inter, data = d)
      av <- anova(m_a, m_i)
      m_pn <- lm(f_pos_neg, data = d)
      av_pn <- anova(m_a, m_pn)
      list(
        rep = k$rep, d = d, m_a = m_a,
        ajuste = tibble(
          rep = k$rep, respuesta = resp, n = nrow(d),
          r2_aditivo = summary(m_a)$r.squared, r2_adj_aditivo = summary(m_a)$adj.r.squared,
          r2_interaccion = summary(m_i)$r.squared, r2_adj_interaccion = summary(m_i)$adj.r.squared,
          F = av$F[2], df1 = av$Df[2], df2 = av$Res.Df[2], p_F = av$`Pr(>F)`[2],
          r2_pos_neg = summary(m_pn)$r.squared, r2_adj_pos_neg = summary(m_pn)$adj.r.squared,
          F_pos_neg = av_pn$F[2], p_F_pos_neg = av_pn$`Pr(>F)`[2]
        ),
        pos_neg = broom::tidy(m_pn, conf.int = TRUE) %>% mutate(rep = k$rep, .before = 1),
        principales = broom::tidy(m_a, conf.int = TRUE) %>% mutate(rep = k$rep, .before = 1),
        interacciones = broom::tidy(m_i, conf.int = TRUE) %>%
          filter(str_detect(term, ":")) %>%
          mutate(rep = k$rep, .before = 1)
      )
    })
}
modelos <- ajustar_modelos(act, "mean")

etiqueta_term <- function(term) {
  term %>%
    str_remove_all("TRUE") %>%
    str_split(":") %>%
    map_chr(~ paste(c(ELEMENTOS, NEGATIVOS)[.x], collapse = " × "))
}

ajuste <- map_dfr(modelos, "ajuste")
tabla_interacciones <- function(modelos) {
  map_dfr(modelos, "interacciones") %>%
    group_by(rep) %>%
    mutate(p_adj = p.adjust(p.value, "BH")) %>%
    ungroup() %>%
    mutate(
      interaccion = etiqueta_term(term), .after = term,
      direccion = case_when(
        p_adj < 0.05 & estimate > 0 ~ "Sinergia (> suma)",
        p_adj < 0.05 & estimate < 0 ~ "Subaditiva (< suma)",
        TRUE ~ "n.s."
      )
    ) %>%
    arrange(rep, p_adj)
}
principales <- map_dfr(modelos, "principales") %>%
  mutate(elemento = if_else(term == "(Intercept)", "(Intercepto)", etiqueta_term(term)), .after = term)
interacciones <- tabla_interacciones(modelos)
save_tsv(ajuste, "p3_modelos_r2_F.tsv")
save_tsv(principales, "p3_coef_aditivo.tsv")

# Modelo aditivo con positivos + negativos: coeficientes
pos_neg <- map_dfr(modelos, "pos_neg") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    elemento = etiqueta_term(term), .after = term,
    clase = if_else(str_remove(term, "TRUE$") %in% names(ELEMENTOS), "Positivo", "Negativo"),
    p_adj = p.adjust(p.value, "BH"), .by = rep
  ) %>%
  relocate(p_adj, .after = last_col())
save_tsv(pos_neg, "p3_coef_aditivo_pos_neg.tsv")
p3c <- pos_neg %>%
  mutate(elemento = factor(elemento, levels = c(rev(NEGATIVOS), rev(ELEMENTOS)))) %>%
  ggplot(aes(estimate, elemento, color = clase)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  facet_wrap(~rep) +
  scale_color_manual(values = c(Positivo = "#D6741F", Negativo = "#216869"), name = NULL) +
  labs(x = "Efecto sobre la actividad media (coeficiente, IC 95%)", y = NULL) +
  th + theme(legend.position = "top")
ggsave(file.path(out_dir, "p3_coef_aditivo_pos_neg.jpg"), p3c, width = 12, height = 6, units = "in")
save_tsv(interacciones, "p3_interacciones.tsv")

# Interacciones: forest plot de los 10 coeficientes por replica
orden_int <- interacciones %>%
  group_by(interaccion) %>%
  summarise(e = mean(estimate)) %>%
  arrange(e) %>%
  pull(interaccion)
p3a <- interacciones %>%
  mutate(interaccion = factor(interaccion, levels = orden_int), direccion = factor(direccion, levels = c("Sinergia (> suma)", "Subaditiva (< suma)", "n.s."))) %>%
  ggplot(aes(estimate, interaccion, color = direccion)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  facet_wrap(~rep) +
  scale_color_manual(values = c("Sinergia (> suma)" = "#D6741F", "Subaditiva (< suma)" = "#216869", "n.s." = "grey60"), drop = TRUE, name = NULL) +
  labs(x = "Coeficiente de interacción (IC 95%)", y = NULL) +
  th + theme(legend.position = "top", panel.spacing.x = unit(2, "lines"))
ggsave(file.path(out_dir, "p3_interacciones.jpg"), p3a, width = 12, height = 7, units = "in")

# Observado vs. predicho por el modelo aditivo, por combinacion
obs_pred <- map_dfr(modelos, function(m) {
  m$d %>%
    mutate(rep = m$rep, pred = fitted(m$m_a)) %>% # group_map quita la columna rep
    group_by(rep, across(all_of(names(ELEMENTOS)))) %>%
    summarise(n = n(), n_elem = first(n_elem), observada = mean(mean), predicha = mean(pred), .groups = "drop")
}) %>%
  mutate(
    combinacion = pmap_chr(across(all_of(names(ELEMENTOS))), function(...) {
      v <- c(...)
      if (!any(v)) "Ninguno" else paste(ELEMENTOS[names(v)[v]], collapse = " + ")
    }),
    residuo = observada - predicha
  ) %>%
  relocate(combinacion, n, .after = rep)
save_tsv(obs_pred %>% arrange(rep, desc(n)), "p3_obs_vs_pred_combinaciones.tsv")

obs_pred_f <- obs_pred %>% filter(n >= MIN_CAT)
rmse <- obs_pred_f %>%
  group_by(rep) %>%
  summarise(
    n_comb = n(), rmse = sqrt(mean(residuo^2)), r = cor(observada, predicha),
    lab = sprintf("%d combinaciones (n ≥ %d)\nr = %.3f, RMSE = %.3f", n_comb, MIN_CAT, r, rmse), .groups = "drop"
  )
save_tsv(rmse %>% select(-lab), "p3_obs_vs_pred_resumen.tsv")
lims <- range(c(obs_pred_f$observada, obs_pred_f$predicha)) + c(-0.1, 0.1)
p3b <- ggplot(obs_pred_f, aes(predicha, observada)) +
  geom_abline(linetype = 2, color = "grey50") +
  geom_point(aes(size = n, fill = factor(n_elem)), shape = 21, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = combinacion), size = 2.8, max.overlaps = 15, min.segment.length = 0, seed = 1) +
  geom_text(data = rmse, aes(x = lims[1], y = lims[2], label = lab), hjust = 0, vjust = 1, size = 3.8, inherit.aes = FALSE) +
  facet_wrap(~rep) +
  coord_equal(xlim = lims, ylim = lims) +
  scale_fill_viridis_d(option = "mako", direction = -1, name = "N.º de elementos") +
  scale_size_area(max_size = 8, name = "Promotores") +
  labs(x = "Actividad media predicha (modelo aditivo)", y = "Actividad media observada") +
  th + theme(legend.position = "right")
ggsave(file.path(out_dir, "p3_obs_vs_pred_aditivo.jpg"), p3b, width = 15, height = 7.5, units = "in")

# --- Mismos modelos con el percentil de actividad como respuesta -------------
# Si las sinergias fueran un efecto de piso de la escala de bins (muchos
# promotores en el minimo), deberian achicarse o desaparecer al pasar a
# percentiles. Coeficientes en puntos de percentil.
modelos_pct <- ajustar_modelos(act, "percentil")
ajuste_pct <- map_dfr(modelos_pct, "ajuste")
save_tsv(bind_rows(ajuste, ajuste_pct), "p3_modelos_r2_F_media_vs_percentil.tsv")
interacciones_pct <- tabla_interacciones(modelos_pct)
save_tsv(interacciones_pct, "p3_interacciones_percentil.tsv")
comp_int <- bind_rows(
  interacciones %>% mutate(respuesta = "Actividad media"),
  interacciones_pct %>% mutate(respuesta = "Percentil")
) %>%
  select(respuesta, rep, interaccion, estimate, conf.low, conf.high, p.value, p_adj, direccion) %>%
  pivot_wider(names_from = respuesta, values_from = c(estimate, conf.low, conf.high, p.value, p_adj, direccion), names_glue = "{.value}_{respuesta}") %>%
  arrange(rep, `p_adj_Actividad media`)
save_tsv(comp_int, "p3_interacciones_media_vs_percentil.tsv")
pos_neg_pct <- map_dfr(modelos_pct, "pos_neg") %>%
  filter(term != "(Intercept)") %>%
  mutate(elemento = etiqueta_term(term), .after = term, p_adj = p.adjust(p.value, "BH"), .by = rep) %>%
  relocate(p_adj, .after = last_col())
save_tsv(pos_neg_pct, "p3_coef_aditivo_pos_neg_percentil.tsv")

# --- Techo de varianza explicable: acuerdo entre replicas --------------------
techo <- act %>%
  select(seq_id, rep, mean, percentil) %>%
  pivot_wider(names_from = rep, values_from = c(mean, percentil)) %>%
  drop_na()
techo <- tibble(
  escala = c("Actividad media", "Percentil"),
  n_ambas = nrow(techo),
  pearson_r = c(cor(techo$`mean_Rep 1`, techo$`mean_Rep 2`), cor(techo$`percentil_Rep 1`, techo$`percentil_Rep 2`)),
  spearman_rho = cor(techo$`mean_Rep 1`, techo$`mean_Rep 2`, method = "spearman")
) %>%
  mutate(r2 = pearson_r^2)
save_tsv(techo, "p3_techo_acuerdo_replicas.tsv")

# =============================================================================
# PARTE 4 - caracteristicas negativas y puntaje neto
# =============================================================================

cortes_neg <- agrupar_conteo(act, "n_neg")
act <- act %>% mutate(cat_neg = aplicar_cat(n_neg, cortes_neg))
r4 <- resumen_conteo(act, "cat_neg", "n_neg")
save_tsv(r4$conteo %>% rename(categoria = cat_neg), "p4_conteo_n_negativos.tsv")
save_tsv(r4$hl, "p4_hl_vs_0_negativos.tsv")
save_tsv(r4$spearman, "p4_spearman_n_negativos.tsv")

# Grilla positivos x negativos
grilla <- act %>%
  group_by(rep, cat_elem, cat_neg) %>%
  summarise(n = n(), mediana = median(mean), media = mean(mean), .groups = "drop")
save_tsv(grilla, "p4_grilla_pos_x_neg.tsv")
# Un heatmap por replica, cada uno con su propia escala de color (las
# replicas tienen niveles de actividad distintos); se omiten celdas con
# menos de MIN_CELDA promotores.
heatmap_rep <- function(r) {
  grilla %>%
    filter(n >= MIN_CELDA, rep == r) %>%
    ggplot(aes(cat_elem, cat_neg, fill = mediana)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f
n=%d", mediana, n)), size = 4.2) +
    facet_wrap(~rep) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Actividad
mediana") +
    labs(x = "Elementos positivos
(isla CpG, TATA, CCAAT, YY1, CRE)", y = "Características negativas
(TSS no canónico, E-box, poli(dA:dT))") +
    th + theme(legend.position = "right", panel.grid = element_blank())
}
p4a <- heatmap_rep("Rep 1") + (heatmap_rep("Rep 2") + theme(axis.title.y = element_blank()))
save_fig(p4a, "grilla_pos_x_neg", 15, 6.5)

# Figura completa: violines (A) arriba, heatmaps (B) abajo
figura <- wrap_elements(full = p1) / wrap_elements(full = p4a) +
  plot_layout(heights = c(6, 6.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 26, face = "bold"))
save_fig(figura, "aditividad_elementos_promotor", 15, 13)

# Puntaje neto = positivos - negativos
cortes_punt <- agrupar_conteo(act, "puntaje")
act <- act %>% mutate(cat_puntaje = aplicar_cat(puntaje, cortes_punt))
ref_punt <- cortes_punt[["0"]]
r4p <- resumen_conteo(act, "cat_puntaje", "puntaje", ref_cat = ref_punt)
save_tsv(r4p$conteo %>% rename(categoria = cat_puntaje), "p4_conteo_puntaje_neto.tsv")
save_tsv(r4p$hl, "p4_hl_vs_0_puntaje_neto.tsv")
save_tsv(r4p$spearman, "p4_spearman_puntaje_neto.tsv")
p4b <- violines(act, "cat_puntaje", . ~ rep, "Puntaje neto (elementos positivos − características negativas)", "#216869")
ggsave(file.path(out_dir, "p4_violines_puntaje_neto.jpg"), p4b, width = 13, height = 6, units = "in")

message("Salidas en ", out_dir)
