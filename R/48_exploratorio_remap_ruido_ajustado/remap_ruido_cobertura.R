# Exploratorio complementario de remap_ruido_ajustado.R (punto 7 del mismo
# pedido), en la misma carpeta de salida: control por amplitud de union.
#
# El efecto de cada factor sobre el ruido correlaciona con cuantos promotores
# une, y los miembros de MLL1/2 unen muchos mas promotores que el factor
# tipico. Un factor que une el 70-80% de la library define un grupo no unido
# que es el 20-30% mas atipico (sin modulo cis-regulatorio, sin actividad
# endogena, sin isla CpG), justo los promotores mas ruidosos. Controles:
#
# p7.1 Ventana de cobertura (mismo procedimiento que nfy_especificidad,
#   R/16_tf_contexto_endogeno/nfy_especificidad.R): para cada miembro de
#   MLL1/2 (GO:0044665), posicion de su coeficiente dentro de los otros 981
#   factores y dentro de los que unen un n de promotores a +-25% del suyo.
#   Rango 1 = coeficiente mas negativo (mayor reduccion del ruido); percentil
#   = % de factores comparables con coeficiente <= el del focal (bajo = mas
#   reduccion). Segunda referencia que excluye ademas a todos los miembros de
#   MLL1/2. Se calcula para todos los modelos; el pedido es sobre `ajustado`.
# p7.2 Figura: coeficiente ajustado vs. n unidos (log10), por replica, una
#   columna por miembro del nucleo COMPASS con su ventana sombreada.
# p7.3 Ocupacion total como covariable: modelo ajustado + log2(n + 1), con n
#   = factores distintos del catalogo ReMap completo (1205) con pico sobre el
#   promotor (incluye al propio factor; en escala log2 sobre ~400 factores su
#   aporte es despreciable).
# p7.4 Referencia restringida: modelo ajustado excluyendo los promotores del
#   decil inferior de ocupacion (y, mas estricto, del cuartil inferior).
# GSEA de MLL1/2 complex (mismos parametros que remap_ruido_ajustado.R) sobre
#   los rankings nuevos de p7.3 y p7.4.
#
# Salida en figures/exploratorio_remap_ruido_ajustado/ con prefijo p7 (solo
# jpg). Requiere data/processed/{remap_tf_hits,activity_stats_highconf,
# prom_df}.tsv, remap_tf_noise_auc.tsv y p1_coeficientes_por_factor.tsv
# (remap_ruido_ajustado.R). Tarda unos minutos (GSEA). Run from the TesisDoc
# repo root.

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter

BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/exploratorio_remap_ruido_ajustado"
hs <- org.Hs.eg.db

cov_window <- 0.25
mll_id <- "GO:0044665"
compass <- c("KMT2A", "KMT2B", "WDR5", "RBBP5", "ASH2L")

# --- Datos ---------------------------------------------------------------

auc_orig <- read_tsv(file.path(fig_dir("remap_tf_noise_gsea"), "remap_tf_noise_auc.tsv"), show_col_types = FALSE)
tfs <- unique(auc_orig$feature)
coef_p1 <- read_tsv(file.path(out_dir, "p1_coeficientes_por_factor.tsv"), show_col_types = FALSE)
mll_tfs <- read_tsv(file.path(out_dir, "p4_mll1_2_miembros.tsv"), show_col_types = FALSE) %>% pull(TF) %>% unique()
stopifnot(all(compass %in% mll_tfs))

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE) %>% distinct(seq_id, TF)
ocupacion <- remap_hits %>% count(seq_id, name = "n_full")
sets <- remap_hits %>% filter(TF %in% tfs) %>% with(split(seq_id, TF))

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  left_join(ocupacion, by = "seq_id") %>%
  mutate(n_full = replace_na(n_full, 0L), log_hek = log10(hek_tpm + 1), log_n = log2(n_full + 1)) %>%
  select(seq_id, rep, var_rank_sw, log_hek, mean_dnase, n_full, log_n)
stopifnot(!anyNA(data))
repname <- unique(data$rep)

# --- p7.3 / p7.4: modelos nuevos para los 982 factores ----------------------

# cortes de ocupacion por replica (sobre los promotores de esa replica); se
# excluyen los promotores con n_full <= corte
cortes <- data %>%
  group_by(rep) %>%
  summarise(decil = quantile(n_full, 0.1), cuartil = quantile(n_full, 0.25))

variantes <- list(
  ajustado_ocupacion = list(fo = var_rank_sw ~ unido + log_hek + mean_dnase + log_n, corte = NULL),
  ajustado_sin_decil = list(fo = var_rank_sw ~ unido + log_hek + mean_dnase, corte = "decil"),
  ajustado_sin_cuartil = list(fo = var_rank_sw ~ unido + log_hek + mean_dnase, corte = "cuartil")
)

n_promotores <- list()
coef_nuevos <- map(repname, function(r) {
  d_rep <- data %>% filter(rep == r)
  imap(variantes, function(v, m) {
    d <- d_rep
    if (!is.null(v$corte)) d <- d %>% filter(n_full > cortes[[v$corte]][cortes$rep == r])
    n_promotores[[length(n_promotores) + 1]] <<- tibble(rep = r, modelo = m, n_promotores = nrow(d), n_full_min = min(d$n_full))
    map(tfs, function(tf) {
      d$unido <- as.numeric(d$seq_id %in% sets[[tf]])
      fit <- lm(v$fo, data = d)
      co <- summary(fit)$coefficients
      ci <- confint(fit)["unido", ]
      tibble(
        TF = tf, rep = r, modelo = m, n_unidos = sum(d$unido), n_no_unidos = sum(d$unido == 0),
        estimate = co["unido", 1], ci_low = ci[1], ci_high = ci[2], pval = co["unido", 4],
        coef_log_n = if ("log_n" %in% rownames(co)) co["log_n", 1] else NA_real_
      )
    }) %>% list_rbind()
  }) %>% list_rbind()
}) %>%
  list_rbind() %>%
  group_by(rep, modelo) %>%
  mutate(pval_corr = p.adjust(pval, method = "BH")) %>%
  ungroup()
n_promotores <- bind_rows(n_promotores) %>% left_join(cortes, by = "rep")
write_tsv(n_promotores, file.path(out_dir, "p7_4_n_promotores.tsv"))
print(n_promotores)

# todos los modelos en formato largo (los de p1 + los nuevos)
coef_all <- bind_rows(
  coef_p1 %>%
    select(TF, rep, n_unidos, matches("^(estimate|ci_low|ci_high|pval_corr)_(crudo|solo_actividad|ajustado)$")) %>%
    pivot_longer(-c(TF, rep, n_unidos), names_to = c(".value", "modelo"), names_pattern = "^(estimate|ci_low|ci_high|pval_corr)_(.*)$"),
  coef_nuevos %>% select(TF, rep, modelo, n_unidos, n_no_unidos, estimate, ci_low, ci_high, pval_corr, coef_log_n)
)
write_tsv(coef_nuevos, file.path(out_dir, "p7_3_4_coeficientes_por_factor.tsv"))

modelo_lab <- c(
  crudo = "Crudo", solo_actividad = "Ajustado por actividad endógena",
  ajustado = "Ajustado por actividad endógena y accesibilidad",
  ajustado_ocupacion = "Ajustado + ocupación total (log2)",
  ajustado_sin_decil = "Ajustado, sin el decil inferior de ocupación",
  ajustado_sin_cuartil = "Ajustado, sin el cuartil inferior de ocupación"
)

# resumen por modelo: mediana de los 982, consistentes y dependencia de la cobertura
resumen <- coef_all %>%
  group_by(modelo, rep) %>%
  summarise(
    mediana_982 = median(estimate),
    rho_n_unidos_coef = cor(n_unidos, estimate, method = "spearman"),
    mediana_mll = median(estimate[TF %in% mll_tfs]),
    mediana_compass = median(estimate[TF %in% compass]),
    coef_log_n_mediano = if (all(is.na(coef_log_n))) NA_real_ else median(coef_log_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    coef_all %>%
      group_by(modelo, TF) %>%
      summarise(cons = all(pval_corr < 0.05) & n_distinct(sign(estimate)) == 1, signo = mean(estimate) < 0, .groups = "drop") %>%
      group_by(modelo) %>%
      summarise(consistentes_menor_ruido = sum(cons & signo), consistentes_mayor_ruido = sum(cons & !signo)),
    by = "modelo"
  ) %>%
  mutate(modelo = factor(modelo, levels = names(modelo_lab))) %>%
  arrange(modelo, rep)
write_tsv(resumen, file.path(out_dir, "p7_3_resumen_modelos.tsv"))
print(as.data.frame(resumen), digits = 3)

# --- p7.1: ventana de cobertura ---------------------------------------------

pct_of <- function(x, ref) 100 * mean(ref <= x)
rank_of <- function(x, ref) sum(ref < x) + 1L

posicion <- coef_all %>%
  filter(TF %in% mll_tfs) %>%
  select(modelo, rep, TF, n_unidos, estimate) %>%
  pmap(function(modelo, rep, TF, n_unidos, estimate) {
    ref <- coef_all %>% filter(modelo == !!modelo, rep == !!rep, TF != !!TF)
    win <- ref %>% filter(abs(n_unidos - !!n_unidos) <= cov_window * !!n_unidos)
    win_sm <- win %>% filter(!TF %in% mll_tfs)
    tibble(
      modelo, rep, TF, nucleo_compass = TF %in% compass, n_unidos, estimate,
      n_universo = nrow(ref), rango_universo = rank_of(estimate, ref$estimate), percentil_universo = pct_of(estimate, ref$estimate),
      win_min = (1 - cov_window) * n_unidos, win_max = (1 + cov_window) * n_unidos,
      n_ventana = nrow(win), rango_ventana = rank_of(estimate, win$estimate), percentil_ventana = pct_of(estimate, win$estimate),
      mediana_ventana = median(win$estimate),
      n_ventana_sin_mll = nrow(win_sm), rango_ventana_sin_mll = rank_of(estimate, win_sm$estimate), percentil_ventana_sin_mll = pct_of(estimate, win_sm$estimate)
    )
  }) %>%
  list_rbind() %>%
  mutate(modelo = factor(modelo, levels = names(modelo_lab))) %>%
  arrange(modelo, desc(nucleo_compass), TF, rep)
write_tsv(posicion, file.path(out_dir, "p7_1_ventana_cobertura_mll.tsv"))
print(posicion %>% filter(nucleo_compass) %>% select(modelo, rep, TF, n_unidos, estimate, rango_universo, n_universo, percentil_universo, rango_ventana, n_ventana, percentil_ventana, percentil_ventana_sin_mll) %>% as.data.frame(), digits = 3)

# --- p7.2: figura -------------------------------------------------------------

compass_clr <- "#AD343E"
foc <- posicion %>%
  filter(modelo == "ajustado", nucleo_compass) %>%
  mutate(
    TF = factor(TF, levels = compass),
    label = sprintf("p%.0f universo (%d/%d)\np%.0f ventana (%d/%d)", percentil_universo, rango_universo, n_universo + 1L, percentil_ventana, rango_ventana, n_ventana + 1L)
  )
ref_plot <- coef_all %>%
  filter(modelo == "ajustado") %>%
  select(rep, TF_ref = TF, n_unidos, estimate) %>%
  crossing(TF = factor(compass, levels = compass)) %>%
  filter(TF_ref != TF)
p_cov <- ggplot(ref_plot, aes(n_unidos, estimate)) +
  geom_rect(data = foc, aes(xmin = win_min, xmax = win_max, ymin = -Inf, ymax = Inf), fill = compass_clr, alpha = 0.15, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_point(color = "grey60", size = 0.9, alpha = 0.6) +
  geom_point(data = foc, shape = 21, color = compass_clr, fill = compass_clr, size = 3.2) +
  geom_label(data = foc, aes(x = min(ref_plot$n_unidos), y = Inf, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2, label.size = 0, fill = alpha("white", 0.85)) +
  facet_grid(rep ~ TF) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(
    x = "Promotores unidos (ReMap, escala log10)", y = "Efecto ajustado de la unión sobre\nel rango de varianza (1-100)",
    caption = "Punto rojo: miembro del núcleo COMPASS de cada columna. Sombreado: ventana de cobertura ±25%. Rango 1 = mayor reducción del ruido; posición / n de factores comparados."
  ) +
  ggpubr::theme_pubr(base_size = 14) +
  theme(plot.caption = element_text(size = 10), strip.background = element_blank())
ggsave(file.path(out_dir, "p7_2_efecto_vs_cobertura_compass.jpg"), p_cov, width = 20, height = 8.5, units = "in")

# p7.3/p7.4: coeficiente de los miembros de MLL1/2 en los modelos de control
# mismo orden que p4; nucleo COMPASS en negrita
mll_orden <- coef_p1 %>% filter(TF %in% mll_tfs, rep == repname[1]) %>% arrange(estimate_crudo) %>% pull(TF)
p_mll <- coef_all %>%
  filter(TF %in% mll_tfs, modelo %in% c("ajustado", names(variantes))) %>%
  mutate(
    modelo = factor(modelo_lab[modelo], levels = modelo_lab[c("ajustado", names(variantes))]),
    TF = factor(TF, levels = mll_orden),
    sig = pval_corr < 0.05
  ) %>%
  ggplot(aes(estimate, TF, col = modelo, shape = sig)) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, position = position_dodge(0.75)) +
  geom_point(size = 2.3, position = position_dodge(0.75)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p corregido < 0.05", `FALSE` = "p corregido ≥ 0.05")) +
  scale_color_manual(values = c("#AD343E", "#1B8C8E", "#D6741F", "#0D2C54")) +
  facet_wrap(~rep) +
  ggpubr::theme_pubr(base_size = 13) +
  theme(legend.position = "top", legend.direction = "vertical", strip.background = element_blank(),
    axis.text.y = element_text(face = ifelse(mll_orden %in% compass, "bold", "plain"))) +
  labs(x = "Efecto de la unión sobre el rango de varianza (IC 95%)", y = "Miembros de MLL1/2 complex en ReMap", col = NULL, shape = NULL)
ggsave(file.path(out_dir, "p7_3_4_mll1_2_miembros.jpg"), p_mll, width = 13, height = 10, units = "in")

# --- GSEA de MLL1/2 sobre los rankings nuevos -------------------------------

ids <- AnnotationDbi::select(hs, keys = tfs, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% distinct()
run_gsea <- function(m, r) {
  pre <- coef_nuevos %>%
    filter(modelo == m, rep == r) %>%
    arrange(desc(estimate)) %>%
    left_join(ids, by = c("TF" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  gseGO(setNames(pre$estimate, pre$ENTREZID), ont = "CC", OrgDb = "org.Hs.eg.db", pvalueCutoff = 1, seed = TRUE)
}
gsea_mll <- map(names(variantes), function(m) {
  set.seed(2024)
  map(repname, function(r) {
    as.data.frame(run_gsea(m, r)) %>%
      filter(ID == mll_id) %>%
      transmute(ranking = m, rep = r, ID, Description, setSize, NES, pvalue, p.adjust)
  }) %>% list_rbind()
}) %>% list_rbind()
# referencia: los rankings de p3
gsea_p3 <- read_tsv(file.path(out_dir, "p3_gsea_todos_los_terminos.tsv"), show_col_types = FALSE) %>%
  filter(ID == mll_id) %>%
  select(ranking, rep, ID, Description, setSize, NES, pvalue, p.adjust)
gsea_mll <- bind_rows(gsea_p3, gsea_mll)
write_tsv(gsea_mll, file.path(out_dir, "p7_gsea_mll1_2.tsv"))
print(as.data.frame(gsea_mll), digits = 3)

message("Exploratorio guardado en ", out_dir)
