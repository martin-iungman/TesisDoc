# ruido_contexto_endogeno (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Figura que cierra el bloque de contexto endogeno de la seccion de ruido.
# Fusiona remap_tf_noise_gsea y chipatlas_histonas y suma el gradiente de
# ocupacion de figures/exploratorio_ausencia_senal/ (panel_ocupacion):
#   A  AUC de la union de cada factor de ReMap sobre el ruido (los
#      consistentes, formato original de remap_tf_noise_gsea).
#   B  GSEA (GO CC) de MLL1/2 complex en Rep 1 (formato original).
#   C  AUC de las marcas de histona de ChIP-Atlas consistentes (formato
#      original de chipatlas_histonas), H3K4me3 resaltada, con la
#      prevalencia de cada marca en una columna adjunta.
#   D  Rango de varianza vs. bin de actividad segun H3K4me3, solo
#      promotores sin isla CpG (CGI_frag), bins y rangos recalculados
#      dentro del subconjunto.
#   E  Rango de varianza medio por numero de factores de ReMap distintos
#      (bin 0 + deciles del resto), por isla CpG.
# Replicas por alpha en A-C (formato original) y en facetas en D-E.
#
# Los calculos reproducen los de R/21_remap_noise_gsea (mismos criterios,
# misma semilla del GSEA), R/22_chipatlas_histonas (+ p de Wald DeLong y
# BH, como ruido_summary) y R/54_exploratorio_ausencia_senal (gradiente de
# ocupacion; aca sobre el set completo, sin excluir los 4 / 5 promotores
# sin clase de FANTOM5, que ese analisis excluia por usar FANTOM5).
# Ademas guarda en numeros_texto.tsv los valores que cita el texto.
#
# Requiere: data/processed/{activity_stats_highconf,prom_df,remap_tf_hits}.tsv,
# data/external/allPeaks_chipatlas_counted.tsv y la experimentList de
# ChIP-Atlas (path_chipatlas_explist, R/00_prom_features/heavy_data_paths.R).
# ~2 min (AUC de ~1000 factores y GSEA). Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(pROC)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
BiocParallel::register(BiocParallel::SerialParam())

slug <- "ruido_contexto_endogeno"
out_dir <- fig_dir(slug)

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}
num <- function(d = 2) function(x) formatC(x, format = "f", digits = d, big.mark = "")
tema <- theme_pubr(base_size = 16) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, face = "bold"), legend.position = "top", panel.spacing = unit(1.5, "lines"))
col_efecto <- c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800", "n.s." = "grey70")
col_estrato <- c("Con isla CpG" = "#1B8C8E", "Sin isla CpG" = "#D6741F")
col_set <- "darkred"

# AUC (P(ruido con > ruido sin)), IC y p de Wald con la varianza de DeLong,
# igual que ruido_summary. `x` logico.
auc_bool <- function(x, score) {
  resp <- factor(x, levels = c("TRUE", "FALSE"))
  roc_obj <- pROC::roc(resp, score, ci = TRUE, direction = ">", quiet = TRUE)
  auc_val <- as.numeric(roc_obj$ci[2])
  se <- sqrt(pROC::var(roc_obj, method = "delong"))
  tibble(n_con = sum(x), n_sin = sum(!x), AUC = auc_val, ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]), pval = 2 * pnorm(-abs((auc_val - 0.5) / se)))
}
# significativo (BH < 0.05) y mismo sentido en ambas replicas
estado_consistente <- function(df) {
  df %>%
    group_by(feature) %>%
    mutate(efecto = if (all(pval_corr < 0.05) && n() == 2 && n_distinct(AUC > 0.5) == 1) if_else(AUC > 0.5, "Ruido alto", "Ruido bajo") else "n.s.") %>%
    ungroup()
}

# --- Datos --------------------------------------------------------------------
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>% add_noise_rank()
repname <- sort(unique(data$rep))
n_rep <- count(data, rep, name = "N")

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE) %>% distinct(seq_id, TF)

chipatlas <- read_tsv("data/external/allPeaks_chipatlas_counted.tsv", col_names = c("seq_id", "feature", "n_samples"), show_col_types = FALSE)
explist <- read_tsv(path_chipatlas_explist, col_names = FALSE, show_col_types = FALSE) %>%
  select(X3, X4) %>%
  distinct(X4, .keep_all = TRUE)
histones <- chipatlas %>%
  left_join(explist, by = c("feature" = "X4")) %>%
  filter(X3 == "Histone") %>%
  distinct(seq_id, feature)

# --- A. ReMap: AUC por factor ---------------------------------------------------
# Factores con >100 promotores unidos en ambas replicas, sin marcas de
# histona (mismo criterio que remap_tf_noise_gsea).
tf_n <- data %>%
  select(seq_id, rep) %>%
  inner_join(remap_hits, by = "seq_id", relationship = "many-to-many") %>%
  count(rep, TF)
vbles_tf <- tf_n %>%
  filter(n > 100, !str_detect(TF, "^H\\d")) %>%
  count(TF) %>%
  filter(n == 2) %>%
  pull(TF)
hits_por_tf <- split(remap_hits$seq_id, remap_hits$TF)

auc_tf <- map(repname, function(r) {
  df <- filter(data, rep == r)
  map(vbles_tf, ~ auc_bool(df$seq_id %in% hits_por_tf[[.x]], df$var_rank_sw) %>% mutate(feature = .x)) %>%
    list_rbind() %>%
    mutate(rep = r)
}) %>%
  list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, method = "BH")) %>%
  estado_consistente() %>%
  relocate(feature, rep)
write_tsv(auc_tf, file.path(out_dir, "panel_a_remap_auc.tsv"))

# --- B. GSEA (GO CC) sobre la lista de factores ordenada por AUC ---------------
# Igual que remap_tf_noise_gsea (misma semilla).
ids <- AnnotationDbi::select(org.Hs.eg.db, keys = vbles_tf, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% distinct()
pregsea <- function(r) {
  auc_tf %>%
    filter(rep == r) %>%
    arrange(desc(AUC)) %>%
    left_join(ids, by = c("feature" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
}
set.seed(2024)
gse <- map(repname, function(r) {
  p <- pregsea(r)
  gseGO(setNames(p$AUC, p$ENTREZID), ont = "CC", OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.1, seed = TRUE)
})
names(gse) <- repname
mll_id <- "GO:0044665"
gsea_mll <- map(repname, ~ as.data.frame(gse[[.x]]) %>% filter(ID == mll_id) %>% mutate(rep = .x)) %>% list_rbind()
stopifnot(nrow(gsea_mll) == 2)
write_tsv(gsea_mll %>% select(rep, ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, core_enrichment), file.path(out_dir, "panel_b_gsea_mll.tsv"))
mll_tf <- pregsea(repname[1]) %>% filter(ENTREZID %in% gse[[1]]@geneSets[[mll_id]]) %>% pull(feature)


# --- C. ChIP-Atlas: AUC por marca + prevalencia --------------------------------
hist_por_marca <- split(histones$seq_id, histones$feature)
marca_n <- data %>%
  select(seq_id, rep) %>%
  inner_join(histones, by = "seq_id", relationship = "many-to-many") %>%
  count(rep, feature)
vbles_hist <- marca_n %>% filter(n > 100) %>% count(feature) %>% filter(n == 2) %>% pull(feature)

auc_hist <- map(repname, function(r) {
  df <- filter(data, rep == r)
  map(vbles_hist, ~ auc_bool(df$seq_id %in% hist_por_marca[[.x]], df$var_rank_sw) %>% mutate(feature = .x)) %>%
    list_rbind() %>%
    mutate(rep = r)
}) %>%
  list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, method = "BH"), prevalencia = n_con / (n_con + n_sin)) %>%
  estado_consistente() %>%
  relocate(feature, rep)
write_tsv(auc_hist, file.path(out_dir, "panel_c_chipatlas_auc.tsv"))

# --- D / numeros: H3K4me3 por isla CpG ------------------------------------------
# Version principal (la del panel D): bins y rango recalculados dentro de
# cada estrato. Sensibilidad: rango global restringido al estrato.
h3k4me3 <- hist_por_marca[["H3K4me3"]]
data_h3 <- data %>% mutate(H3K4me3 = seq_id %in% h3k4me3, estrato = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG"))
data_h3_est <- data_h3 %>%
  select(seq_id, rep, mean, var, estrato, H3K4me3) %>%
  group_split(estrato) %>%
  map(add_noise_rank) %>%
  list_rbind()
auc_h3 <- bind_rows(
  data_h3 %>% group_by(rep) %>% group_modify(~ auc_bool(.x$H3K4me3, .x$var_rank_sw)) %>% mutate(estrato = "Todos", rango = "global"),
  data_h3_est %>% group_by(rep, estrato) %>% group_modify(~ auc_bool(.x$H3K4me3, .x$var_rank_sw)) %>% mutate(rango = "recalculado en el estrato"),
  data_h3 %>% group_by(rep, estrato) %>% group_modify(~ auc_bool(.x$H3K4me3, .x$var_rank_sw)) %>% mutate(rango = "global restringido al estrato")
) %>%
  ungroup() %>%
  mutate(prevalencia = n_con / (n_con + n_sin)) %>%
  relocate(estrato, rango, rep)
write_tsv(auc_h3, file.path(out_dir, "h3k4me3_auc_por_isla.tsv"))

# --- E. Gradiente de ocupacion ---------------------------------------------------
ocup <- data %>%
  left_join(remap_hits %>% count(seq_id, name = "n_full"), by = "seq_id") %>%
  mutate(n_full = replace_na(n_full, 0L), estrato = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")) %>%
  group_by(rep) %>%
  mutate(decil = if_else(n_full == 0, 0L, as.integer(cut(n_full, unique(quantile(n_full[n_full > 0], 0:10 / 10)), include.lowest = TRUE, labels = FALSE)))) %>%
  ungroup() %>%
  mutate(bin = factor(if_else(decil == 0, "0", paste0("D", decil)), c("0", paste0("D", 1:10))))
resumen_ic <- function(df) {
  df %>% summarise(
    n = n(), media = mean(var_rank_sw), se = sd(var_rank_sw) / sqrt(n()),
    ic_inf = media - qt(0.975, pmax(n() - 1, 1)) * se, ic_sup = media + qt(0.975, pmax(n() - 1, 1)) * se,
    ocup_min = min(n_full), ocup_max = max(n_full), .groups = "drop"
  )
}
grad <- bind_rows(
  ocup %>% group_by(rep, bin) %>% resumen_ic() %>% mutate(estrato = "Todos"),
  ocup %>% group_by(rep, estrato, bin) %>% resumen_ic()
) %>% relocate(rep, estrato, bin)
rho <- ocup %>%
  bind_rows(mutate(ocup, estrato = "Todos")) %>%
  group_by(rep, estrato) %>%
  summarise(rho = cor(n_full, var_rank_sw, method = "spearman"), p = cor.test(n_full, var_rank_sw, method = "spearman", exact = FALSE)$p.value, .groups = "drop")
write_tsv(grad %>% left_join(rho, by = c("rep", "estrato")), file.path(out_dir, "panel_e_ocupacion.tsv"))

# Meseta: medias por tramos fijos de ocupacion y recta quebrada (quiebre
# por grilla de 50 a 600, sin bootstrap - el IC esta en
# exploratorio_ausencia_senal/p2d).
tramos <- ocup %>%
  mutate(tramo = cut(n_full, c(-1, 0, 100, 200, 300, 400, 500, 600, Inf), labels = c("0", "1-100", "101-200", "201-300", "301-400", "401-500", "501-600", ">600"))) %>%
  group_by(rep, tramo) %>%
  resumen_ic()
quiebre <- ocup %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    grilla <- seq(50, 600, 25)
    rss <- map_dbl(grilla, ~ deviance(lm(var_rank_sw ~ n_full + pmax(n_full - .x, 0), data = d)))
    q <- grilla[which.min(rss)]
    fit <- lm(var_rank_sw ~ n_full + pmax(n_full - q, 0), data = d)
    tibble(quiebre = q, pendiente_antes_x100 = 100 * coef(fit)[2], pendiente_despues_x100 = 100 * (coef(fit)[2] + coef(fit)[3]))
  }) %>%
  ungroup()
write_tsv(tramos, file.path(out_dir, "meseta_tramos_fijos.tsv"))
write_tsv(quiebre, file.path(out_dir, "meseta_quiebre.tsv"))

# --- Numeros que cita el texto ---------------------------------------------------
cons_tf <- auc_tf %>% distinct(feature, efecto)
mll <- gsea_mll %>% select(rep, NES, p.adjust)
h3_glob <- auc_hist %>% filter(feature == "H3K4me3")
orden_hist <- auc_hist %>% group_by(rep) %>% arrange(AUC, .by_group = TRUE) %>% slice_head(n = 5) %>% summarise(top5 = paste(sprintf("%s (%.3f)", feature, AUC), collapse = "; "))
h3_est <- auc_h3 %>% filter(rango == "recalculado en el estrato")
cero <- grad %>% filter(estrato == "Todos", bin == "0")
txt <- function(...) tibble(...)
numeros <- bind_rows(
  txt(item = "1. Marcas de histona evaluadas (>100 promotores con la marca en ambas réplicas)", rep = "ambas", valor = as.character(length(vbles_hist))),
  txt(item = "1b. Marcas con efecto consistente (BH < 0.05, mismo sentido)", rep = "ambas", valor = paste(count(distinct(auc_hist, feature, efecto), efecto) %>% transmute(s = paste0(efecto, ": ", n)) %>% pull(s), collapse = "; ")),
  map_dfr(repname, function(r) {
    a <- filter(h3_glob, rep == r)
    txt(item = "2. Prevalencia de H3K4me3 en la library", rep = r, valor = sprintf("%.1f%% (%d de %d)", 100 * a$prevalencia, a$n_con, a$n_con + a$n_sin))
  }),
  map_dfr(c("Sin isla CpG", "Con isla CpG"), function(e) map_dfr(repname, function(r) {
    a <- filter(h3_est, rep == r, estrato == e)
    txt(item = paste0(if_else(e == "Sin isla CpG", "3", "4"), ". AUC de H3K4me3, ", e, " (rango recalculado en el estrato)"), rep = r,
      valor = sprintf("%.3f [%.3f, %.3f], p %.2g; con H3K4me3 %d, sin %d (%.1f%% con)", a$AUC, a$ci2.5, a$ci97.5, a$pval, a$n_con, a$n_sin, 100 * a$prevalencia))
  })),
  txt(item = "5a. Factores de ReMap evaluados", rep = "ambas", valor = as.character(length(vbles_tf))),
  txt(item = "5a. Factores con efecto consistente", rep = "ambas", valor = paste(count(cons_tf, efecto) %>% transmute(s = paste0(efecto, ": ", n)) %>% pull(s), collapse = "; ")),
  map_dfr(repname, function(r) txt(item = "5b. MLL1/2 complex (GO:0044665): NES, p ajustado", rep = r, valor = sprintf("NES %.2f, p ajustado %.3f", mll$NES[mll$rep == r], mll$p.adjust[mll$rep == r]))),
  map_dfr(repname, function(r) txt(item = "5c. Cinco marcas de menor AUC", rep = r, valor = orden_hist$top5[orden_hist$rep == r])),
  map_dfr(repname, function(r) {
    a <- filter(cero, rep == r)
    txt(item = "5d. Rango de varianza medio en promotores sin ningún pico de ReMap", rep = r, valor = sprintf("%.1f [%.1f, %.1f], n %d; media del conjunto %.1f", a$media, a$ic_inf, a$ic_sup, a$n, mean(data$var_rank_sw[data$rep == r])))
  }),
  map_dfr(repname, function(r) {
    q <- filter(quiebre, rep == r)
    tr <- filter(tramos, rep == r)
    txt(item = "5e. Meseta: quiebre de la recta quebrada y medias por tramo", rep = r, valor = sprintf("quiebre %d factores; tramos: %s", q$quiebre, paste(sprintf("%s %.1f", tr$tramo, tr$media), collapse = ", ")))
  })
)
write_tsv(numeros, file.path(out_dir, "numeros_texto.tsv"))
print(numeros, n = Inf, width = Inf)

# --- Paneles ----------------------------------------------------------------------
# A, B y C en el formato de las figuras originales (remap_tf_noise_gsea,
# chipatlas_histonas): barras de AUC - 0.5 de plot_noise_auc_summary(),
# solo lo significativo y del mismo sentido en ambas replicas, replicas por
# alpha; D y E con replicas en facetas. Mismo base_size en todos.
col_destacado <- "#AD343E"

# A: factores de ReMap consistentes (sin etiquetas ni IC: son 514)
panel_a <- plot_noise_auc_summary(auc_tf, show_labels = FALSE, show_errorbars = FALSE, base_size = 16) +
  labs(y = "Factores de transcripción (ReMap)")
save_jpg_pdf(panel_a, "panel_a_remap_auc", 9, 6.75)

# B: GSEA de MLL1/2 complex en Rep 1 (curva + todos los factores de la lista
# en el mismo orden, los del complejo resaltados)
mll_desc <- gsea_mll$Description[1]
pregsea_rep1 <- pregsea(repname[1])
p_b1 <- enrichplot::gseaplot2(gse[[repname[1]]], mll_id, title = mll_desc, subplots = 1, color = col_set, base_size = 16) +
  labs(y = "Puntaje de\nenriquecimiento")
p_b2 <- pregsea_rep1 %>%
  mutate(noise = case_when(feature %in% mll_tf ~ "set", AUC > 0.5 ~ "Ruido alto", TRUE ~ "Ruido bajo")) %>%
  ggplot(aes(x = AUC - 0.5, y = fct_inorder(feature))) +
  geom_col(orientation = "y", aes(fill = noise, alpha = noise)) +
  scale_x_continuous(labels = function(x) x + 0.5) +
  theme_pubr(base_size = 16) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "AUC (efecto\nsobre el ruido)", y = "Factores de transcripción") +
  scale_fill_manual(values = c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800", "set" = col_set)) +
  scale_alpha_manual(values = c("Ruido alto" = 0.4, "Ruido bajo" = 0.4, "set" = 1)) +
  coord_flip()
panel_b <- ggplotify::as.ggplot(aplot::gglist(list(p_b1, p_b2), ncol = 1, heights = c(0.75, 0.5)))
save_jpg_pdf(panel_b, "panel_b_gsea_mll", 10, 7)

# C: marcas consistentes, H3K4me3 resaltada en color; columna adjunta con
# el % de promotores con la marca (Rep 1 / Rep 2), en el mismo orden
destacar <- function(x) if_else(x == "H3K4me3", paste0("<b style='color:", col_destacado, "'>", x, "</b>"), x)
panel_c_auc <- plot_noise_auc_summary(auc_hist, base_size = 16) +
  scale_y_discrete(labels = destacar) +
  theme(axis.text.y = ggtext::element_markdown()) +
  labs(y = "Marcas de histona (ChIP-Atlas)")
niveles_c <- levels(panel_c_auc$data$feature)
d_c_prev <- auc_hist %>%
  filter(feature %in% niveles_c) %>%
  group_by(feature) %>%
  summarise(lab = paste(sprintf("%.1f", 100 * prevalencia[order(rep)]), collapse = " / "), .groups = "drop") %>%
  mutate(lab = if_else(feature == "H3K4me3", paste0("<b style='color:", col_destacado, "'>", lab, "</b>"), lab), feature = factor(feature, niveles_c))
panel_c_prev <- ggplot(d_c_prev, aes(x = 0, y = feature, label = lab)) +
  ggtext::geom_richtext(size = 4.5, fill = NA, label.color = NA) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "% con la marca\n(Rep 1 / Rep 2)", y = NULL) +
  theme_pubr(base_size = 16) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
panel_c <- (panel_c_auc | panel_c_prev) + plot_layout(widths = c(4, 1.2), guides = "keep")
save_jpg_pdf(panel_c, "panel_c_chipatlas_auc", 10, 2.5 + 0.3 * length(niveles_c))

# D: scatter de ruido vs. actividad segun H3K4me3, sin isla CpG
d_d <- data_h3_est %>% filter(estrato == "Sin isla CpG") %>% mutate(H3K4me3 = factor(if_else(H3K4me3, "Con H3K4me3", "Sin H3K4me3"), c("Sin H3K4me3", "Con H3K4me3")))
lab_d <- h3_est %>%
  filter(estrato == "Sin isla CpG") %>%
  transmute(rep, lab = sprintf("Con: n = %d\nSin: n = %d\nAUC = %.2f", n_con, n_sin, AUC))
panel_d <- ggplot(d_d, aes(mean_sw, var_rank_sw)) +
  geom_point(aes(col = H3K4me3), size = 0.2, alpha = 0.5) +
  scale_color_manual(values = c("Sin H3K4me3" = "grey65", "Con H3K4me3" = thesis_clr), name = NULL, guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  ggnewscale::new_scale_color() +
  geom_smooth(aes(col = H3K4me3), se = TRUE, linewidth = 1.2, show.legend = FALSE, method = "loess", formula = y ~ x) +
  scale_color_manual(values = c("Sin H3K4me3" = "#525252", "Con H3K4me3" = "#136869")) +
  geom_label(data = lab_d, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.05, size = 4, label.size = 0, fill = alpha("white", 0.8), inherit.aes = FALSE) +
  facet_wrap(~rep, scales = "free_x") +
  scale_x_continuous(labels = num(0)) +
  coord_cartesian(ylim = c(1, 100)) +
  labs(x = "Bin de actividad", y = "Rango de varianza") +
  tema
save_jpg_pdf(panel_d, "panel_d_h3k4me3_sin_isla", 11, 5.5)

# E: gradiente de ocupacion por isla CpG, n debajo de cada punto
d_e <- grad %>% filter(estrato != "Todos") %>% mutate(estrato = factor(estrato, names(col_estrato)), x = as.integer(bin) - 1)
ye <- range(c(d_e$ic_inf, d_e$ic_sup, 50))
paso <- diff(ye) * 0.09
d_e <- d_e %>% mutate(y_n = if_else(estrato == "Con isla CpG", ye[1] - paso, ye[1] - 2 * paso))
rho_lab <- rho %>%
  filter(estrato != "Todos") %>%
  mutate(estrato = factor(estrato, names(col_estrato)), label = sprintf('rho == "%s"', formatC(rho, format = "f", digits = 2)), y = if_else(estrato == "Con isla CpG", ye[2] + 1.4 * paso, ye[2] + 0.5 * paso))
panel_e <- ggplot(d_e, aes(x, media, color = estrato)) +
  geom_hline(yintercept = 50.5, linetype = "dotted", col = "grey60") +
  geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0, position = position_dodge(width = 0.45), linewidth = 0.8) +
  geom_line(position = position_dodge(width = 0.45), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.45), size = 2.5) +
  geom_text(aes(y = y_n, label = n), size = 3, show.legend = FALSE) +
  annotate("text", x = -0.75, y = ye[1] - 1.5 * paso, label = "n", size = 4, fontface = "italic") +
  geom_text(data = rho_lab, aes(x = 10.4, y = y, label = label), parse = TRUE, hjust = 1, size = 4.5, show.legend = FALSE) +
  scale_color_manual(values = col_estrato, name = NULL) +
  scale_x_continuous(breaks = 0:10, labels = c("0", paste0("D", 1:10)), limits = c(-0.9, 10.5)) +
  scale_y_continuous(labels = num(0)) +
  coord_cartesian(ylim = c(ye[1] - 2.5 * paso, ye[2] + 2 * paso), clip = "off") +
  facet_wrap(~rep) +
  labs(x = "Factores de ReMap unidos (0 y deciles)", y = "Rango de varianza\nmedio") +
  tema +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_jpg_pdf(panel_e, "panel_e_ocupacion", 11, 5.5)

# --- Figura compuesta ---------------------------------------------------------------
# Orden de lectura: A | B, C | D, E abajo a lo ancho.
h <- c(7, 7.5, 6)
figura <- wrap_elements(panel_a) + wrap_elements(panel_b) + wrap_elements(panel_c) + wrap_elements(panel_d) + wrap_elements(panel_e) +
  plot_layout(design = "AB\nCD\nEE", heights = h)
figura <- cowplot::ggdraw(figura) +
  cowplot::draw_plot_label(c("A", "B", "C", "D", "E"), x = c(0, 0.5, 0, 0.5, 0), y = 1 - c(0, 0, h[1], h[1], h[1] + h[2]) / sum(h), size = 26)
save_jpg_pdf(figura, "ruido_contexto_endogeno", 20, sum(h))

message("Figuras guardadas en ", out_dir)
