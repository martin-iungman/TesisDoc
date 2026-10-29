# Exploratorio de remap_tf_noise_gsea (R/21_remap_noise_gsea): control por
# actividad endogena y accesibilidad de la cromatina.
#
# La union de un factor en ReMap marca en buena medida que el promotor esta
# activo y accesible en su contexto endogeno, y ambas cosas se asocian a
# menor ruido (ruido_summary). Aca, para cada uno de los 982 factores de
# remap_tf_noise_auc.tsv y por replica, tres modelos lineales sobre
# var_rank_sw (rango de varianza dentro del bin de actividad, 1-100):
#   crudo:         var_rank_sw ~ unido
#   solo_actividad: var_rank_sw ~ unido + log10(hek_tpm + 1)
#   ajustado:      var_rank_sw ~ unido + log10(hek_tpm + 1) + mean_dnase
# hek_tpm (CAGE en HEK293) y mean_dnase (DNase-seq) son las versiones
# continuas de prom_df (el panorama usa sus binarizaciones). No hay NAs en
# el set highconf, asi que los tres modelos usan los mismos promotores.
# BH dentro de cada replica y modelo. El modelo ajustado es conservador: la
# accesibilidad es en parte consecuencia de la union, no solo un confusor
# previo; solo_actividad ajusta por lo que no depende de la union local.
#
# GSEA (GO Cellular Component) con los mismos parametros que
# remap_tf_noise_gsea (gseGO, ont = "CC", set.seed(2024), seed = TRUE,
# parametros por defecto de fgsea) sobre cuatro rankings por replica: la AUC
# original (reproduce el GSEA actual) y el coeficiente de `unido` de cada
# modelo. pvalueCutoff = 1 para tener NES/p ajustado de todos los terminos
# (el p ajustado es el mismo que con 0.1; solo cambia que terminos se listan).
#
# Salida en figures/exploratorio_remap_ruido_ajustado/ (no es una figura
# numerada, solo jpg), con prefijo por punto del pedido (p1-p4) y
# handoff.md. Requiere data/processed/remap_tf_hits.tsv,
# activity_stats_highconf.tsv, prom_df.tsv y la tabla remap_tf_noise_auc.tsv
# de remap_tf_noise_gsea. Tarda unos minutos (GSEA). Run from the TesisDoc
# repo root.

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/exploratorio_remap_ruido_ajustado"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hs <- org.Hs.eg.db

modelos <- c(crudo = "Crudo", solo_actividad = "Ajustado por actividad endógena", ajustado = "Ajustado por actividad endógena y accesibilidad")

# --- Datos ---------------------------------------------------------------

auc_orig <- read_tsv(file.path(fig_dir("remap_tf_noise_gsea"), "remap_tf_noise_auc.tsv"), show_col_types = FALSE)
tfs <- unique(auc_orig$feature)

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE) %>%
  filter(TF %in% tfs) %>%
  distinct(seq_id, TF)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  mutate(log_hek = log10(hek_tpm + 1)) %>%
  select(seq_id, rep, var_rank_sw, log_hek, mean_dnase)
stopifnot(!anyNA(data))
repname <- unique(data$rep)

# correlacion de las covariables con el ruido (contexto para el handoff)
cov_cor <- data %>%
  group_by(rep) %>%
  summarise(
    rho_hek_ruido = cor(log_hek, var_rank_sw, method = "spearman"),
    rho_dnase_ruido = cor(mean_dnase, var_rank_sw, method = "spearman"),
    rho_hek_dnase = cor(log_hek, mean_dnase, method = "spearman"),
    n = n()
  )
write_tsv(cov_cor, file.path(out_dir, "p1_covariables_correlacion.tsv"))

sets <- split(remap_hits$seq_id, remap_hits$TF)

# --- p1: coeficiente de `unido` por factor, replica y modelo --------------

formulas <- list(
  crudo = var_rank_sw ~ unido,
  solo_actividad = var_rank_sw ~ unido + log_hek,
  ajustado = var_rank_sw ~ unido + log_hek + mean_dnase
)

coef_df <- map(repname, function(r) {
  d <- data %>% filter(rep == r)
  map(tfs, function(tf) {
    d$unido <- as.numeric(d$seq_id %in% sets[[tf]])
    imap(formulas, function(fo, m) {
      fit <- lm(fo, data = d)
      co <- summary(fit)$coefficients["unido", ]
      ci <- confint(fit)["unido", ]
      tibble(TF = tf, rep = r, modelo = m, n_unidos = sum(d$unido), estimate = co[1], ci_low = ci[1], ci_high = ci[2], pval = co[4])
    }) %>% list_rbind()
  }) %>% list_rbind()
}) %>%
  list_rbind() %>%
  group_by(rep, modelo) %>%
  mutate(pval_corr = p.adjust(pval, method = "BH")) %>%
  ungroup()

coef_wide <- coef_df %>%
  pivot_wider(id_cols = c(TF, rep, n_unidos), names_from = modelo, values_from = c(estimate, ci_low, ci_high, pval_corr), names_glue = "{.value}_{modelo}") %>%
  left_join(auc_orig %>% select(TF = feature, rep, AUC, pval_corr_AUC = pval_corr), by = c("TF", "rep"))
write_tsv(coef_wide, file.path(out_dir, "p1_coeficientes_por_factor.tsv"))

# --- p2: factores significativos y consistentes, reduccion del efecto ------

consistencia <- coef_df %>%
  group_by(TF, modelo) %>%
  summarise(
    consistente = all(pval_corr < 0.05) & n_distinct(sign(estimate)) == 1,
    signo = ifelse(mean(estimate) > 0, "mayor ruido", "menor ruido"),
    .groups = "drop"
  )
conteo <- consistencia %>%
  filter(consistente) %>%
  count(modelo, signo) %>%
  complete(modelo = names(modelos), signo = c("menor ruido", "mayor ruido"), fill = list(n = 0)) %>%
  bind_rows(
    auc_orig %>%
      filter(pval_corr < 0.05) %>%
      group_by(feature) %>%
      filter(n() == 2, n_distinct(AUC > 0.5) == 1) %>%
      summarise(signo = ifelse(mean(AUC) > 0.5, "mayor ruido", "menor ruido")) %>%
      count(signo) %>%
      complete(signo = c("menor ruido", "mayor ruido"), fill = list(n = 0)) %>%
      mutate(modelo = "AUC_original")
  )

# reduccion: entre los factores consistentes en el modelo crudo, mediana del
# coeficiente y mediana de 1 - ajustado / crudo
crudos_cons <- consistencia %>% filter(modelo == "crudo", consistente) %>% pull(TF)
reduccion <- coef_wide %>%
  filter(TF %in% crudos_cons) %>%
  group_by(rep) %>%
  summarise(
    n = n(),
    mediana_crudo = median(estimate_crudo),
    mediana_solo_actividad = median(estimate_solo_actividad),
    mediana_ajustado = median(estimate_ajustado),
    reduccion_mediana_solo_actividad = median(1 - estimate_solo_actividad / estimate_crudo),
    reduccion_mediana_ajustado = median(1 - estimate_ajustado / estimate_crudo)
  )
write_tsv(conteo, file.path(out_dir, "p2_conteo_consistentes.tsv"))
write_tsv(reduccion, file.path(out_dir, "p2_reduccion_efecto.tsv"))
print(conteo)
print(reduccion)

# --- p3: GSEA sobre cada ranking -------------------------------------------

ids <- AnnotationDbi::select(hs, keys = tfs, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% distinct()
rankings <- bind_rows(
  auc_orig %>% transmute(TF = feature, rep, ranking = "AUC_original", valor = AUC),
  coef_df %>% transmute(TF, rep, ranking = modelo, valor = estimate)
)

run_gsea <- function(rk, r) {
  pre <- rankings %>%
    filter(ranking == rk, rep == r) %>%
    arrange(desc(valor)) %>%
    left_join(ids, by = c("TF" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  # mismo orden de llamadas que remap_tf_noise_gsea: set.seed(2024) antes de
  # la tanda de replicas
  gseGO(setNames(pre$valor, pre$ENTREZID), ont = "CC", OrgDb = "org.Hs.eg.db", pvalueCutoff = 1, seed = TRUE)
}
gsea <- map(set_names(unique(rankings$ranking)), function(rk) {
  set.seed(2024)
  map(set_names(repname), ~ run_gsea(rk, .x))
})
gsea_df <- imap(gsea, function(g, rk) imap(g, ~ as.data.frame(.x) %>% mutate(ranking = rk, rep = .y)) %>% list_rbind()) %>%
  list_rbind() %>%
  as_tibble()
write_tsv(gsea_df %>% select(ranking, rep, ID, Description, setSize, NES, pvalue, p.adjust, core_enrichment), file.path(out_dir, "p3_gsea_todos_los_terminos.tsv"))

terminos <- c(
  "GO:0044665", "GO:0071339", "GO:0070603", "GO:1904949", "GO:0016514", "GO:0090544",
  "GO:0140513", "GO:1902494", "GO:0032991", "GO:0000785", "GO:0005667", "GO:1990234"
)
# terminos significativos (p ajustado < 0.1) en cualquier ranking, ademas de los pedidos
terminos <- union(terminos, gsea_df %>% filter(p.adjust < 0.1) %>% pull(ID))
gsea_comp <- gsea_df %>%
  filter(ID %in% terminos) %>%
  select(ID, Description, setSize, ranking, rep, NES, p.adjust) %>%
  mutate(rep = str_remove(rep, " ")) %>%
  pivot_wider(id_cols = c(ID, Description), names_from = c(ranking, rep), values_from = c(NES, p.adjust), names_glue = "{.value}_{ranking}_{rep}") %>%
  mutate(ID = factor(ID, levels = terminos)) %>%
  arrange(ID)
col_orden <- c("ID", "Description", unlist(map(c("AUC_original", names(modelos)), function(rk) map(c("Rep1", "Rep2"), ~ paste0(c("NES_", "p.adjust_"), rk, "_", .x)))))
gsea_comp <- gsea_comp %>% select(any_of(col_orden))
write_tsv(gsea_comp, file.path(out_dir, "p3_gsea_comparacion.tsv"))

# --- p4: miembros de MLL1/2 complex ---------------------------------------

mll_id <- "GO:0044665"
mll_genes <- suppressMessages(AnnotationDbi::select(hs, keys = mll_id, columns = "SYMBOL", keytype = "GOALL"))$SYMBOL %>% unique()
mll_tfs <- intersect(mll_genes, tfs)
# core enrichment (leading edge) del GSEA de MLL1/2 en cada ranking y replica
leading <- gsea_df %>%
  filter(ID == mll_id) %>%
  transmute(ranking, rep, sym = map(str_split(core_enrichment, "/"), ~ AnnotationDbi::mapIds(hs, .x, "SYMBOL", "ENTREZID"))) %>%
  unnest(sym)
mll_miembros <- coef_wide %>%
  filter(TF %in% mll_tfs) %>%
  mutate(
    leading_edge_AUC = map2_lgl(TF, rep, ~ .x %in% leading$sym[leading$ranking == "AUC_original" & leading$rep == .y]),
    leading_edge_ajustado = map2_lgl(TF, rep, ~ .x %in% leading$sym[leading$ranking == "ajustado" & leading$rep == .y])
  ) %>%
  arrange(rep, estimate_crudo)
write_tsv(mll_miembros, file.path(out_dir, "p4_mll1_2_miembros.tsv"))

# --- Figuras ----------------------------------------------------------------

# p1: dispersion crudo vs ajustado (los dos modelos ajustados), MLL1/2 destacados
disp <- coef_wide %>%
  pivot_longer(c(estimate_solo_actividad, estimate_ajustado), names_to = "modelo", values_to = "estimate_aj") %>%
  mutate(
    modelo = factor(modelos[str_remove(modelo, "estimate_")], levels = modelos[-1]),
    mll = TF %in% mll_tfs
  )
p_disp <- disp %>%
  ggplot(aes(estimate_crudo, estimate_aj)) +
  geom_hline(yintercept = 0, col = "grey60") +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_abline(linetype = "dashed", col = "grey40") +
  geom_point(data = ~ filter(.x, !mll), col = "grey70", alpha = 0.5, size = 1) +
  geom_point(data = ~ filter(.x, mll), col = "darkred", size = 2.2) +
  geom_text_repel(data = ~ filter(.x, mll), aes(label = TF), col = "darkred", size = 3, max.overlaps = 30, min.segment.length = 0) +
  facet_grid(modelo ~ rep, labeller = labeller(modelo = label_wrap_gen(30))) +
  coord_equal() +
  ggpubr::theme_pubr(base_size = 14) +
  theme(strip.background = element_blank()) +
  labs(
    x = "Efecto crudo de la unión\n(cambio en el rango de varianza, 1-100)", y = "Efecto ajustado de la unión",
    title = "Miembros de MLL1/2 complex (GO:0044665) en rojo"
  )
ggsave(file.path(out_dir, "p1_crudo_vs_ajustado.jpg"), p_disp, width = 11, height = 11, units = "in")

# p3: NES de los terminos pedidos en los cuatro rankings
p_nes <- gsea_df %>%
  filter(ID %in% terminos[1:12]) %>%
  mutate(
    ranking = factor(ranking, levels = c("AUC_original", names(modelos)), labels = c("AUC original", modelos)),
    Description = fct_reorder(paste0(Description, " (", ID, ")"), NES),
    sig = p.adjust < 0.1
  ) %>%
  ggplot(aes(NES, Description, col = ranking, shape = sig)) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_point(size = 3, position = position_dodge(0.6)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p ajustado < 0.1", `FALSE` = "p ajustado ≥ 0.1")) +
  scale_color_manual(values = c("black", "#1B8C8E", "#D6741F", "#AD343E")) +
  facet_wrap(~rep) +
  ggpubr::theme_pubr(base_size = 13) +
  theme(legend.position = "top", legend.direction = "vertical", strip.background = element_blank()) +
  labs(x = "NES (negativo = factores asociados a menor ruido)", y = NULL, col = "Ranking", shape = NULL)
ggsave(file.path(out_dir, "p3_gsea_nes_comparacion.jpg"), p_nes, width = 14, height = 8, units = "in")

# p3: curva de enriquecimiento de MLL1/2 sobre el ranking ajustado, por replica
for (r in repname) {
  p <- enrichplot::gseaplot2(gsea$ajustado[[r]], mll_id, title = paste0("MLL1/2 complex - ranking ajustado - ", r), subplots = 1:2, base_size = 14)
  ggsave(file.path(out_dir, paste0("p3_gsea_mll1_2_ajustado_", str_remove(r, " "), ".jpg")), p, width = 9, height = 6, units = "in")
}

# p4: efecto de cada miembro de MLL1/2 en los tres modelos
p_mll <- coef_df %>%
  filter(TF %in% mll_tfs) %>%
  mutate(
    modelo = factor(modelos[modelo], levels = modelos),
    TF = factor(TF, levels = mll_miembros %>% filter(rep == repname[1]) %>% arrange(estimate_crudo) %>% pull(TF)),
    sig = pval_corr < 0.05
  ) %>%
  ggplot(aes(estimate, TF, col = modelo, shape = sig)) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, position = position_dodge(0.7)) +
  geom_point(size = 2.5, position = position_dodge(0.7)) +
  geom_text(data = ~ distinct(.x, TF, rep, n_unidos), aes(x = Inf, y = TF, label = n_unidos), inherit.aes = FALSE, hjust = 1.1, size = 3.2, col = "grey30") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p corregido < 0.05", `FALSE` = "p corregido ≥ 0.05")) +
  scale_color_manual(values = c("black", "#D6741F", "#AD343E")) +
  facet_wrap(~rep) +
  ggpubr::theme_pubr(base_size = 13) +
  theme(legend.position = "top", legend.direction = "vertical", strip.background = element_blank()) +
  labs(x = "Efecto de la unión sobre el rango de varianza (IC 95%)", y = "Miembros de MLL1/2 complex en ReMap", col = NULL, shape = NULL,
    caption = "A la derecha: n de promotores unidos")
ggsave(file.path(out_dir, "p4_mll1_2_miembros.jpg"), p_mll, width = 13, height = 9, units = "in")

message("Exploratorio guardado en ", out_dir)
