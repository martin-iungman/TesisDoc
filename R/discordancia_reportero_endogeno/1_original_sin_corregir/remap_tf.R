# discordancia_reportero_endogeno / 1_original_sin_corregir / remap_tf
#
# VERSION ORIGINAL (sin corregir) del analisis de "diferencia de rango"
# aplicado a TFs de ReMap. Se mantiene aca (2026-07-30) como referencia
# de comparacion: el panel remap_rank_dif_vs_activity.jpg de este mismo
# script muestra el diagnostico que motivo la correccion - el efecto de
# cada TF sobre esta "diferencia de rango" correlaciona r~0.8 con su
# efecto sobre la actividad cruda (R2.6), pese a controlar por
# rank_reporter_scaled - es decir, el analisis practicamente no aporta
# nada mas alla de re-detectar que TFs predicen actividad. Ver
# 2_dif_signed_bland_altman/remap_tf.R para el diseño corregido (r~0.27)
# y 3_doble_glm_dispersion/remap_tf.R para el intento de aislar la
# magnitud/incertidumbre (que tampoco se pudo desconfundir). No escribe
# via fig_dir()/mapping_figuras.csv, ver README.md de la carpeta padre.
#
# Analogo al modelo LMM de prom_df_features.R (misma carpeta) pero con la
# union de cada factor de transcripcion (ChIP-seq ReMap2022, cualquier
# linea celular) como feature en vez de las ~32 features curadas del
# promotor: para cada TF, dif_rank_rank ~ rank_reporter_scaled + TF_bound
# + (1|replica), con correccion BH sobre el efecto fijo del TF. GSEA
# (terminos GO) sobre esa lista de TFs ordenada por efecto, y ORA
# (sobre-representacion) en el set de TFs con efecto negativo.
#
# Requiere: data/processed/remap_tf_hits.tsv (ver
# R/00_prom_features/build_remap_tf_hits.R), data/processed/
# activity_stats_highconf.tsv, data/processed/prom_df.tsv,
# data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R) y el output de
# R/17_remap_activity_gsea.R (figures/R2.6_.../remap_tf_activity_
# wilcoxon.tsv, para el panel de comparacion). Run from the TesisDoc
# repo root, despues de haber corrido R2.6 al menos una vez. Tarda varios
# minutos: ~986 ajustes LMM (uno por TF, ~0.2s c/u) + GSEA sobre
# anotaciones GO.

library(tidyverse)
library(lmerTest)
library(broom.mixed)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(writexl)
source("R/functions/plot_helpers.R")

# clusterProfiler/AnnotationDbi definen genericos S4 que tapan
# dplyr::select/filter/rename; fijarlos a las versiones de dplyr (mismo
# guard que R2.6/R3.3).
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/discordancia_reportero_endogeno/1_original_sin_corregir"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hs <- org.Hs.eg.db

# --- Data prep: diferencia de rango x union a cada TF ----------------------

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  mutate(max_tpm = replace_na(max_tpm, 0)) %>%
  group_by(rep) %>%
  mutate(
    rank_reporter = dense_rank(mean),
    rank_endo = row_number(max_tpm),
    dif_rank_rank = rank(abs(rank_endo - rank_reporter)),
    rank_reporter_scaled = scale(rank_reporter)[, 1]
  ) %>%
  ungroup()

binary_df <- remap_hits %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, dif_rank_rank, rank_reporter_scaled) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, dif_rank_rank, rank_reporter_scaled), ~ replace_na(.x, FALSE)))

# TFs con >100 promotores unidos en AMBAS replicas, excluyendo marcas de
# histona (mismo criterio que R2.6/R3.3).
nTF <- data_TF %>%
  group_by(rep) %>%
  summarise(across(where(is.logical), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(-rep) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  filter(!str_detect(name, "^H\\d")) %>%
  distinct(name)

data_TF <- data_TF %>%
  mutate(
    rep = as.factor(rep),
    across(all_of(nTF$name), ~ factor(.x, levels = c("TRUE", "FALSE")))
  )

# --- LMM: efecto de la union a cada TF sobre la diferencia de rango -------

resultados_lmm <- map(nTF$name, function(var) {
  formula <- as.formula(paste(
    "dif_rank_rank ~ rank_reporter_scaled +",
    paste0("`", var, "`"), "+ (1 | rep)"
  ))
  fit <- lmerTest::lmer(formula, data = data_TF, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = var)
}) %>%
  list_rbind() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(resultados_lmm, file.path(out_dir, "remap_tf_rank_dif_lmm.tsv"))

# --- Plots: efecto de union sobre la diferencia de rango -------------------

sig_colors <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")
add_sig <- function(df) {
  df %>% mutate(sig = case_when(
    p_adj < 0.001 ~ "FDR < 0.001",
    p_adj < 0.01 ~ "FDR < 0.01",
    p_adj < 0.05 ~ "FDR < 0.05",
    TRUE ~ "ns"
  ) %>% factor(levels = names(sig_colors)))
}

panel_all <- resultados_lmm %>%
  add_sig() %>%
  mutate(variable = fct_reorder(variable, estimate)) %>%
  ggplot(aes(x = estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "[ORIGINAL, SIN CORREGIR] Efecto de la unión de cada TF (ReMap) sobre la diferencia\nde rango entre actividad endógena y del reportero",
    subtitle = "Diferencia de rango ~ TF + actividad del reportero + (1|réplica)",
    x = "Estimado (β)", y = "Factores de transcripción (ReMap, n=986)", color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(out_dir, "remap_rank_dif_lmm.jpg"), panel_all, width = 9, height = 6.75, units = "in")

panel_menor <- resultados_lmm %>%
  filter(p_adj < 0.05, estimate < 0) %>%
  add_sig() %>%
  mutate(variable = fct_reorder(variable, estimate)) %>%
  ggplot(aes(x = estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, alpha = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "TFs con menor diferencia de rango (mejor concordancia reportero/endógena)",
    x = "Estimado (β) con IC 95%", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "remap_rank_dif_lmm_menor.jpg"), panel_menor, width = 9, height = 6, units = "in")

# --- Comparacion con el efecto sobre actividad cruda (R2.6) ---------------
# ESTE ES EL DIAGNOSTICO CLAVE: si el analisis de "diferencia de rango"
# aportara algo mas alla de actividad cruda, no deberia correlacionar
# fuerte con el efecto de Wilcoxon de R2.6. Da r~0.8 - el analisis esta
# dominado por el mismo efecto que ya mide R2.6.

wilcox_activity <- read_tsv("figures/R2.6_remap_tf_activity_gsea/remap_tf_activity_wilcoxon.tsv", show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(activity_estimate = mean(estimate), .groups = "drop")

comparacion <- resultados_lmm %>% inner_join(wilcox_activity, by = c("variable" = "feature"))
r_val <- cor(comparacion$activity_estimate, comparacion$estimate)
write_tsv(comparacion, file.path(out_dir, "remap_rank_dif_vs_activity.tsv"))

p_comparacion <- comparacion %>%
  ggplot(aes(activity_estimate, estimate)) +
  geom_point(alpha = 0.3, col = "#358AAA") +
  geom_smooth(method = "lm", col = "#216869") +
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, label = paste0("r = ", round(r_val, 2)), size = 6, label.size = 0) +
  labs(
    x = "Efecto sobre la actividad cruda (Wilcoxon, R2.6)",
    y = "Efecto sobre la diferencia de rango (LMM, sin corregir)",
    title = "[ORIGINAL, SIN CORREGIR] ¿Aporta algo mas alla de actividad cruda? (spoiler: no - r alto)"
  ) +
  theme_bw(base_size = 14)
ggsave(file.path(out_dir, "remap_rank_dif_vs_activity.jpg"), p_comparacion, width = 9, height = 6.75, units = "in")
message("Correlacion entre efecto de actividad (R2.6) y efecto de diferencia de rango (SIN CORREGIR): r=", round(r_val, 3))

# --- GSEA sobre la lista de TFs ordenada por el estimado del LMM ----------

ids <- AnnotationDbi::select(hs, keys = nTF$name, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct()
pregsea <- resultados_lmm %>%
  arrange(desc(estimate)) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))
ordered_genes <- pregsea$estimate
names(ordered_genes) <- pregsea$ENTREZID

gse_go <- gseGO(ordered_genes, ont = "all", OrgDb = "org.Hs.eg.db")
gse_go_df <- as.data.frame(gse_go) %>% arrange(desc(NES))
write_xlsx(list(GSEA = gse_go_df), file.path(out_dir, "GSEA_remap_rank_dif.xlsx"))

gse_go_sig <- gse_go_df %>% filter(p.adjust < 0.05)
if (nrow(gse_go_sig) > 0) {
  p_terms <- gse_go_sig %>%
    slice_max(order_by = abs(NES), n = 20) %>%
    mutate(
      dif = ifelse(NES > 0, "Mayor diferencia de rango", "Menor diferencia de rango"),
      Description = fct_reorder(Description, NES)
    ) %>%
    ggplot(aes(x = NES, y = Description, fill = dif)) +
    geom_col() +
    ggpubr::theme_pubr(base_size = 12) +
    theme(legend.position = "top") +
    labs(fill = "Efecto", x = "NES", y = NULL, title = "[ORIGINAL] Términos GO enriquecidos (GSEA sobre efecto de union de TFs)") +
    scale_fill_manual(values = c("Mayor diferencia de rango" = "#D6741F", "Menor diferencia de rango" = "#7FB800"))
  ggsave(file.path(out_dir, "GSEA_terminos_enriquecidos_remap_rank_dif.jpg"), p_terms, width = 9, height = 8, units = "in")
} else {
  message("GSEA no devolvio terminos significativos (p.adjust<0.05) - no se genero el panel de terminos enriquecidos.")
}

# --- ORA: enriquecimiento GO en los TFs con efecto negativo ---------------

neg_tf_entrez <- resultados_lmm %>%
  filter(p_adj < 0.05, estimate < 0) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID)

ora_neg <- enrichGO(
  gene = neg_tf_entrez,
  universe = ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "all",
  pAdjustMethod = "BH"
)
ora_neg_df <- as.data.frame(ora_neg) %>% arrange(p.adjust)
write_xlsx(list(ORA_efecto_negativo = ora_neg_df), file.path(out_dir, "ORA_remap_rank_dif_efecto_negativo.xlsx"))

if (nrow(ora_neg_df) > 0) {
  p_ora <- ora_neg_df %>%
    slice_min(order_by = p.adjust, n = 20) %>%
    mutate(Description = fct_reorder(Description, -p.adjust)) %>%
    ggplot(aes(x = -log10(p.adjust), y = Description, fill = Count)) +
    geom_col() +
    ggpubr::theme_pubr(base_size = 12) +
    labs(
      x = "-log10(p ajustado)", y = NULL, fill = "N° de TFs",
      title = paste0("[ORIGINAL] Enriquecimiento GO (ORA) en TFs con efecto negativo (n=", length(neg_tf_entrez), ")")
    ) +
    scale_fill_gradient(low = "#7FB800", high = "#0D2C54")
  ggsave(file.path(out_dir, "ORA_remap_rank_dif_efecto_negativo.jpg"), p_ora, width = 9, height = 8, units = "in")
  message(nrow(ora_neg_df), " terminos GO significativos (BH) en el ORA de TFs con efecto negativo (n=", length(neg_tf_entrez), " TFs).")
} else {
  message("ORA no devolvio terminos significativos para el set de TFs con efecto negativo (n=", length(neg_tf_entrez), " TFs).")
}

message("Figuras guardadas en ", out_dir)
