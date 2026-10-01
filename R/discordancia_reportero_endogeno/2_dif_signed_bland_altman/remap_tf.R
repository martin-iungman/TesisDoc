# discordancia_reportero_endogeno / 2_dif_signed_bland_altman / remap_tf
#
# SACADO DEL PIPELINE NUMERADO 2026-07-30 (era R4.7/R/28_remap_rank_dif_
# lmm_gsea) y puesto en esta carpeta de staging junto con las otras 2
# variantes del analisis (1_original_sin_corregir, 3_doble_glm_
# dispersion) mientras el autor decide cual usar para la tesis - no
# escribe via fig_dir()/mapping_figuras.csv, ver README.md de la carpeta
# padre.
#
# Analogo a 2_dif_signed_bland_altman/prom_df_features.R pero con la union de cada
# transcripcion (ChIP-seq ReMap2022, cualquier linea celular) como
# feature en vez de las ~32 features curadas del promotor: para cada TF,
# dif_signed ~ avg_rank + TF_bound + (1|replica), con correccion BH sobre
# el efecto fijo del TF. GSEA (terminos GO) sobre esa lista de TFs
# ordenada por efecto - mismo patron que R2.6 (remap_tf_activity_gsea) y
# R3.3 (remap_tf_noise_gsea) pero con el estimado del LMM en vez de
# Wilcoxon/AUC como metrica de ranking.
#
# dif_signed = rr - re (rangos de reportero/endogena normalizados 0-1
# por replica; convencion elegida por el autor 2026-07-31): positivo =
# el reportero SOBREESTIMA respecto a la actividad endogena, negativo =
# el reportero SUBESTIMA.
#
# SIGNO DEL COEFICIENTE POR TF (`estimate` en resultados_lmm/los tsv):
# cada TF es un factor(levels=c("TRUE","FALSE")), TRUE=referencia, asi
# que lmer/broom devuelven crudo el contraste FALSE-menos-TRUE
# (contraintuitivo para un grafico rotulado por nombre de TF). Por eso,
# justo despues de armar resultados_lmm, invertimos el signo (y el IC)
# para que quede como el contraste directo TRUE-menos-FALSE: de ahi en
# mas (incluyendo el resto de este script - panel_neg, GSEA, ORA, la
# comparacion con actividad cruda), `estimate` positivo = el TF unido
# predice dif_signed mas alto = el reportero SOBREESTIMA para ese TF;
# negativo = SUBESTIMA (direccion dominante, 878/986 TFs). Verificado 4
# formas distintas (medias crudas por grupo, modelo con contraste por
# defecto FALSE-como-referencia, el modelo tal cual esta en el pipeline,
# y releveled a mano) - ver sesion 2026-07-31.
#
# avg_rank (promedio de ambos rangos) es el control de nivel de
# actividad - ver R/27_rank_dif_lmm/rank_dif_lmm.R para la discusion
# completa de por que reemplaza a controlar por rank_reporter solo
# (REDISEÑADO 2026-07-29: la version anterior de este script, con
# dif_rank_rank ~ rank_reporter_scaled + TF, mostraba una correlacion de
# r~0.8 entre el efecto de cada TF sobre la "diferencia de rango" y su
# efecto sobre la actividad cruda - ver remap_rank_dif_vs_activity.jpg
# mas abajo para el valor actualizado con el diseño corregido, deberia
# rondar r~0.27).
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
#
# No portado de un script del repo viejo (no existe un equivalente de
# rank_dif_lmm.R para TFs de ReMap alli) - nuevo analisis, combinando
# R/27_rank_dif_lmm/rank_dif_lmm.R (el modelo LMM) con la estructura de
# datos/GSEA de R/17_remap_activity_gsea.R y R/21_remap_noise_gsea.R.

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

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hs <- org.Hs.eg.db

# --- Data prep: discordancia con signo x union a cada TF ------------------

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance()

binary_df <- remap_hits %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, dif_signed, avg_rank) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, dif_signed, avg_rank), ~ replace_na(.x, FALSE)))

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

# --- LMM: efecto de la union a cada TF sobre la discordancia con signo ---

resultados_lmm <- map(nTF$name, function(var) {
  formula <- as.formula(paste(
    "dif_signed ~ avg_rank +",
    paste0("`", var, "`"), "+ (1 | rep)"
  ))
  fit <- lmerTest::lmer(formula, data = data_TF, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = var)
}) %>%
  list_rbind() %>%
  mutate(
    conf.low_true = -conf.high,
    conf.high_true = -conf.low,
    estimate = -estimate,
    conf.low = conf.low_true,
    conf.high = conf.high_true
  ) %>%
  select(-conf.low_true, -conf.high_true) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(resultados_lmm, file.path(out_dir, "remap_tf_rank_dif_lmm.tsv"))

# --- Plots: efecto de union sobre la discordancia con signo ---------------
# Misma estructura (dot+errorbar, color por categoria de FDR) que R4.6
# (rank_dif_lmm.R) y el panel de histonas (chipatlas_histonas_rank_dif_
# lmm.R) - los 3 paneles de "discordancia" (prom_df, TFs, histonas) se
# leen como una sola familia de figuras.

sig_colors <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")
add_sig <- function(df) {
  df %>% mutate(sig = case_when(
    p_adj < 0.001 ~ "FDR < 0.001",
    p_adj < 0.01 ~ "FDR < 0.01",
    p_adj < 0.05 ~ "FDR < 0.05",
    TRUE ~ "ns"
  ) %>% factor(levels = names(sig_colors)))
}

# Panel con los 986 TFs testeados (sin labels - a diferencia de R4.6/
# histonas, acá no entran ~1000 nombres en el eje y).
panel_all <- resultados_lmm %>%
  add_sig() %>%
  mutate(variable = fct_reorder(variable, estimate)) %>%
  ggplot(aes(x = estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Efecto de la unión de cada TF (ReMap) sobre la discordancia (con signo)\nentre actividad del reportero y endógena",
    subtitle = "Discordancia (reportero - endo) ~ TF + nivel de actividad + (1|réplica)",
    x = "Estimado (β) — positivo: el reportero sobreestima", y = "Factores de transcripción (ReMap, n=986)", color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(out_dir, "remap_rank_dif_lmm.jpg"), panel_all, width = 9, height = 6.75, units = "in")

# Panel con los TOP 20 TFs de efecto negativo (subestima, por
# |estimate|), no todos los significativos (878 - con el diseño
# corregido, "subestima" es la direccion dominante, no un subconjunto
# chico, asi que graficar todos con labels seria ilegible). `estimate`
# ya esta invertido a contraste TRUE-menos-FALSE (ver nota de signo mas
# arriba): negativo = el TF unido predice dif_signed mas bajo = SUBESTIMA.
panel_neg <- resultados_lmm %>%
  filter(p_adj < 0.05, estimate < 0) %>%
  add_sig() %>%
  slice_max(order_by = abs(estimate), n = 20) %>%
  mutate(variable = fct_reorder(variable, estimate)) %>%
  ggplot(aes(x = estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, alpha = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Top 20 TFs donde el reportero subestima más la actividad respecto a la endógena",
    x = "Estimado (β) con IC 95%", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "remap_rank_dif_lmm_top_subestima.jpg"), panel_neg, width = 9, height = 7, units = "in")

# --- Comparacion con el efecto sobre actividad cruda (R2.6) ---------------
# Con el diseño anterior (control por rank_reporter solo) esto daba
# r~0.8 - ver R/27_rank_dif_lmm/rank_dif_lmm.R para la discusion completa
# de por que, y por que este diseño (avg_rank + dif_signed, tipo
# Bland-Altman) lo reduce a ~0.27.

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
    y = "Efecto sobre la discordancia con signo (LMM, controlando avg_rank)",
    title = "Efecto de cada TF: actividad cruda vs. discordancia reportero/endógena"
  ) +
  theme_bw(base_size = 14)
ggsave(file.path(out_dir, "remap_rank_dif_vs_activity.jpg"), p_comparacion, width = 9, height = 6.75, units = "in")
message("Correlacion entre efecto de actividad (R2.6) y efecto de discordancia con signo: r=", round(r_val, 3))

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

# --- Panel: terminos GO enriquecidos (top por NES, p.adjust<0.05) --------

gse_go_sig <- gse_go_df %>% filter(p.adjust < 0.05)
if (nrow(gse_go_sig) > 0) {
  p_terms <- gse_go_sig %>%
    slice_max(order_by = abs(NES), n = 20) %>%
    mutate(
      # NES>0 = extremo de mayor `estimate` (ya invertido a TRUE-menos-
      # FALSE, ver nota de signo mas arriba) = grupo TRUE con dif_signed
      # MAS ALTO = reportero SOBREESTIMA.
      dif = ifelse(NES > 0, "Reportero sobreestima", "Reportero subestima"),
      Description = fct_reorder(Description, NES)
    ) %>%
    ggplot(aes(x = NES, y = Description, fill = dif)) +
    geom_col() +
    ggpubr::theme_pubr(base_size = 12) +
    theme(legend.position = "top") +
    labs(fill = "Efecto", x = "NES", y = NULL, title = "Términos GO enriquecidos (GSEA sobre efecto de union de TFs)") +
    scale_fill_manual(values = c("Reportero sobreestima" = "#D6741F", "Reportero subestima" = "#7FB800"))
  ggsave(file.path(out_dir, "GSEA_terminos_enriquecidos_remap_rank_dif.jpg"), p_terms, width = 9, height = 8, units = "in")
} else {
  message("GSEA no devolvio terminos significativos (p.adjust<0.05) - no se genero el panel de terminos enriquecidos. Ver ", file.path(out_dir, "GSEA_remap_rank_dif.xlsx"), ".")
}

# --- ORA: enriquecimiento GO en los TFs con efecto POSITIVO (el
# reportero SOBREESTIMA respecto a la actividad endogena, p_adj<0.05)
# frente al resto de los TFs testeados (universo = todos los TFs con
# mapeo a ENTREZID) - a diferencia del GSEA de arriba (ranking continuo
# por estimado), esto es una prueba de sobre-representacion (hipergeo-
# metrica) sobre un set discreto. Corre sobre el set chico (sobreestima,
# `estimate`>0 ya invertido a TRUE-menos-FALSE - ver nota de signo mas
# arriba), no el negativo/subestima, que con el diseño corregido es la
# enorme mayoria - 878/986 TFs - y por eso deja de ser un subconjunto
# informativo para ORA: casi todo el universo cae adentro del "gene
# set", sin contraste real. El set sobreestima es chico (n=6, ver mas
# abajo) - el ORA tiene poca potencia con tan pocos genes, pero es la
# comparacion correcta.
# ---------------------------------------------------------------------

pos_tf_entrez <- resultados_lmm %>%
  filter(p_adj < 0.05, estimate > 0) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID)

ora_pos <- enrichGO(
  gene = pos_tf_entrez,
  universe = ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "all",
  pAdjustMethod = "BH"
)
ora_pos_df <- as.data.frame(ora_pos) %>% arrange(p.adjust)
write_xlsx(list(ORA_efecto_positivo = ora_pos_df), file.path(out_dir, "ORA_remap_rank_dif_efecto_positivo.xlsx"))

if (nrow(ora_pos_df) > 0) {
  p_ora <- ora_pos_df %>%
    slice_min(order_by = p.adjust, n = 20) %>%
    mutate(Description = fct_reorder(Description, -p.adjust)) %>%
    ggplot(aes(x = -log10(p.adjust), y = Description, fill = Count)) +
    geom_col() +
    ggpubr::theme_pubr(base_size = 12) +
    labs(
      x = "-log10(p ajustado)", y = NULL, fill = "N° de TFs",
      title = paste0("Enriquecimiento GO (ORA) en TFs donde el reportero sobreestima\nrespecto a la actividad endógena vs. el resto de los TFs testeados (n=", length(pos_tf_entrez), ")")
    ) +
    scale_fill_gradient(low = "#7FB800", high = "#0D2C54")
  ggsave(file.path(out_dir, "ORA_remap_rank_dif_efecto_positivo.jpg"), p_ora, width = 9, height = 8, units = "in")
  message(nrow(ora_pos_df), " terminos GO significativos (BH) en el ORA de TFs donde el reportero sobreestima (n=", length(pos_tf_entrez), " TFs).")
} else {
  message("ORA no devolvio terminos significativos para el set de TFs con efecto positivo (n=", length(pos_tf_entrez), " TFs) - ver ", file.path(out_dir, "ORA_remap_rank_dif_efecto_positivo.xlsx"), " (vacio).")
}

message("Figuras guardadas en ", out_dir)
