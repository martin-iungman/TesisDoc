# discordancia_reportero_endogeno / 3_doble_glm_dispersion / remap_tf
#
# Analogo a prom_df_features.R (misma carpeta) pero con la union de cada
# factor de transcripcion (ChIP-seq ReMap2022) como feature en vez de
# las ~32 features curadas del promotor. Doble-GLM (glmmTMB con
# dispformula) que ajusta simultaneamente la MEDIA y la DISPERSION de
# dif_signed - ver README.md de la carpeta padre y el header de
# prom_df_features.R para la discusion completa de por que se prueba
# este modelo (la magnitud de discordancia resistio 5 intentos previos
# de desconfundirla de la actividad cruda).
#
# GSEA sobre la lista de TFs ordenada por el efecto de DISPERSION (no de
# media - ya sabemos que la media da r~0.27 con activity, razonable; acá
# se explora si el ranking de incertidumbre/dispersion tiene una
# estructura funcional GO propia, aunque el r~-0.8 con actividad cruda
# sugiere que va a estar dominado por el mismo eje).
#
# Requiere: data/processed/remap_tf_hits.tsv, data/processed/
# activity_stats_highconf.tsv, data/processed/prom_df.tsv,
# data/processed/fantom_endo_activity_summary.tsv (EXCEPCION, ver
# R/00_prom_features/analysis_tables_exceptions.R) y el output de R2.6
# (figures/R2.6_.../remap_tf_activity_wilcoxon.tsv). Run from the
# TesisDoc repo root. Tarda ~15-20 minutos: ~986 ajustes doble-GLM (uno
# por TF, ~0.7s c/u) + GSEA sobre anotaciones GO.

library(tidyverse)
library(glmmTMB)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(writexl)
source("R/functions/plot_helpers.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/discordancia_reportero_endogeno/3_doble_glm_dispersion"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hs <- org.Hs.eg.db

# --- Data prep -------------------------------------------------------------

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
  select(seq_id, rep, avg_rank, dif_signed) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, avg_rank, dif_signed), ~ replace_na(.x, FALSE)))

nTF <- data_TF %>%
  group_by(rep) %>%
  summarise(across(where(is.logical), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(-rep) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  filter(!str_detect(name, "^H\\d")) %>%
  distinct(name)

data_TF <- data_TF %>% mutate(rep = as.factor(rep))

# --- Doble-GLM: media y dispersion de dif_signed, por TF ------------------

fit_double_glm <- function(var) {
  df <- data_TF
  df$tf_test <- df[[var]]
  fit <- tryCatch(
    glmmTMB(dif_signed ~ avg_rank + rep + tf_test, dispformula = ~ avg_rank + tf_test, data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(variable = var, mean_estimate = NA_real_, mean_pval = NA_real_, disp_estimate = NA_real_, disp_pval = NA_real_))
  }
  s <- summary(fit)
  cc_mean <- s$coefficients$cond
  cc_disp <- s$coefficients$disp
  tibble(
    variable = var,
    mean_estimate = cc_mean["tf_testTRUE", "Estimate"], mean_pval = cc_mean["tf_testTRUE", "Pr(>|z|)"],
    disp_estimate = cc_disp["tf_testTRUE", "Estimate"], disp_pval = cc_disp["tf_testTRUE", "Pr(>|z|)"]
  )
}

message("Ajustando doble-GLM para ", length(nTF$name), " TFs (~15-20 min)...")
resultados <- map(nTF$name, fit_double_glm) %>%
  list_rbind() %>%
  mutate(mean_padj = p.adjust(mean_pval, method = "BH"), disp_padj = p.adjust(disp_pval, method = "BH"))
write_tsv(resultados, file.path(out_dir, "remap_tf_doble_glm.tsv"))

sig_colors <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")
add_sig <- function(df, padj_col) {
  df %>% mutate(sig = case_when(
    .data[[padj_col]] < 0.001 ~ "FDR < 0.001", .data[[padj_col]] < 0.01 ~ "FDR < 0.01",
    .data[[padj_col]] < 0.05 ~ "FDR < 0.05", TRUE ~ "ns"
  ) %>% factor(levels = names(sig_colors)))
}

# --- Plots: efecto sobre MEDIA y sobre DISPERSION (986 TFs, sin labels) ---

p_mean <- resultados %>% add_sig("mean_padj") %>% mutate(variable = fct_reorder(variable, mean_estimate)) %>%
  ggplot(aes(x = mean_estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = sig_colors) +
  labs(title = "Doble-GLM: efecto de cada TF (ReMap) sobre la MEDIA de la discordancia", x = "Estimado (β) — media", y = paste0("TFs (n=", nrow(resultados), ")"), color = NULL) +
  theme_bw(base_size = 14) + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(out_dir, "remap_tf_media.jpg"), p_mean, width = 9, height = 6.75, units = "in")

p_disp <- resultados %>% add_sig("disp_padj") %>% mutate(variable = fct_reorder(variable, disp_estimate)) %>%
  ggplot(aes(x = disp_estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = sig_colors) +
  labs(title = "Doble-GLM: efecto de cada TF (ReMap) sobre la DISPERSIÓN de la discordancia", x = "Estimado (β, escala log) — dispersión", y = paste0("TFs (n=", nrow(resultados), ")"), color = NULL) +
  theme_bw(base_size = 14) + theme(legend.position = "bottom", axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(out_dir, "remap_tf_dispersion.jpg"), p_disp, width = 9, height = 6.75, units = "in")

# --- Comparacion con el efecto sobre actividad cruda (R2.6), para media y
# para dispersion ----------------------------------------------------------

wilcox_activity <- read_tsv("figures/R2.6_remap_tf_activity_gsea/remap_tf_activity_wilcoxon.tsv", show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(activity_estimate = mean(estimate), .groups = "drop")

comparacion <- resultados %>% inner_join(wilcox_activity, by = c("variable" = "feature"))
r_mean <- cor(comparacion$mean_estimate, comparacion$activity_estimate, use = "complete.obs")
r_disp <- cor(comparacion$disp_estimate, comparacion$activity_estimate, use = "complete.obs")
write_tsv(comparacion, file.path(out_dir, "remap_tf_doble_glm_vs_actividad.tsv"))

comparacion_top <- comparacion %>%
  mutate(resid_disp = residuals(lm(disp_estimate ~ activity_estimate, data = comparacion))) %>%
  mutate(rank_resid = rank(-abs(resid_disp)), destacado = rank_resid <= 8)

p_comp_disp <- comparacion_top %>%
  ggplot(aes(activity_estimate, disp_estimate)) +
  geom_point(alpha = 0.25, col = "#358AAA") +
  geom_smooth(method = "lm", col = "#216869") +
  geom_text_repel(data = comparacion_top %>% filter(destacado), aes(label = variable), size = 3.5, fontface = "bold", color = "#AD343E", max.overlaps = Inf) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, label = paste0("r = ", round(r_disp, 2)), size = 6, label.size = 0) +
  labs(
    x = "Efecto sobre la actividad cruda (Wilcoxon, R2.6)", y = "Efecto sobre la dispersión (doble-GLM)",
    title = "TFs de ReMap: actividad cruda vs. dispersión reportero/endógeno",
    subtitle = "Etiquetadas las 8 mayores excepciones a la tendencia general"
  ) +
  theme_bw(base_size = 14)
ggsave(file.path(out_dir, "remap_tf_dispersion_vs_actividad.jpg"), p_comp_disp, width = 10, height = 7.5, units = "in")

message("r (media vs actividad): ", round(r_mean, 3), " | r (dispersion vs actividad): ", round(r_disp, 3))

# --- GSEA sobre la lista de TFs ordenada por el efecto de DISPERSION ------

ids <- AnnotationDbi::select(hs, keys = nTF$name, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct()
pregsea <- resultados %>%
  filter(!is.na(disp_estimate)) %>%
  arrange(desc(disp_estimate)) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))
ordered_disp <- pregsea$disp_estimate
names(ordered_disp) <- pregsea$ENTREZID

gse_go_disp <- gseGO(ordered_disp, ont = "all", OrgDb = "org.Hs.eg.db")
gse_go_disp_df <- as.data.frame(gse_go_disp) %>% arrange(desc(NES))
write_xlsx(list(GSEA_dispersion = gse_go_disp_df), file.path(out_dir, "GSEA_remap_dispersion.xlsx"))

gse_sig <- gse_go_disp_df %>% filter(p.adjust < 0.05)
if (nrow(gse_sig) > 0) {
  p_terms <- gse_sig %>%
    slice_max(order_by = abs(NES), n = 20) %>%
    mutate(
      dir = ifelse(NES > 0, "Mayor dispersión", "Menor dispersión"),
      Description = fct_reorder(Description, NES)
    ) %>%
    ggplot(aes(x = NES, y = Description, fill = dir)) +
    geom_col() +
    ggpubr::theme_pubr(base_size = 12) +
    theme(legend.position = "top") +
    labs(fill = "Efecto", x = "NES", y = NULL, title = "Términos GO enriquecidos (GSEA sobre efecto de dispersión de unión de TFs)") +
    scale_fill_manual(values = c("Mayor dispersión" = "#D6741F", "Menor dispersión" = "#7FB800"))
  ggsave(file.path(out_dir, "GSEA_terminos_enriquecidos_dispersion.jpg"), p_terms, width = 9, height = 8, units = "in")
  message(nrow(gse_sig), " terminos GO significativos (BH) en el GSEA de dispersion.")
} else {
  message("GSEA de dispersion no devolvio terminos significativos (p.adjust<0.05).")
}

message("Figuras guardadas en ", out_dir)
