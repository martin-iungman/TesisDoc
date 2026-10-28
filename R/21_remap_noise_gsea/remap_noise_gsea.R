# remap_tf_noise_gsea (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Efecto de la union de cada factor de transcripcion (ChIP-seq
# ReMap2022, cualquier linea celular) sobre el ruido transcripcional
# (var_rank_sw - ver R/19_ruido_cgi_tata), y GSEA (GO Cellular
# Component) sobre esa lista ordenada por AUC.
#
# Requiere: data/processed/remap_tf_hits.tsv (ver
# R/00_prom_features/build_remap_tf_hits.R), data/processed/
# activity_stats_highconf.tsv y data/processed/prom_df.tsv. Run from
# the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("## Remap noise analysis"), analogo a R/17_remap_activity_gsea (remap_tf_activity_gsea)
# pero con AUC/ROC sobre var_rank_sw en vez de Wilcoxon sobre mean.
# El original usa un objeto `tidy_TF_data` que nunca se define en el
# script - claramente pensado como el `data_TF` de la seccion de
# actividad (misma union de TFs) pero con var_rank_sw en vez de mean;
# reconstruido aca desde remap_tf_hits.tsv siguiendo el mismo patron
# que remap_tf_activity_gsea. IMPORTANTE: las columnas booleanas se convierten a
# factor(levels=c("TRUE","FALSE")) antes de pROC::roc(..., direction=
# ">") - igual que build_tidy_features() y Histone_chipatlas.qmd's
# tidy_hist_data (que sí definia esto para el analisis de ruido, a
# diferencia de la seccion de actividad). Sin esto pROC ordena los
# niveles de un logical alfabeticamente (FALSE=control, TRUE=caso) y
# direction=">" queda invertido - confirmado con datos reales de TATA-
# box: AUC=0.551 (correcto, TATA-box tiene mas ruido) vs AUC=0.449 (el
# invertido) para la misma comparacion.

library(tidyverse)
library(pROC)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(writexl)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

# clusterProfiler/AnnotationDbi definen genericos S4 que tapan
# dplyr::select/filter/rename; fijarlos a las versiones de dplyr.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

BiocParallel::register(BiocParallel::SerialParam())

slug <- "remap_tf_noise_gsea"
out_dir <- fig_dir(slug)
hs <- org.Hs.eg.db

# --- Data prep: ruido x union a cada TF (cualquier linea celular) --------

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data <- add_noise_rank(data)

binary_df <- remap_hits %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, var_rank_sw) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, var_rank_sw), ~ replace_na(.x, FALSE)))

# TFs con >100 promotores unidos en AMBAS replicas, excluyendo marcas de
# histona (mismo criterio que remap_tf_activity_gsea).
nTF <- data_TF %>%
  group_by(rep) %>%
  summarise(across(where(is.logical), sum, na.rm = TRUE)) %>%
  pivot_longer(-rep) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  filter(!str_detect(name, "^H\\d")) %>%
  distinct(name)

vbles <- nTF$name
repname <- unique(data_TF$rep)

data_TF <- data_TF %>% mutate(across(all_of(vbles), ~ factor(.x, levels = c("TRUE", "FALSE"))))

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}

# --- AUC de cada TF prediciendo ruido alto --------------------------------
# Significancia igual que ruido_summary: p de Wald con la varianza de DeLong
# de la AUC, corregido por BH sobre TFs x replicas.

auc_df <- map(repname, function(r) {
  df <- data_TF %>% filter(rep == r)
  map(vbles, function(f) {
    roc_obj <- pROC::roc(df[[f]], df[["var_rank_sw"]], ci = TRUE, direction = ">", quiet = TRUE)
    auc_val <- as.numeric(roc_obj$auc)
    se_delong <- sqrt(pROC::var(roc_obj, method = "delong"))
    tibble(
      feature = f, rep = r, AUC = auc_val, ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]),
      pval = 2 * pnorm(-abs((auc_val - 0.5) / se_delong))
    )
  }) %>% list_rbind()
}) %>% list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, method = "BH"))

write_tsv(auc_df, file.path(out_dir, "remap_tf_noise_auc.tsv"))

# TFs con pval_corr < 0.05 y el mismo sentido en ambas replicas (los que
# grafica plot_noise_auc_summary)
n_consistentes <- auc_df %>%
  filter(pval_corr < 0.05) %>%
  group_by(feature) %>%
  filter(n() == 2, n_distinct(AUC > 0.5) == 1) %>%
  ungroup() %>%
  distinct(feature, ruido_alto = AUC > 0.5) %>%
  count(ruido_alto)
print(n_consistentes)

panel_auc <- plot_noise_auc_summary(auc_df, show_labels = FALSE, show_errorbars = FALSE) +
  labs(y = "Factores de transcripción (ReMap)")
save_jpg_pdf(panel_auc, "Remap_noise", 9, 6.75)

# --- GSEA (GO Cellular Component) sobre la lista de TFs por AUC ----------

ids <- AnnotationDbi::select(hs, keys = vbles, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% distinct()
pregsea_noise <- function(rep_id) {
  auc_df %>%
    filter(rep == rep_id) %>%
    arrange(desc(AUC)) %>%
    left_join(ids, by = c("feature" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
}

set.seed(2024)
gse_go_noise <- map(repname, function(r) {
  pregsea <- pregsea_noise(r)
  gseGO(setNames(pregsea$AUC, pregsea$ENTREZID), ont = "CC", OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.1, seed = TRUE)
})
names(gse_go_noise) <- repname
gse_go_noise_df <- map(gse_go_noise, ~ as.data.frame(.x) %>% arrange(p.adjust))

write_xlsx(gse_go_noise_df, file.path(out_dir, "GSEA_remap_noise.xlsx"))

# --- Panel: complejo MLL1/2 (GO:0044665) + TFs que lo integran ------------
# El original resaltaba MLL1 complex (GO:0071339), significativo solo en
# Rep 1; MLL1/2 complex (lo contiene) es significativo en ambas replicas.
# NES negativo: los TFs del complejo se concentran en el extremo de AUC mas
# baja (mayor reduccion del ruido). Se grafica Rep 1.

mll_id <- "GO:0044665"
stopifnot(all(map_lgl(gse_go_noise_df, ~ mll_id %in% .x$ID)))
mll_desc <- gse_go_noise_df[[repname[1]]]$Description[gse_go_noise_df[[repname[1]]]$ID == mll_id]
pregsea_rep1 <- pregsea_noise(repname[1])
mll_tf <- pregsea_rep1$feature[pregsea_rep1$ENTREZID %in% gse_go_noise[[repname[1]]]@geneSets[[mll_id]]]
message(mll_desc, " (Rep 1): ", paste(mll_tf, collapse = ", "))

p1 <- enrichplot::gseaplot2(gse_go_noise[[repname[1]]], mll_id, title = mll_desc, subplots = 1, color = "darkred", base_size = 20) +
  labs(y = "Puntaje de\nenriquecimiento")

# todos los TFs de la lista del GSEA de Rep 1, en el mismo orden
p2 <- pregsea_rep1 %>%
  mutate(noise = case_when(feature %in% mll_tf ~ "set", AUC > 0.5 ~ "Ruido alto", TRUE ~ "Ruido bajo")) %>%
  ggplot(aes(x = AUC - 0.5, y = fct_inorder(feature))) +
  geom_col(orientation = "y", aes(fill = noise, alpha = noise)) +
  scale_x_continuous(labels = function(x) x + 0.5) +
  ggpubr::theme_pubr() +
  theme(text = element_text(size = 20), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "AUC (efecto\nsobre el ruido)", y = "Factores de transcripción") +
  scale_fill_manual(values = c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800", "set" = "darkred")) +
  scale_alpha_manual(values = c("Ruido alto" = 0.4, "Ruido bajo" = 0.4, "set" = 1)) +
  coord_flip()

combined <- aplot::gglist(list(p1, p2), ncol = 1, heights = c(0.75, 0.5))
save_jpg_pdf(combined, "GSEA_MLL1_2", 12, 8)

# Figura completa: AUC de todos los TFs consistentes (A) arriba y el GSEA de
# MLL1/2 complex en Rep 1 (B) abajo, mismo formato que el panel C de
# remap_tf_activity_gsea
h_a <- 6.75
h_b <- 8
figura <- cowplot::plot_grid(panel_auc, ggplotify::as.ggplot(combined), ncol = 1, rel_heights = c(h_a, h_b)) +
  cowplot::draw_plot_label(c("A", "B"), x = 0, y = c(1, h_b / (h_a + h_b)), size = 26)
save_jpg_pdf(figura, "remap_tf_noise_gsea", 12, h_a + h_b)

message("Figuras guardadas en ", out_dir)
