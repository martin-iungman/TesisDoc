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
# ("## Remap noise analysis"), analogo a R/17_remap_activity_gsea (R2.6)
# pero con AUC/ROC sobre var_rank_sw en vez de Wilcoxon sobre mean.
# El original usa un objeto `tidy_TF_data` que nunca se define en el
# script - claramente pensado como el `data_TF` de la seccion de
# actividad (misma union de TFs) pero con var_rank_sw en vez de mean;
# reconstruido aca desde remap_tf_hits.tsv siguiendo el mismo patron
# que R2.6. IMPORTANTE: las columnas booleanas se convierten a
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
# histona (mismo criterio que R2.6).
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

# --- AUC (IC DeLong) de cada TF prediciendo ruido alto --------------------

auc_df <- map(repname, function(r) {
  df <- data_TF %>% filter(rep == r)
  map(vbles, function(f) {
    roc_obj <- pROC::roc(df[[f]], df[["var_rank_sw"]], ci = TRUE, direction = ">", quiet = TRUE)
    tibble(feature = f, rep = r, AUC = as.numeric(roc_obj$ci[2]), ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]))
  }) %>% list_rbind()
}) %>% list_rbind()

write_tsv(auc_df, file.path(out_dir, "remap_tf_noise_auc.tsv"))

panel <- plot_noise_auc_summary(auc_df, show_labels = FALSE)
ggsave(file.path(out_dir, "Remap_noise.jpg"), panel, width = 9, height = 6.75, units = "in")

# --- GSEA (GO Cellular Component) sobre la lista de TFs por AUC ----------

pregsea_noise <- function(rep_id) {
  ids <- AnnotationDbi::select(hs, keys = vbles, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% distinct()
  auc_df %>%
    filter(rep == rep_id) %>%
    arrange(desc(AUC)) %>%
    left_join(ids, by = c("feature" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID))
}

gse_go_noise <- map(repname, function(r) {
  pregsea <- pregsea_noise(r)
  ordered_noise <- pregsea$AUC
  names(ordered_noise) <- pregsea$ENTREZID
  gseGO(ordered_noise, ont = "CC", OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.1)
})
names(gse_go_noise) <- repname
gse_go_noise_df <- map(gse_go_noise, ~ as.data.frame(.x) %>% arrange(desc(NES)))

write_xlsx(gse_go_noise_df, file.path(out_dir, "GSEA_remap_noise.xlsx"))

# --- Panel: complejo MLL1 (GO:0071339) + TFs que lo integran --------------
# El original resalta especificamente este termino (encontrado
# significativo en su GSEA de ruido, Rep 1) - reproducido tal cual si
# aparece en nuestros resultados; si no, se señala en vez de forzarlo.

mll1_id <- "GO:0071339"
if (mll1_id %in% names(gse_go_noise[[repname[1]]]@geneSets) && mll1_id %in% gse_go_noise_df[[repname[1]]]$ID) {
  mll1_desc <- gse_go_noise_df[[repname[1]]]$Description[gse_go_noise_df[[repname[1]]]$ID == mll1_id]
  mll1_entrez <- gse_go_noise[[repname[1]]]@geneSets[[mll1_id]]
  pregsea_rep1 <- pregsea_noise(repname[1])
  mll1_tf <- pregsea_rep1$feature[pregsea_rep1$ENTREZID %in% mll1_entrez]

  p1 <- enrichplot::gseaplot2(gse_go_noise[[repname[1]]], mll1_id, title = mll1_desc, subplots = 1, color = "darkred")

  p2 <- auc_df %>%
    filter(rep == repname[1]) %>%
    mutate(noise = ifelse(AUC > 0.5, "Ruido alto", "Ruido bajo")) %>%
    filter((ci2.5 > 0.5 & ci97.5 > 0.5) | (ci2.5 < 0.5 & ci97.5 < 0.5)) %>%
    mutate(noise = ifelse(feature %in% mll1_tf, "set", noise)) %>%
    arrange(desc(AUC)) %>%
    ggplot(aes(x = AUC, y = fct_inorder(feature))) +
    geom_col(orientation = "y", position = "dodge", aes(fill = noise, alpha = noise)) +
    ggpubr::theme_pubr() +
    theme(text = element_text(size = 10), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(fill = "Efecto", x = "Efecto sobre el ruido", y = "Factores de transcripción") +
    scale_fill_manual(values = c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800", "set" = "darkred")) +
    scale_alpha_manual(values = c("Ruido alto" = 0.4, "Ruido bajo" = 0.4, "set" = 1)) +
    coord_flip()

  combined <- aplot::gglist(list(p1, p2), ncol = 1, heights = c(0.75, 0.5))
  ggsave(file.path(out_dir, "GSEA_MLL1.jpg"), combined, width = 9, height = 6.75, units = "in")
  message("GO:0071339 (MLL1 complex) encontrado y graficado en Rep 1.")
} else {
  message("GO:0071339 (MLL1 complex) NO aparece entre los resultados de GSEA de ruido en Rep 1 con pvalueCutoff=0.1 - no se genero GSEA_MLL1.jpg. Ver ", file.path(out_dir, "GSEA_remap_noise.xlsx"), " para los terminos que si salieron.")
}

message("Figuras guardadas en ", out_dir)
