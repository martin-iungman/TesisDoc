# remap_tf_activity_gsea (ver docs/mapping_figuras.csv para el numero de
# figura vigente - todavia sin asignar)
# Efecto de la union de cada factor de transcripcion (ChIP-seq ReMap2022,
# union en cualquier linea celular) sobre la actividad transcripcional
# media, y enriquecimiento funcional (GSEA sobre terminos GO) de esa
# lista de TFs ordenada por efecto.
#
# Requiere: data/processed/remap_tf_hits.tsv (ver
# R/00_prom_features/build_remap_tf_hits.R), data/processed/
# activity_stats_highconf.tsv y data/processed/prom_df.tsv (ver
# R/00_prom_features y R/01_activity_stats). Run from the TesisDoc repo
# root. Tarda varios minutos: ~500-800 tests de Wilcoxon (uno por TF) +
# GSEA sobre anotaciones GO.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("Whole Remap - activity analysis" y "GSEA Remap activity"). Se
# corrigieron dos bugs del script original en pregsea_act(): (1)
# `filter(rep==rep)` comparaba la columna contra si misma (siempre
# TRUE, no filtraba nada); (2) `pregsea[[rep]]` indexaba con un string
# ("Rep 1") una lista sin nombres devuelta por group_split() (esa
# lista se indexa por posicion, no por nombre) - se reescribio para
# filtrar directamente por rep en vez de group_split()+indexar. Tambien
# se unifico el top pathway usado para el panel POL2-like: el original
# tomaba el ID del objeto gseGO sin ordenar (GSE_GO_rep1$ID[1]) para
# extraer el gene set pero el ID ordenado por NES (GSE_GO_rep1_df$ID[1])
# para el titulo del plot - podian no coincidir. Aca se usa el ordenado
# por NES para ambos.

library(tidyverse)
library(coin)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(writexl)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

# clusterProfiler/AnnotationDbi define S4 generics that shadow dplyr's
# select/filter/rename; pin them to the dplyr versions.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# fgsea (used internally by gseGO) defaults to BiocParallel's multicore
# backend, which spawns separate worker processes and on Windows/this
# machine ran out of memory on the 2nd gseGO() call ("worker evaluation
# failed: memory exhausted"). Force serial execution instead.
BiocParallel::register(BiocParallel::SerialParam())

slug <- "remap_tf_activity_gsea"
out_dir <- fig_dir(slug)
hs <- org.Hs.eg.db

# --- Data prep: actividad x union a cada TF (cualquier linea celular) -----

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

binary_df <- remap_hits %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, mean) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, mean), ~ replace_na(.x, FALSE)))

# TFs con >100 promotores unidos en AMBAS replicas, excluyendo marcas de
# histona (nombres tipo "H3K4me3", no son TFs puntuales).
nTF <- data_TF %>%
  group_by(rep) %>%
  summarise(across(where(is.logical), sum, na.rm = TRUE)) %>%
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
repname <- levels(data_TF$rep)

# --- Wilcoxon: efecto de cada TF sobre la actividad media -----------------

wilcox <- map(data_TF %>% select(all_of(nTF$name)), ~ coin::wilcox_test(formula = mean ~ .x | rep, data = data_TF) %>% pvalue())
wilcox <- tibble(feature = names(wilcox), pval = list_c(wilcox), pval_corr = p.adjust(pval, "BH", length(wilcox)))

wilcox_by_rep <- function(rep_id) {
  map(data_TF %>% filter(rep == rep_id) %>% select(all_of(nTF$name)), ~ coin::wilcox_test(formula = mean ~ .x, data = data_TF %>% filter(rep == rep_id), conf.int = TRUE))
}
wilcox_rep1 <- wilcox_by_rep(repname[1])
wilcox_rep2 <- wilcox_by_rep(repname[2])

wilcox_df <- map2(
  list(wilcox_rep1, wilcox_rep2), repname,
  ~ tibble(
    feature = names(.x),
    estimate = map(.x, ~ confint(.x)$estimate) %>% list_c(),
    P2.5 = map(.x, ~ confint(.x)$conf.int[1]) %>% list_c(),
    P97.5 = map(.x, ~ confint(.x)$conf.int[2]) %>% list_c(),
    rep = .y
  ) %>%
    pivot_longer(c(starts_with("estimate"), starts_with("P2.5"), starts_with("P97.5")), names_to = "val", values_to = "estimate") %>%
    arrange(desc(estimate)) %>%
    mutate(feature = fct_inorder(feature))
) %>%
  list_rbind() %>%
  left_join(wilcox, by = "feature")

write_tsv(wilcox_df, file.path(out_dir, "remap_tf_activity_wilcoxon.tsv"))

# --- Plots: efecto de union sobre la actividad -----------------------------

wilcox_wide <- wilcox_df %>%
  filter(pval_corr < 0.05) %>%
  mutate(act = ifelse(sign(estimate) == 1, "Actividad alta", "Actividad baja")) %>%
  group_by(feature) %>%
  mutate(n_dir = length(unique(act))) %>%
  filter(n_dir == 1) %>%
  ungroup() %>%
  pivot_wider(names_from = val, values_from = estimate)

panel_all <- wilcox_wide %>%
  arrange(desc(rep), estimate) %>%
  ggplot(aes(x = estimate, y = fct_inorder(feature), group = fct_inorder(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = act, alpha = rep)) +
  scale_alpha_manual(values = setNames(c(1, 0.5), repname)) +
  ggpubr::theme_pubr() +
  theme(text = element_text(size = 20), legend.position = "top", axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(fill = "Efecto", x = "Efecto sobre la actividad", y = "Factores de transcripción (ReMap)", alpha = "") +
  scale_fill_manual(values = c("Actividad alta" = "#D6741F", "Actividad baja" = "#7FB800"))
ggsave(file.path(out_dir, "remap_act.jpg"), panel_all, width = 9, height = 6.75, units = "in")

panel_low <- wilcox_wide %>%
  filter(act == "Actividad baja") %>%
  arrange(desc(rep), estimate) %>%
  ggplot(aes(x = estimate, y = fct_inorder(feature), group = fct_inorder(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = act, alpha = rep)) +
  scale_alpha_manual(values = setNames(c(1, 0.5), repname)) +
  geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5) +
  ggpubr::theme_pubr() +
  theme(text = element_text(size = 20), legend.position = "top") +
  labs(fill = "Efecto", x = "Efecto sobre la actividad", y = "Factores de transcripción (ReMap)", alpha = "") +
  scale_fill_manual(values = c("Actividad baja" = "#7FB800"))
ggsave(file.path(out_dir, "remap_act_low.jpg"), panel_low, width = 9, height = 8, units = "in")

# --- Scatter del TF con efecto mas negativo (analogo a los de R8) ---------
# Union a PHF19 (subunidad de PRC2) vs. actividad, en bins de 100
# promotores ordenados por rango - mismo formato que los scatters de
# R/16_tf_contexto_endogeno para NFYA/SP1/SP2.

top_neg_tf <- wilcox_df %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(min_estimate = min(estimate), .groups = "drop") %>%
  slice_min(min_estimate, n = 1) %>%
  pull(feature)

data_TF_top_neg <- add_mean_sw_bins(data %>% select(seq_id, rep, mean)) %>%
  left_join(remap_hits %>% filter(TF == top_neg_tf) %>% distinct(seq_id) %>% mutate(bound = TRUE), by = "seq_id") %>%
  mutate(bound = replace_na(bound, FALSE))

ggsave(
  file.path(out_dir, paste0(top_neg_tf, "_scatter.jpg")),
  tf_scatter(data_TF_top_neg, "bound", paste0("Unión a ", top_neg_tf, " (mayor efecto negativo sobre actividad)")),
  width = 9, height = 6.75, units = "in"
)

# --- GSEA sobre la lista de TFs ordenada por efecto en actividad ----------

pregsea_act <- function(rep_id) {
  feats <- wilcox_df %>% filter(rep == rep_id, val == "estimate") %>% pull(feature) %>% unique()
  ids <- AnnotationDbi::select(hs, keys = feats, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
    distinct()
  pregsea <- wilcox_df %>%
    filter(rep == rep_id, val == "estimate") %>%
    arrange(desc(estimate)) %>%
    left_join(ids, by = c("feature" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID))
  ordered_genes <- pregsea$estimate
  names(ordered_genes) <- pregsea$ENTREZID
  ordered_genes
}

gse_go <- map(repname, ~ gseGO(pregsea_act(.x), ont = "all", OrgDb = "org.Hs.eg.db"))
names(gse_go) <- repname
gse_go_df <- map(gse_go, ~ as.data.frame(.x) %>% arrange(desc(NES)))

write_xlsx(gse_go_df, file.path(out_dir, "GSEA_activity_remap.xlsx"))

# --- Panel: top pathway (mayor NES, Rep 1) + TFs que lo integran ----------

top_id <- gse_go_df[[repname[1]]]$ID[1]
top_desc <- gse_go_df[[repname[1]]]$Description[1]
top_entrez <- gse_go[[repname[1]]]@geneSets[[top_id]]
feats_rep1 <- wilcox_df %>% filter(rep == repname[1], val == "estimate") %>% pull(feature) %>% unique()
ids_rep1 <- AnnotationDbi::select(hs, keys = feats_rep1, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% distinct()
top_tf <- ids_rep1 %>% filter(ENTREZID %in% top_entrez) %>% pull(SYMBOL) %>% unique()

p1 <- enrichplot::gseaplot2(gse_go[[repname[1]]], top_id, title = top_desc, subplots = 1, color = "darkred")

p2 <- wilcox_wide %>%
  filter(rep == repname[1]) %>%
  mutate(act = ifelse(feature %in% top_tf, "set", act)) %>%
  arrange(desc(estimate)) %>%
  ggplot(aes(x = estimate, y = fct_inorder(feature))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = act, alpha = act)) +
  ggpubr::theme_pubr() +
  theme(text = element_text(size = 10), legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(fill = "Efecto", x = "Efecto sobre la actividad", y = "Factores de transcripción") +
  scale_fill_manual(values = c("Actividad alta" = "#D6741F", "Actividad baja" = "#7FB800", "set" = "darkred")) +
  scale_alpha_manual(values = c("Actividad alta" = 0.4, "Actividad baja" = 0.4, "set" = 1)) +
  coord_flip()

combined <- aplot::gglist(list(p1, p2), ncol = 1, heights = c(0.75, 0.5))
ggsave(file.path(out_dir, "GSEA_top_pathway_rep1.jpg"), combined, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
