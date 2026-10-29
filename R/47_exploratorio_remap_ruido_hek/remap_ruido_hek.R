# Exploratorio de remap_tf_noise_gsea (R/21_remap_noise_gsea) usando solo
# los picos de ChIP-seq de ReMap2022 en HEK293 (str_detect "HEK293":
# HEK293, HEK293T, HEK293-FT, igual que tf_especificidad.R), en vez de
# cualquier linea celular. Mismo analisis: AUC de la union de cada TF
# prediciendo var_rank_sw alto (p de Wald con varianza de DeLong, BH sobre
# TFs x replicas), mismo criterio de TFs (>100 promotores unidos en ambas
# replicas, sin histonas) y GSEA (GO Cellular Component) sobre la lista
# ordenada por AUC. Ademas compara la AUC de cada TF con la del analisis
# con todas las lineas.
#
# Salida en figures/exploratorio_remap_ruido_hek/ (no es una figura
# numerada, solo jpg). Requiere data/processed/remap_tf_hits.tsv,
# activity_stats_highconf.tsv y prom_df.tsv, y la tabla
# remap_tf_noise_auc.tsv de remap_tf_noise_gsea para la comparacion.
# Run from the TesisDoc repo root.

library(tidyverse)
library(pROC)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(writexl)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/exploratorio_remap_ruido_hek"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hs <- org.Hs.eg.db

# --- Data prep: ruido x union a cada TF en HEK293 ---------------------------

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE) %>%
  filter(str_detect(sample, "HEK293"))
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>% add_noise_rank()

binary_df <- remap_hits %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, var_rank_sw) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, var_rank_sw), ~ replace_na(.x, FALSE)))

vbles <- data_TF %>%
  group_by(rep) %>%
  summarise(across(where(is.logical), sum)) %>%
  pivot_longer(-rep) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2, !str_detect(name, "^H[0-9]")) %>%
  distinct(name) %>%
  pull(name)
repname <- unique(data_TF$rep)
message(length(vbles), " TFs con >100 promotores unidos en HEK293 en ambas replicas")

data_TF <- data_TF %>% mutate(across(all_of(vbles), ~ factor(.x, levels = c("TRUE", "FALSE"))))

# --- AUC de cada TF prediciendo ruido alto --------------------------------

auc_df <- map(repname, function(r) {
  df <- data_TF %>% filter(rep == r)
  map(vbles, function(f) {
    roc_obj <- pROC::roc(df[[f]], df[["var_rank_sw"]], ci = TRUE, direction = ">", quiet = TRUE)
    auc_val <- as.numeric(roc_obj$auc)
    se_delong <- sqrt(pROC::var(roc_obj, method = "delong"))
    tibble(
      feature = f, rep = r, n_unidos = sum(df[[f]] == "TRUE"), AUC = auc_val,
      ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]),
      pval = 2 * pnorm(-abs((auc_val - 0.5) / se_delong))
    )
  }) %>% list_rbind()
}) %>% list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, method = "BH"))
write_tsv(auc_df, file.path(out_dir, "remap_hek_noise_auc.tsv"))

consistentes <- auc_df %>%
  filter(pval_corr < 0.05) %>%
  group_by(feature) %>%
  filter(n() == 2, n_distinct(AUC > 0.5) == 1) %>%
  summarise(AUC = mean(AUC), .groups = "drop") %>%
  arrange(AUC)
message("Consistentes: ", sum(consistentes$AUC > 0.5), " con ruido alto, ", sum(consistentes$AUC < 0.5), " con ruido bajo")
message("Ruido bajo (top 15): ", paste(head(consistentes$feature[consistentes$AUC < 0.5], 15), collapse = ", "))
message("Ruido alto (top 15): ", paste(head(rev(consistentes$feature[consistentes$AUC > 0.5]), 15), collapse = ", "))

panel_auc <- plot_noise_auc_summary(auc_df, base_size = 16) +
  labs(y = "Factores de transcripción (ReMap, HEK293)")
ggsave(file.path(out_dir, "remap_hek_noise.jpg"), panel_auc, width = 12, height = max(6.75, 2 + 0.25 * nrow(consistentes)), units = "in", limitsize = FALSE)

# --- Comparacion con todas las lineas --------------------------------------

auc_todas <- read_tsv(file.path(fig_dir("remap_tf_noise_gsea"), "remap_tf_noise_auc.tsv"), show_col_types = FALSE)
comp <- inner_join(auc_df, auc_todas, by = c("feature", "rep"), suffix = c("_hek", "_todas"))
cors <- comp %>% group_by(rep) %>% summarise(rho = cor(AUC_hek, AUC_todas, method = "spearman"), n = n())
print(cors)
p_comp <- comp %>%
  ggplot(aes(AUC_todas, AUC_hek)) +
  geom_hline(yintercept = 0.5, col = "grey60") +
  geom_vline(xintercept = 0.5, col = "grey60") +
  geom_abline(linetype = "dashed", col = "grey40") +
  geom_point(alpha = 0.5) +
  facet_wrap(~rep) +
  geom_text(data = cors, aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f\nn = %d", rho, n)), hjust = -0.1, vjust = 1.2, size = 5) +
  ggpubr::theme_pubr(base_size = 16) +
  labs(x = "AUC (picos en cualquier línea)", y = "AUC (solo picos en HEK293)")
ggsave(file.path(out_dir, "auc_hek_vs_todas.jpg"), p_comp, width = 11, height = 6, units = "in")

# --- GSEA (GO Cellular Component) -----------------------------------------

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
gse_go <- map(repname, function(r) {
  pregsea <- pregsea_noise(r)
  gseGO(setNames(pregsea$AUC, pregsea$ENTREZID), ont = "CC", OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.1, seed = TRUE)
})
names(gse_go) <- repname
gse_go_df <- map(gse_go, ~ as.data.frame(.x) %>% arrange(p.adjust))
write_xlsx(gse_go_df, file.path(out_dir, "GSEA_remap_hek_noise.xlsx"))
walk(repname, function(r) {
  message("GSEA ", r, " (", nrow(gse_go_df[[r]]), " terminos con p.adjust < 0.1):")
  print(head(select(gse_go_df[[r]], ID, Description, setSize, NES, p.adjust), 15))
})

# MLL1/2 complex (el termino destacado en remap_tf_noise_gsea), si entra
mll_id <- "GO:0044665"
for (r in repname) {
  if (mll_id %in% gse_go_df[[r]]$ID) {
    p <- enrichplot::gseaplot2(gse_go[[r]], mll_id, title = paste0(gse_go_df[[r]]$Description[gse_go_df[[r]]$ID == mll_id], " - ", r), subplots = 1:2, base_size = 16)
    ggsave(file.path(out_dir, paste0("GSEA_MLL1_2_", str_replace(r, " ", ""), ".jpg")), p, width = 10, height = 7, units = "in")
  } else {
    message("MLL1/2 complex no significativo en ", r)
  }
}

message("Exploratorio guardado en ", out_dir)
