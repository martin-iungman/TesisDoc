# GSEA + ORA for the hek_tpm-based TF screen (remap_tf_hek_endo.R, same
# folder) - same design as remap_tf.R's GSEA/ORA sections, applied to
# remap_tf_rank_dif_lmm_hek_endo_excl_inactive.tsv (HEK293-active
# promoters only) instead of the max_tpm version. GSEA over the full
# ranked list (by estimate); ORA over the sobreestima set (estimate>0,
# p_adj<0.05, n=32 here vs n=6 in the max_tpm version - better powered).
#
# Does not refit the 954 LMMs - reads the already-computed TSV.

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(writexl)

# clusterProfiler/AnnotationDbi define S4 generics that shadow
# dplyr::select/filter/rename - same guard as remap_tf.R.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
hs <- org.Hs.eg.db

resultados_lmm <- read_tsv(file.path(out_dir, "remap_tf_rank_dif_lmm_hek_endo_excl_inactive.tsv"), show_col_types = FALSE)

ids <- AnnotationDbi::select(hs, keys = resultados_lmm$variable, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct()

# --- GSEA sobre el ranking completo por estimado ---------------------------

pregsea <- resultados_lmm %>%
  arrange(desc(estimate)) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))
ordered_genes <- pregsea$estimate
names(ordered_genes) <- pregsea$ENTREZID

gse_go <- gseGO(ordered_genes, ont = "all", OrgDb = "org.Hs.eg.db")
gse_go_df <- as.data.frame(gse_go) %>% arrange(desc(NES))
write_xlsx(list(GSEA = gse_go_df), file.path(out_dir, "GSEA_remap_rank_dif_hek_endo.xlsx"))

gse_go_sig <- gse_go_df %>% filter(p.adjust < 0.05)
if (nrow(gse_go_sig) > 0) {
  p_terms <- gse_go_sig %>%
    slice_max(order_by = abs(NES), n = 20) %>%
    mutate(
      dif = ifelse(NES > 0, "Reportero sobreestima", "Reportero subestima"),
      Description = fct_reorder(Description, NES)
    ) %>%
    ggplot(aes(x = NES, y = Description, fill = dif)) +
    geom_col() +
    ggpubr::theme_pubr(base_size = 12) +
    theme(legend.position = "top") +
    labs(fill = "Efecto", x = "NES", y = NULL, title = "Términos GO enriquecidos (GSEA, efecto de unión de TFs vs. HEK293)") +
    scale_fill_manual(values = c("Reportero sobreestima" = "#D6741F", "Reportero subestima" = "#7FB800"))
  ggsave(file.path(out_dir, "GSEA_terminos_enriquecidos_remap_rank_dif_hek_endo.jpg"), p_terms, width = 9, height = 8, units = "in")
  message(nrow(gse_go_sig), " terminos GO significativos (GSEA, p.adjust<0.05)")
} else {
  message("GSEA no devolvio terminos significativos - ver ", file.path(out_dir, "GSEA_remap_rank_dif_hek_endo.xlsx"))
}

# --- ORA sobre el set sobreestima (estimate>0, p_adj<0.05) -----------------

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
write_xlsx(list(ORA_sobreestima = ora_pos_df), file.path(out_dir, "ORA_remap_rank_dif_hek_endo_sobreestima.xlsx"))

if (nrow(ora_pos_df) > 0) {
  p_ora <- ora_pos_df %>%
    slice_min(order_by = p.adjust, n = 20) %>%
    mutate(Description = fct_reorder(Description, -p.adjust)) %>%
    ggplot(aes(x = -log10(p.adjust), y = Description, fill = Count)) +
    geom_col() +
    ggpubr::theme_pubr(base_size = 12) +
    labs(
      x = "-log10(p ajustado)", y = NULL, fill = "N° de TFs",
      title = paste0("Enriquecimiento GO (ORA) en TFs donde el reportero sobreestima\nvs. HEK293 (n=", length(pos_tf_entrez), ")")
    ) +
    scale_fill_gradient(low = "#7FB800", high = "#0D2C54")
  ggsave(file.path(out_dir, "ORA_remap_rank_dif_hek_endo_sobreestima.jpg"), p_ora, width = 9, height = 8, units = "in")
  message(nrow(ora_pos_df), " terminos GO significativos (ORA, BH) en el set sobreestima (n=", length(pos_tf_entrez), " TFs)")
} else {
  message("ORA no devolvio terminos significativos para el set sobreestima (n=", length(pos_tf_entrez), " TFs)")
}

message("Guardado en ", out_dir)
