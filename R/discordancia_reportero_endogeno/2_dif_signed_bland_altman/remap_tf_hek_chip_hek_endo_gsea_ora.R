# GSEA + ORA for the fully matched-context TF screen
# (remap_tf_hek_chip_hek_endo.R, same folder): ChIP-seq restricted to
# HEK293/HEK293T, endogenous side also HEK293 (hek_tpm). Same design as
# remap_tf.R / remap_tf_hek_endo_gsea_ora.R. ORA over the sobreestima set
# is expected to be underpowered here (n=10, thematically mixed - see
# session 2026-09-11 discussion) but cheap to check.

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(writexl)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
hs <- org.Hs.eg.db

resultados_lmm <- read_tsv(file.path(out_dir, "remap_tf_rank_dif_lmm_hek_chip_hek_endo.tsv"), show_col_types = FALSE)

ids <- AnnotationDbi::select(hs, keys = resultados_lmm$variable, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct()

# --- GSEA ------------------------------------------------------------------

pregsea <- resultados_lmm %>%
  arrange(desc(estimate)) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))
ordered_genes <- pregsea$estimate
names(ordered_genes) <- pregsea$ENTREZID

gse_go <- gseGO(ordered_genes, ont = "all", OrgDb = "org.Hs.eg.db")
gse_go_df <- as.data.frame(gse_go) %>% arrange(desc(NES))
write_xlsx(list(GSEA = gse_go_df), file.path(out_dir, "GSEA_remap_rank_dif_hek_chip_hek_endo.xlsx"))
message(sum(gse_go_df$p.adjust < 0.05), " terminos GO significativos (GSEA, p.adjust<0.05) de ", nrow(gse_go_df), " testeados")

# --- ORA sobre el set sobreestima -------------------------------------------

pos_tf_entrez <- resultados_lmm %>%
  filter(p_adj < 0.05, estimate > 0) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID)

message(length(pos_tf_entrez), " de los 10 TFs sobreestima mapean a ENTREZID")

ora_pos <- enrichGO(
  gene = pos_tf_entrez,
  universe = ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "all",
  pAdjustMethod = "BH"
)
ora_pos_df <- as.data.frame(ora_pos) %>% arrange(p.adjust)
write_xlsx(list(ORA_sobreestima = ora_pos_df), file.path(out_dir, "ORA_remap_rank_dif_hek_chip_hek_endo_sobreestima.xlsx"))
message(nrow(ora_pos_df), " terminos GO significativos (ORA, BH) en el set sobreestima")
if (nrow(ora_pos_df) > 0) print(ora_pos_df %>% select(Description, ONTOLOGY, Count, p.adjust) %>% head(10))

message("Guardado en ", out_dir)
