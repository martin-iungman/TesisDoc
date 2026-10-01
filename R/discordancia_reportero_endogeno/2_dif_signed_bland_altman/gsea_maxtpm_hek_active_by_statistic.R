# Robustness check for remap_tf_maxtpm_hek_active_gsea_ora.R's GSEA
# (same folder): re-runs it ranking by the t-statistic (estimate/SE)
# instead of the raw estimate. GSEA doesn't use each gene's own p-value/
# significance at all - it ranks purely by whatever score you feed it,
# so TFs with few bound promoters (near the >100 threshold, wide CIs)
# get the same weight as well-powered ones if ranked by raw estimate.
# The t-statistic naturally down-weights imprecise estimates.

library(tidyverse)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

BiocParallel::register(BiocParallel::SerialParam())

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
hs <- org.Hs.eg.db

resultados_lmm <- read_tsv(file.path(out_dir, "remap_tf_rank_dif_lmm_maxtpm_hek_active.tsv"), show_col_types = FALSE)

# IMPORTANT: `statistic` in this table is NOT flipped the way `estimate`
# is (remap_tf_maxtpm_hek_active.R only negates estimate/conf.low/
# conf.high, not statistic) - it's still on the raw FALSE-vs-TRUE
# convention, i.e. opposite sign from the already-flipped `estimate`
# (checked directly: SETDB1 has estimate=+0.076 but statistic=-11.77).
# Negate it here so ranking direction matches estimate/p_adj.
ids <- AnnotationDbi::select(hs, keys = resultados_lmm$variable, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct()

pregsea <- resultados_lmm %>%
  mutate(statistic_flipped = -statistic) %>%
  arrange(desc(statistic_flipped)) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))
ordered_genes <- pregsea$statistic_flipped
names(ordered_genes) <- pregsea$ENTREZID

gse_go <- gseGO(ordered_genes, ont = "all", OrgDb = "org.Hs.eg.db")
gse_go_df <- as.data.frame(gse_go) %>% arrange(desc(NES))
write_xlsx(list(GSEA = gse_go_df), file.path(out_dir, "GSEA_remap_rank_dif_maxtpm_hek_active_by_statistic.xlsx"))

message(sum(gse_go_df$p.adjust < 0.05), " terminos GO significativos de ", nrow(gse_go_df), " testeados (ranking por t-estadistico)")
print(gse_go_df %>% filter(p.adjust < 0.05) %>% arrange(desc(abs(NES))) %>% select(Description, ONTOLOGY, NES, p.adjust, setSize), n = 40)

message("\nComparacion rapida con el ranking por estimate crudo (mismos terminos con p.adjust<0.05?):")
gse_go_estimate <- read_xlsx(file.path(out_dir, "GSEA_remap_rank_dif_maxtpm_hek_active.xlsx")) %>% filter(p.adjust < 0.05) %>% pull(Description)
gse_go_stat_sig <- gse_go_df %>% filter(p.adjust < 0.05) %>% pull(Description)
message("En ambos: ", length(intersect(gse_go_estimate, gse_go_stat_sig)))
message("Solo en ranking por estimate: ", length(setdiff(gse_go_estimate, gse_go_stat_sig)))
message("Solo en ranking por statistic: ", length(setdiff(gse_go_stat_sig, gse_go_estimate)))

message("Guardado en ", out_dir)
