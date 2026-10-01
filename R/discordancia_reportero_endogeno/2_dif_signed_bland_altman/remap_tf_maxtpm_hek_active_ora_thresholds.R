# ORA over the max_tpm + all ReMap ChIP + HEK293-active screen
# (remap_tf_maxtpm_hek_active.R, same folder), comparing FDR<0.05 vs.
# FDR<0.001 as the gene-set-of-interest threshold, for BOTH sobreestima
# and subestima. Companion to remap_tf_maxtpm_hek_active_gsea_ora.R
# (which only ran ORA on sobreestima @ FDR<0.05).

library(tidyverse)
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

ids <- AnnotationDbi::select(hs, keys = resultados_lmm$variable, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct()

run_ora <- function(fdr_cutoff, direction) {
  genes <- resultados_lmm %>%
    filter(p_adj < fdr_cutoff, if (direction == "sobreestima") estimate > 0 else estimate < 0) %>%
    left_join(ids, by = c("variable" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    pull(ENTREZID)

  n_genes <- length(genes)
  if (n_genes == 0) {
    return(list(n = 0, df = tibble()))
  }

  ora <- enrichGO(
    gene = genes,
    universe = ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "all",
    pAdjustMethod = "BH"
  )
  df <- as.data.frame(ora) %>% arrange(p.adjust) %>% as_tibble()
  list(n = n_genes, df = df)
}

combos <- expand_grid(direction = c("sobreestima", "subestima"), fdr = c(0.05, 0.001))

resultados <- combos %>%
  mutate(res = map2(fdr, direction, run_ora)) %>%
  mutate(n_genes = map_dbl(res, "n"), df = map(res, "df"), n_sig_terms = map_dbl(df, nrow)) %>%
  select(-res)

message("--- Resumen ---")
resultados %>%
  select(direction, fdr, n_genes, n_sig_terms) %>%
  print(n = 10)

write_xlsx(
  set_names(resultados$df, paste0(resultados$direction, "_fdr", resultados$fdr)),
  file.path(out_dir, "ORA_remap_rank_dif_maxtpm_hek_active_thresholds.xlsx")
)

walk2(resultados$direction, resultados$fdr, function(dir, fdr) {
  row <- resultados %>% filter(direction == dir, fdr == !!fdr)
  message("\n=== ", dir, " @ FDR<", fdr, " (n=", row$n_genes, " TFs, ", row$n_sig_terms, " terminos GO sig.) ===")
  if (nrow(row$df[[1]]) > 0) {
    print(row$df[[1]] %>% select(Description, ONTOLOGY, Count, p.adjust) %>% head(15))
  }
})

message("\nGuardado en ", out_dir)
