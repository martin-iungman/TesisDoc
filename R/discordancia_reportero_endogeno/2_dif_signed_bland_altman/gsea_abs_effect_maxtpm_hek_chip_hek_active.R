# GSEA over |estimate| (magnitude of the signed discordance effect,
# regardless of direction), for the max_tpm + HEK293/HEK293T-restricted
# ChIP + HEK293-active screen (remap_tf_maxtpm_hek_chip_hek_active.R,
# same folder). Companion to gsea_abs_effect_maxtpm_hek_active.R (the
# all-ChIP version).
#
# CAVEAT (session 2026-09-14): |estimate| is not the same fully-
# confounded "magnitude/dispersion" question established earlier in the
# project (avg_rank is already controlled for in the mean here), but the
# signed effect still carries SOME residual correlation with each TF's
# raw activity effect (same mechanism as the prom_df/remap_tf "vs.
# activity" comparison plots elsewhere in this folder) - so a GSEA on
# |estimate| risks partly re-detecting "TFs with big activity effects in
# general" rather than something specific to discordance. This script
# also reports that residual correlation as a sanity check.

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

resultados_lmm <- read_tsv(file.path(out_dir, "remap_tf_rank_dif_lmm_maxtpm_hek_chip_hek_active.tsv"), show_col_types = FALSE)

# --- Sanity check: |estimate| vs. raw activity effect magnitude -----------

wilcox_activity <- read_tsv("figures/R2.6_remap_tf_activity_gsea/remap_tf_activity_wilcoxon.tsv", show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  group_by(feature) %>%
  summarise(activity_estimate = mean(estimate), .groups = "drop")

chequeo <- resultados_lmm %>%
  inner_join(wilcox_activity, by = c("variable" = "feature")) %>%
  mutate(abs_estimate = abs(estimate), abs_activity = abs(activity_estimate))
r_val <- cor(chequeo$abs_estimate, chequeo$abs_activity, method = "spearman")
message("Correlacion (Spearman) entre |efecto sobre discordancia| y |efecto sobre actividad cruda|: r=", round(r_val, 3), " (n=", nrow(chequeo), ")")

# --- GSEA sobre |estimate| --------------------------------------------------

ids <- AnnotationDbi::select(hs, keys = resultados_lmm$variable, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct()

pregsea <- resultados_lmm %>%
  mutate(abs_estimate = abs(estimate)) %>%
  arrange(desc(abs_estimate)) %>%
  left_join(ids, by = c("variable" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))
ordered_genes <- pregsea$abs_estimate
names(ordered_genes) <- pregsea$ENTREZID

gse_go <- gseGO(ordered_genes, ont = "all", OrgDb = "org.Hs.eg.db", scoreType = "pos")
gse_go_df <- as.data.frame(gse_go) %>% arrange(desc(NES))
write_xlsx(list(GSEA = gse_go_df), file.path(out_dir, "GSEA_remap_rank_dif_maxtpm_hek_chip_hek_active_abs.xlsx"))

message(sum(gse_go_df$p.adjust < 0.05), " terminos GO significativos de ", nrow(gse_go_df), " testeados (ranking por |estimate|)")

gse_sig <- gse_go_df %>% filter(p.adjust < 0.05) %>% arrange(desc(NES)) %>% as_tibble()
print(gse_sig %>% select(Description, ONTOLOGY, NES, p.adjust, setSize), n = 40)

if (nrow(gse_sig) > 0) {
  p_terms <- gse_sig %>%
    slice_max(order_by = NES, n = 20) %>%
    mutate(Description = fct_reorder(Description, NES)) %>%
    ggplot(aes(x = NES, y = Description)) +
    geom_col(fill = "#358AAA") +
    ggpubr::theme_pubr(base_size = 12) +
    labs(x = "NES", y = NULL, title = "GSEA sobre |efecto| (magnitud, sin signo) - max_tpm + ChIP HEK293/HEK293T, solo HEK-activos") +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "GSEA_terminos_enriquecidos_remap_rank_dif_maxtpm_hek_chip_hek_active_abs.jpg"), p_terms, width = 9, height = 8, units = "in")
}

message("Guardado en ", out_dir)
