# Resume por promotor la divergencia de Kimura de las copias de TE que lo
# solapan (data/processed/te_hits_kimura.tsv, generado por
# R/00_prom_features/filter_rmsk_align.py a partir de hg38.fa.align.gz).
# Misma logica que build_te_age.R: TE dominante (mayor superposicion con el
# fragmento de 252 pb) y promedio ponderado por superposicion, para
# kimura_cpg (K2P con correccion CpG) y kimura_nocpg (K2P sin correccion).
#
# Output: data/processed/te_age_kimura.tsv (seq_id, TE_name_k,
# TE_class_k, kimura_cpg_dom, kimura_nocpg_dom, kimura_cpg_wmean,
# kimura_nocpg_wmean, perc_div_hdr_dom). Las copias con K2P saturada (NA)
# se descartan del promedio. Run from the TesisDoc repo root.

library(tidyverse)

hits <- read_tsv("data/processed/te_hits_kimura.tsv", show_col_types = FALSE, na = "NA")
message("Filas (copia x promotor): ", nrow(hits), " | promotores: ", n_distinct(hits$seq_id))

te_age_kimura <- hits %>%
  group_by(seq_id) %>%
  summarise(
    dom = which.max(overlap),
    TE_name_k = name[dom], TE_class_k = class[dom],
    kimura_cpg_dom = kimura_cpg[dom], kimura_nocpg_dom = kimura_nocpg[dom],
    perc_div_hdr_dom = perc_div_hdr[dom],
    kimura_cpg_wmean = weighted.mean(kimura_cpg, overlap, na.rm = TRUE),
    kimura_nocpg_wmean = weighted.mean(kimura_nocpg, overlap, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(-dom)

write_tsv(te_age_kimura, "data/processed/te_age_kimura.tsv")
message("Wrote data/processed/te_age_kimura.tsv with ", nrow(te_age_kimura), " promoters")
print(summary(te_age_kimura %>% select(starts_with("kimura"), perc_div_hdr_dom)))
