# Shared intermediate: per-promoter x per-TF ReMap2022 ChIP-seq peak
# overlaps (long format - one row per seq_id/TF/sample combination),
# built once from the raw non-redundant peaks bed (5.2GB) so downstream
# figure scripts (R/16_tf_contexto_endogeno and the ReMap GSEA figure)
# don't each have to re-read it.
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("ReMap dataset" section) and TF_remap_funciona.R ("EL POSTA" section).
# Run from the TesisDoc repo root. Output: data/processed/remap_tf_hits.tsv

library(rtracklayer)
library(tidyverse)

source("R/00_prom_features/heavy_data_paths.R")

remap <- rtracklayer::import.bed(path_remap_nr)
lib <- rtracklayer::import.bed("data/raw/library.bed")
remap_lib <- plyranges::join_overlap_inner(remap, lib)

remap_tf_hits <- values(remap_lib) %>%
  as_tibble() %>%
  filter(str_detect(name.y, "^FP.{6}_")) %>%
  separate(name.x, into = c("TF", "sample"), sep = ":") %>%
  select(seq_id = name.y, TF, sample) %>%
  distinct()

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(remap_tf_hits, "data/processed/remap_tf_hits.tsv")
message("Wrote data/processed/remap_tf_hits.tsv with ", nrow(remap_tf_hits), " rows")
