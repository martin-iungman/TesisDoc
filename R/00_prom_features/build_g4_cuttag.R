# Senal de G4 plegado (G4-CUT&Tag con Flag-BG4, HEK293T) por fragmento FP
# de la library. GSE178668 (PMID 34400476): BG4 rep1 (GSM5395699), BG4
# rep2 (GSM5401689) e IgG (GSM5395700); bigWig de cobertura en RPM sobre
# hg38. GEO no tiene archivos de picos, asi que se usa la cobertura media
# (ponderada por pb) en los 252 pb del fragmento y el enriquecimiento
# log2((BG4 + PSEUDO) / (IgG + PSEUDO)) de cada replica. CUT&Tag no tiene
# informacion de hebra. Es senal del locus endogeno, no del plasmido.
#
# Tambien la cobertura media de ATAC-seq de HEK293T de la misma serie
# (GSM5395701), para controlar la accesibilidad.
#
# Salida: data/processed/g4_cuttag_hek293t_v1.tsv (seq_id, media de cada
# bigWig y log2 de enriquecimiento por replica).
# Requiere data/raw/library.bed y los bigWig de heavy_data_paths.R.
# Run from the TesisDoc repo root.

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
source("R/00_prom_features/heavy_data_paths.R")

PSEUDO <- 0.1 # RPM
out_path <- "data/processed/g4_cuttag_hek293t_v1.tsv"

lib <- rtracklayer::import.bed("data/raw/library.bed")
lib_fp <- lib[startsWith(lib$name, "FP")]
stopifnot(all(width(lib_fp) == 252L))
strand(lib_fp) <- "*"

# cobertura media por pb en cada region (pb sin dato = 0)
bw_mean <- function(path, regions) {
  cov <- rtracklayer::import(BigWigFile(path), which = regions, as = "RleList")
  vals <- numeric(length(regions))
  for (chr in unique(as.character(seqnames(regions)))) {
    i <- which(as.character(seqnames(regions)) == chr)
    if (!chr %in% names(cov)) next
    vals[i] <- mean(Views(cov[[chr]], ranges(regions)[i]))
  }
  vals
}

res <- tibble(
  seq_id = lib_fp$name,
  bg4_rep1 = bw_mean(path_g4_cuttag_bg4_r1, lib_fp),
  bg4_rep2 = bw_mean(path_g4_cuttag_bg4_r2, lib_fp),
  igg = bw_mean(path_g4_cuttag_igg, lib_fp),
  atac = bw_mean(path_atac_hek293t, lib_fp)
) %>%
  mutate(
    log2enr_rep1 = log2((bg4_rep1 + PSEUDO) / (igg + PSEUDO)),
    log2enr_rep2 = log2((bg4_rep2 + PSEUDO) / (igg + PSEUDO))
  )

print(summary(res[, -1]))
print(cor(select(res, log2enr_rep1, log2enr_rep2, bg4_rep1, igg, atac), method = "spearman"))
write_tsv(res, out_path)
message("Guardado ", out_path)
