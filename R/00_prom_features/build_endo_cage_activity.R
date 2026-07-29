# Shared intermediate: actividad CAGE endogena (FANTOM5, procesada con
# CAGEr) para HeLa y musculo esqueletico, analogo a hek_tpm en prom_df.tsv
# (que ya viene de build_prom_features.R, ver esa seccion "Endogenous
# HEK293 activity"). Separado en su propio output en vez de meterlo en
# prom_df.tsv para no tener que recorrer todo ese pipeline (~10 min) por
# 4 columnas nuevas.
#
# Ported from transcriptional_library/Analysis/scripts/cell_line_cage.qmd
# ("HeLa CAGE data" y "skeletal muscle" secciones). Se uso el patron de
# hek_cage.qmd (sum(score.y, ...) tras el join) en vez del que tenia el
# bloque de HeLa en el qmd original, que sumaba una columna "hela_tpm"
# que en ese punto todavia no existe (bug - probablemente una copia
# incompleta del patron ya usado, correctamente, en el bloque de
# musculo mas abajo en el mismo script).
#
# Requiere: data/raw/library.bed. Run from the TesisDoc repo root.
# Output: data/processed/endo_cage_activity.tsv (seq_id, hela_tpm,
# hela_interq_width, muscle_tpm, muscle_interq_width).

library(CAGEr)
library(tidyverse)
source("R/00_prom_features/heavy_data_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

lib <- rtracklayer::import.bed("data/raw/library.bed")
lib$gene_name <- lib$name

process_cage <- function(files, sample_labels, merged_label, alpha) {
  cell_line <- CAGEr::CAGEexp(
    genomeName = "BSgenome.Hsapiens.UCSC.hg38",
    inputFiles = files,
    inputFilesType = "bedScore",
    sampleLabels = sample_labels
  )
  cell_line <- CAGEr::getCTSS(cell_line)
  cell_line <- CAGEr::mergeSamples(cell_line, mergeIndex = rep(1, length(sample_labels)), mergedSampleLabels = merged_label)
  cell_line <- annotateCTSS(cell_line, lib, upstream = 0, downstream = 252)
  cell_line <- CAGEr::normalizeTagCount(cell_line, method = "powerLaw", fitInRange = c(5, 1000), alpha = alpha, T = 1e6)
  cell_line <- filterLowExpCTSS(cell_line, thresholdIsTpm = TRUE, nrPassThreshold = 1, threshold = 0.5)
  cell_line <- paraclu(cell_line, keepSingletonsAbove = 1, nrCores = 5)
  cell_line <- cumulativeCTSSdistribution(cell_line, clusters = "tagClusters", useMulticore = TRUE)
  cell_line <- quantilePositions(cell_line, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
  tagClustersGR(cell_line, 1, qLow = 0.1, qUp = 0.9) %>%
    plyranges::join_overlap_left_directed(lib, .) %>%
    as_tibble()
}

# --- HeLa (3 replicas ENCODE, FANTOM5) -------------------------------------

hela_df <- process_cage(path_fantom5_hela_raw, paste0("Hela_", 1:3), "Hela", alpha = 1.14)
hela_summ <- hela_df %>%
  group_by(name) %>%
  summarise(hela_tpm = sum(score.y, na.rm = TRUE), hela_interq_width = sum(interquantile_width))

# --- Musculo esqueletico (4 muestras de tejido, FANTOM5) -------------------

tissue_ontology <- readxl::read_xls(path_fantom5_tissue_ontology) %>% rename(sample_id = `Sample ID`)
tissue_files <- list.files(path_fantom5_tissue_dir, full.names = TRUE, pattern = ".nobarcode.ctss.bed.gz$")
tissue_lookup <- tibble(filename = tissue_files, sample_id = str_extract(tissue_files, "CNhs.{5}")) %>%
  inner_join(tissue_ontology, by = "sample_id")
muscle_samples <- tissue_lookup %>% filter(`Facet ontology term` == "skeletal muscle tissue")

muscle_df <- process_cage(muscle_samples$filename, muscle_samples$sample_id, "muscle", alpha = 1.05)
muscle_summ <- muscle_df %>%
  group_by(name) %>%
  summarise(muscle_tpm = sum(score.y, na.rm = TRUE), muscle_interq_width = sum(interquantile_width))

# --- Output -----------------------------------------------------------

endo_cage_activity <- full_join(hela_summ, muscle_summ, by = "name") %>% rename(seq_id = name)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(endo_cage_activity, "data/processed/endo_cage_activity.tsv")
message("Wrote data/processed/endo_cage_activity.tsv with ", nrow(endo_cage_activity), " rows")
