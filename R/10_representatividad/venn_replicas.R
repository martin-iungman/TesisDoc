# representatividad_library_venn (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Diagrama de Venn de secuencias detectadas por replica vs. la library
# completa (sin filtros de calidad - "cuanto de la library aparece
# secuenciado, aunque sea con 1 read").
#
# Requiere: data/processed/data_long.tsv (ver
# R/01_activity_stats/build_activity_stats.R) y data/raw/seq_data.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
source("R/functions/fig_paths.R")

slug <- "representatividad_library_venn"
out_dir <- fig_dir(slug)

data_long <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE)
all_seq <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE) %>% filter(!seq_id %in% dupl_id)
spike_id <- c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")

VennDiagram::venn.diagram(
  list(
    "Rep 1" = data_long %>% filter(counts > 0, rep == "Rep 1", !seq_id %in% spike_id) %>% pull(seq_id) %>% unique(),
    "Rep 2" = data_long %>% filter(counts > 0, rep == "Rep 2", !seq_id %in% spike_id) %>% pull(seq_id) %>% unique(),
    "Library" = all_seq$seq_id
  ),
  file.path(out_dir, paste0(slug, ".png")),
  imagetype = "png",
  disable.logging = TRUE,
  alpha = c(.5, .5, .5),
  fill = c("#AADAD4", "#358AAA", "#dfedf6"),
  cat.cex = c(0, 0, 0)
)

message("Figura guardada en ", out_dir)
