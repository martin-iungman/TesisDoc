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
source("R/functions/venn_library.R")

slug <- "representatividad_library_venn"
out_dir <- fig_dir(slug)

data_long <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE)
all_seq <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE) %>% filter(!seq_id %in% dupl_id)
spike_id <- c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")

detected <- function(r) {
  data_long %>% filter(counts > 0, rep == r, !seq_id %in% spike_id) %>% pull(seq_id) %>% unique()
}

p <- venn_library(detected("Rep 1"), detected("Rep 2"), all_seq$seq_id, "Réplica 1", "Réplica 2")

for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0(slug, ".", ext)), p, width = 4.3, height = 4.3, units = "in", dpi = 600, bg = "white")
}

message("Figura guardada en ", out_dir)
