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

detected <- function(r) {
  data_long %>% filter(counts > 0, rep == r, !seq_id %in% spike_id) %>% pull(seq_id) %>% unique()
}
rep1 <- detected("Rep 1")
rep2 <- detected("Rep 2")
library_ids <- all_seq$seq_id

fuera <- setdiff(union(rep1, rep2), library_ids)
if (length(fuera) > 0) warning(length(fuera), " secuencias detectadas no estan en la library")

n <- c(
  rep1_solo = length(setdiff(rep1, rep2)),
  ambas = length(intersect(rep1, rep2)),
  rep2_solo = length(setdiff(rep2, rep1)),
  library_solo = length(setdiff(library_ids, union(rep1, rep2)))
)
print(n)

# Geometria fija (no proporcional al tamano de los conjuntos): library
# como circulo exterior y las dos replicas superpuestas adentro.
circle <- function(x0, y0, r, id) {
  t <- seq(0, 2 * pi, length.out = 361)
  tibble(x = x0 + r * cos(t), y = y0 + r * sin(t), id = id)
}
circles <- bind_rows(
  circle(0, 0, 1, "library"),
  circle(-0.325, -0.05, 0.525, "rep1"),
  circle(0.325, -0.05, 0.525, "rep2")
)
fills <- c(library = "#EFF6FB", rep1 = "#358AAA", rep2 = "#AADAD4")
alphas <- c(library = 1, rep1 = 0.3, rep2 = 0.5)

labels <- tribble(
  ~x,     ~y,     ~label,                         ~size,
  0,      0.79,   "Library",                      5,
  0,      0.66,   format(n[["library_solo"]]),    4.2,
  -0.47,  0.15,   "Réplica 1",                    5,
  -0.47,  -0.07,  format(n[["rep1_solo"]]),       4.2,
  0,      -0.07,  format(n[["ambas"]]),           4.2,
  0.59,   0.15,   "Réplica 2",                    5,
  0.59,   -0.07,  format(n[["rep2_solo"]]),       4.2
)

p <- ggplot() +
  map(names(fills), function(i) {
    geom_polygon(
      data = filter(circles, id == i), aes(x, y),
      fill = fills[[i]], alpha = alphas[[i]], colour = NA
    )
  }) +
  # bordes en una capa aparte para que la transparencia no los atenue
  geom_path(data = circles, aes(x, y, group = id), colour = "grey15", linewidth = 0.4) +
  geom_text(data = labels, aes(x, y, label = label, size = size), colour = "grey35") +
  scale_size_identity() +
  coord_equal(xlim = c(-1.02, 1.02), ylim = c(-1.02, 1.02), expand = FALSE) +
  theme_void()

for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0(slug, ".", ext)), p, width = 4.3, height = 4.3, units = "in", dpi = 600, bg = "white")
}

message("Figura guardada en ", out_dir)
