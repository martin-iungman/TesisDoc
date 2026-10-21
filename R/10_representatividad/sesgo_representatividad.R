# sesgo_representatividad (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Panel A: venn de comparacion con el reportero de splicing pre-sorting
#   (proyecto 2022) + test de Fisher de asociacion.
# Panel B: violin de G+C para promotores detectados vs no detectados.
# Panel C: sesgo de contenido G+C en las lecturas secuenciadas
#   (2 representaciones alternativas: lineal y de barras).
#
# Requiere: data/processed/prom_df.tsv, data/processed/activity_stats_full.tsv,
# data/processed/data_long.tsv (ver R/00_prom_features y R/01_activity_stats)
# y data/external/2022_presort_counts.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "sesgo_representatividad"
out_dir <- fig_dir(slug)
clrs <- ghibli::ghibli_palette(name = "PonyoMedium")

all_seq <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE) %>% filter(!seq_id %in% dupl_id)
data_long <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
full_stats <- read_tsv("data/processed/activity_stats_full.tsv", show_col_types = FALSE)
spike_id <- c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")

# --- Panel A: venn contra reportero de splicing 2022 (pre-sorting) -------

presort_2022 <- read_tsv("data/external/2022_presort_counts.tsv", show_col_types = FALSE)
venn_sets <- list(
  "Reportero de splicing\n(pre-sorting)" = presort_2022 %>%
    filter(original_count > 0, seqname %in% all_seq$seq_id) %>%
    pull(seqname),
  "Reportero\ntranscripcional" = data_long %>%
    filter(counts > 0, !seq_id %in% spike_id) %>%
    pull(seq_id) %>%
    unique(),
  "Library" = all_seq$seq_id
)

transcr <- venn_sets$`Reportero\ntranscripcional`
splicing <- venn_sets$`Reportero de splicing\n(pre-sorting)`
fuera <- setdiff(union(transcr, splicing), venn_sets$Library)
if (length(fuera) > 0) warning(length(fuera), " secuencias detectadas no estan en la library")
n_venn <- c(
  transcr_solo = length(setdiff(transcr, splicing)),
  ambas = length(intersect(transcr, splicing)),
  splicing_solo = length(setdiff(splicing, transcr)),
  library_solo = length(setdiff(venn_sets$Library, union(transcr, splicing)))
)
print(n_venn)

# Mismo estilo que representatividad_library_venn: geometria fija (no
# proporcional), library como circulo exterior.
circle <- function(x0, y0, r, id) {
  t <- seq(0, 2 * pi, length.out = 361)
  tibble(x = x0 + r * cos(t), y = y0 + r * sin(t), id = id)
}
circles <- bind_rows(
  circle(0, 0, 1, "library"),
  circle(-0.28, -0.03, 0.56, "transcr"),
  circle(0.28, -0.03, 0.56, "splicing")
)
fills <- c(library = "#EFF6FB", transcr = "#358AAA", splicing = "#AADAD4")
alphas <- c(library = 1, transcr = 0.3, splicing = 0.5)

venn_labels <- tribble(
  ~x,     ~y,     ~label,                             ~size,
  0,      0.80,   "Library",                          4.6,
  0,      0.69,   format(n_venn[["library_solo"]]),   4.2,
  -0.52,  0.19,   "MPRA\ntranscripción",              4.4,
  -0.52,  -0.02,  format(n_venn[["transcr_solo"]]),   4.2,
  0,      -0.02,  format(n_venn[["ambas"]]),          4.2,
  0.55,   0.19,   "MPRA\nsplicing",                   4.4,
  0.55,   -0.02,  format(n_venn[["splicing_solo"]]),  4.2
)

panel_a <- ggplot() +
  map(names(fills), function(i) {
    geom_polygon(
      data = filter(circles, id == i), aes(x, y),
      fill = fills[[i]], alpha = alphas[[i]], colour = NA
    )
  }) +
  # bordes en una capa aparte para que la transparencia no los atenue
  geom_path(data = circles, aes(x, y, group = id), colour = "grey15", linewidth = 0.4) +
  geom_text(data = venn_labels, aes(x, y, label = label, size = size), colour = "grey35", lineheight = 0.9) +
  scale_size_identity() +
  coord_equal(xlim = c(-1.02, 1.02), ylim = c(-1.02, 1.02), expand = FALSE) +
  theme_void()
for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0("panel_a_venn_presort2022.", ext)), panel_a, width = 4.3, height = 4.3, units = "in", dpi = 600, bg = "white")
}

fisher_result <- tibble(
  all = venn_sets$Library,
  splicing = venn_sets$Library %in% venn_sets$`Reportero de splicing\n(pre-sorting)`,
  transcr = venn_sets$Library %in% venn_sets$`Reportero\ntranscripcional`
) %>%
  select(-all) %>%
  table() %>%
  fisher.test()
message("Fisher test (splicing x transcripcional): p = ", format.pval(fisher_result$p.value), ", OR = ", round(fisher_result$estimate, 3))

# --- Panel C: sesgo de G+C en lecturas secuenciadas -----------------------

tmp <- bind_rows(prom_df %>% mutate(rep = "Rep 1"), prom_df %>% mutate(rep = "Rep 2")) %>%
  left_join(full_stats, by = c("seq_id", "rep")) %>%
  filter(!seq_id %in% spike_id) %>%
  mutate(sum_norm_counts = replace_na(sum_norm_counts, 0))
tmp$bin_gc <- as.numeric(as.character(cut_interval(tmp$g_c,
  length = 0.02,
  labels = seq(min(tmp$g_c) - .01, max(tmp$g_c) + .01, by = 0.02)
)))
tmp <- tmp %>%
  group_by(bin_gc, rep) %>%
  summarise(sum = sum(sum_norm_counts), n = n(), .groups = "drop")
tmp <- tmp %>%
  group_by(rep) %>%
  mutate(sum_rel = sum / sum(sum), n_rel = n / sum(n)) %>%
  select(sum_rel, bin_gc, n_rel, rep)

panel_c_lineal <- tmp %>%
  pivot_longer(c(sum_rel, n_rel), values_to = "counts", names_to = "data") %>%
  ggplot(aes(bin_gc, counts, col = data)) +
  geom_line() +
  facet_wrap(~rep) +
  xlab("Contenido de [G+C]") +
  ylab("Proporción de lecturas secuenciadas") +
  scale_color_manual(values = clrs[c(3, 5)], labels = c("Distribución uniforme", "Observada")) +
  theme_pubr() +
  theme(legend.title = element_blank(), legend.position = "top", panel.spacing = unit(3, "lines"))
ggsave(file.path(out_dir, "panel_c_gc_bias_lineal.jpg"), panel_c_lineal, width = 9, height = 6.75, units = "in")

panel_c_bars <- tmp %>%
  ggplot(aes(x = bin_gc)) +
  geom_col(aes(y = n_rel, fill = "Distribución uniforme")) +
  geom_col(aes(y = sum_rel, fill = "Distribución observada"), alpha = 0.5) +
  facet_wrap(~rep) +
  scale_fill_manual(
    values = c("Distribución uniforme" = "#AADAD4", "Distribución observada" = "#358AAA"),
    breaks = c("Distribución uniforme", "Distribución observada")
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = c(1, 0.5)))) +
  xlab("Contenido de [G+C]") +
  ylab("Proporción") +
  theme_pubr() +
  theme(
    text = element_text(size = 25),
    # separacion suficiente para que "1.0" y "0.2" de paneles vecinos no se pisen
    panel.spacing = unit(1, "in"),
    # leyenda adentro, arriba a la derecha (sobre el panel de Rep 2)
    legend.position = "inside",
    legend.position.inside = c(0.99, 0.97),
    legend.justification = c(1, 1),
    legend.title = element_blank(),
    legend.background = element_blank(),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("panel_c_gc_bias_bars.", ext)), panel_c_bars, width = 13, height = 6.75, units = "in")
}

# --- Panel B: violin de G+C, secuenciados vs no --------------------------

panel_b <- prom_df %>%
  select(seq_id, g_c) %>%
  unique() %>%
  left_join(full_stats, by = "seq_id") %>%
  ggplot(aes(!is.na(sum_counts), g_c)) +
  geom_violin(draw_quantiles = 0.5, fill = "#AADAD4") +
  scale_x_discrete(labels = c("FALSE" = "Promotores\nno detectados", "TRUE" = "Promotores\ndetectados")) +
  labs(x = NULL, y = "Contenido de [G+C]") +
  theme_pubr() +
  theme(text = element_text(size = 25))
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("panel_b_gc_bias_violin.", ext)), panel_b, width = 6, height = 6.75, units = "in")
}

message("Paneles guardados en ", out_dir, ". Panel C tiene 2 versiones alternativas (lineal/bars) - definir cual usar.")
