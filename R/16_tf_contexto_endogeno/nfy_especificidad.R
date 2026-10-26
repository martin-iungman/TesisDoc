# nfy_especificidad (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Especificidad de NF-YA frente al resto de los factores de ReMap2022
# (universo "cualquier linea": 986 TFs con >100 promotores unidos en ambas
# replicas, sin marcas de histona). NF-YA = interseccion HeLa-S3 / Hep-G2 /
# K-562 (GM12878 excluida, igual que tf_contexto_endogeno.R) -> "NF-YA ∩".
# La union de NFYA en cualquier linea se excluye de la referencia (mismo
# factor) y no se grafica.
# Paneles:
#   nfya_scatter_intersect.jpg  - proporcion de promotores unidos por NF-YA ∩
#                                 por bin de actividad (= NFYA_scatter_intersect)
#   tf_efecto_vs_cobertura.jpg  - efecto (HL) vs. n de promotores unidos de
#                                 todos los TFs, ventana de cobertura +-25%
#   tf_especificidad_ccaat.jpg  - efecto con - sin CCAAT-box de todos los TFs;
#                                 control de cobertura: percentil de NF-YA ∩
#                                 tambien entre los TFs de su ventana +-25%
# Tablas: nfy_percentiles.tsv (percentil/rango en universo y ventana, para
# efecto y para con - sin CCAAT-box), nfy_top10_ccaat.tsv, nfy_top10_ventana.tsv.
#
# NO recalcula efectos: reutiliza las tablas de tf_especificidad.R (misma
# carpeta de scripts), guardadas en figures/exploratorio_tf_contexto_endogeno/:
# tf_efecto_todos_intersect.tsv y tf_motivo_especificidad_intersect.tsv.
# Correr tf_especificidad.R antes si no existen. Solo el scatter necesita
# remap_tf_hits.tsv / activity_stats_highconf.tsv / prom_df.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

out_dir <- fig_dir("nfy_especificidad")
src_dir <- "figures/exploratorio_tf_contexto_endogeno"

cov_window <- 0.25
nfy_defs <- c("NFYA_intersect", "NFYA") # ambas definiciones de NF-YA, fuera de la referencia
nfy_clr <- "#D6741F"

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in",
      device = if (ext == "pdf") cairo_pdf else NULL)
  }
}

pct_of <- function(x, ref) 100 * mean(ref <= x)
rank_of <- function(x, ref) sum(ref > x) + 1L

# =========================================================================
# Panel: scatter de NF-YA ∩ por bin de actividad
# =========================================================================

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>% add_mean_sw_bins()

nfya_lines <- remap_hits %>%
  filter(TF == "NFYA") %>%
  separate_rows(sample, sep = ",") %>%
  filter(sample != "GM12878") %>%
  group_by(sample) %>%
  summarise(fps = list(unique(seq_id)), .groups = "drop") %>%
  deframe()
stopifnot(setequal(names(nfya_lines), c("HeLa-S3", "Hep-G2", "K-562")))
nfya_intersect <- Reduce(intersect, nfya_lines)
rm(remap_hits)

data_nfya <- data %>%
  select(seq_id, rep, mean, mean_sw) %>%
  mutate(NFYA = seq_id %in% nfya_intersect)
p_scatter <- tf_scatter(data_nfya, "NFYA", "Unión a NF-YA (HeLa-S3 ∩ Hep-G2 ∩ K-562)")
save_jpg_pdf(p_scatter, "nfya_scatter_intersect", 9, 6.75)

# =========================================================================
# Efecto vs. cobertura (universo cualquier linea)
# =========================================================================
# Referencia = todos los TFs salvo NF-YA (tambien se excluye la union de
# NFYA en cualquier linea: es el mismo factor)

eff <- read_tsv(file.path(src_dir, "tf_efecto_todos_intersect.tsv"), show_col_types = FALSE) %>%
  filter(universo == "Cualquier línea") %>%
  select(rep, TF, n_unidos, estimate)
ref_eff <- eff %>% filter(!TF %in% nfy_defs)
message("Universo cualquier linea: ", n_distinct(ref_eff$TF), " TFs de referencia")

focal_eff <- eff %>%
  filter(TF == "NFYA_intersect") %>%
  pmap(function(rep, TF, n_unidos, estimate) {
    ref <- ref_eff %>% filter(rep == !!rep)
    win <- ref %>% filter(abs(n_unidos - !!n_unidos) <= cov_window * !!n_unidos)
    tibble(
      rep, TF, n_unidos, estimate,
      n_universo = nrow(ref), percentil_universo = pct_of(estimate, ref$estimate), rango_universo = rank_of(estimate, ref$estimate),
      win_min = (1 - cov_window) * n_unidos, win_max = (1 + cov_window) * n_unidos,
      n_ventana = nrow(win), percentil_ventana = pct_of(estimate, win$estimate), rango_ventana = rank_of(estimate, win$estimate)
    )
  }) %>%
  list_rbind()

# TFs de la ventana de cobertura de NF-YA ∩ (por replica), reusada en el
# panel de especificidad como control de cobertura
win_tfs <- focal_eff %>%
  select(rep, win_min, win_max) %>%
  inner_join(ref_eff, by = "rep") %>%
  filter(n_unidos >= win_min, n_unidos <= win_max) %>%
  select(rep, TF, n_unidos, estimate)

top10_win <- win_tfs %>%
  group_by(rep) %>%
  slice_max(estimate, n = 10) %>%
  mutate(posicion = row_number(), .after = rep) %>%
  ungroup()
write_tsv(top10_win, file.path(out_dir, "nfy_top10_ventana.tsv"))

ann_cov <- focal_eff %>%
  mutate(label = sprintf(
    "NF-YA ∩: p%.0f universo (%d/%d) · p%.0f ventana (%d/%d)",
    percentil_universo, rango_universo, n_universo + 1L, percentil_ventana, rango_ventana, n_ventana + 1L
  ))

p_cov <- ggplot(ref_eff, aes(n_unidos, estimate)) +
  geom_rect(data = focal_eff, aes(xmin = win_min, xmax = win_max, ymin = -Inf, ymax = Inf), fill = nfy_clr, alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_point(color = "grey60", size = 1.2, alpha = 0.7) +
  geom_point(data = focal_eff, shape = 21, color = nfy_clr, fill = nfy_clr, size = 3.5) +
  geom_label(data = ann_cov, aes(x = min(ref_eff$n_unidos), y = Inf, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5, label.size = 0, fill = alpha("white", 0.85)) +
  facet_wrap(~rep) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    x = "Promotores unidos (ReMap, escala log10)", y = "Efecto de la unión sobre\nla actividad (Hodges-Lehmann)",
    caption = "Punto naranja: NF-YA ∩ (HeLa-S3 ∩ Hep-G2 ∩ K-562). Sombreado: ventana de cobertura ±25%. Rangos: posición / n de factores comparados."
  ) +
  theme_pubr(base_size = 16) +
  theme(plot.caption = element_text(size = 10))
save_jpg_pdf(p_cov, "tf_efecto_vs_cobertura", 12, 6)

# =========================================================================
# Especificidad por CCAAT-box (universo cualquier linea)
# =========================================================================
# Control de cobertura: ademas del percentil en todo el universo, percentil
# entre los TFs de la ventana de cobertura de NF-YA ∩ (misma replica).

spec <- read_tsv(file.path(src_dir, "tf_motivo_especificidad_intersect.tsv"), show_col_types = FALSE) %>%
  filter(universo == "Cualquier línea", motivo == "CCAAT-box") # ya filtrado a >=20 unidos y >=20 no unidos por estrato
ref_spec <- spec %>% filter(!TF %in% nfy_defs)

focal_spec <- spec %>%
  filter(TF == "NFYA_intersect") %>%
  select(rep, dif_con_sin) %>%
  pmap(function(rep, dif_con_sin) {
    ref <- ref_spec %>% filter(rep == !!rep)
    win <- ref %>% semi_join(win_tfs %>% filter(rep == !!rep), by = "TF")
    tibble(
      rep, dif_con_sin,
      n_universo = nrow(ref), percentil_universo = pct_of(dif_con_sin, ref$dif_con_sin), rango_universo = rank_of(dif_con_sin, ref$dif_con_sin),
      n_ventana = nrow(win), percentil_ventana = pct_of(dif_con_sin, win$dif_con_sin), rango_ventana = rank_of(dif_con_sin, win$dif_con_sin)
    )
  }) %>%
  list_rbind()

top10_ccaat <- ref_spec %>%
  group_by(rep) %>%
  slice_max(dif_con_sin, n = 10) %>%
  mutate(posicion = row_number(), .after = rep) %>%
  ungroup() %>%
  select(rep, posicion, TF, dif_con_sin, estimate_con, estimate_sin, n_unidos_con, n_unidos_sin)
write_tsv(top10_ccaat, file.path(out_dir, "nfy_top10_ccaat.tsv"))

spec_plot <- ref_spec %>%
  left_join(win_tfs %>% transmute(rep, TF, en_ventana = TRUE), by = c("rep", "TF")) %>%
  mutate(grupo = ifelse(is.na(en_ventana), "Resto de los factores", "Ventana de cobertura ±25%"))
lab_spec <- focal_spec %>%
  mutate(label = sprintf(
    "NF-YA ∩\np%.0f universo (%d/%d)\np%.0f ventana (%d/%d)",
    percentil_universo, rango_universo, n_universo + 1L, percentil_ventana, rango_ventana, n_ventana + 1L
  ))

p_spec <- ggplot(spec_plot, aes(dif_con_sin, fill = grupo)) +
  geom_histogram(bins = 60, color = NA, position = "stack") +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
  geom_vline(data = lab_spec, aes(xintercept = dif_con_sin), color = nfy_clr, linewidth = 1) +
  geom_label(data = lab_spec, aes(x = dif_con_sin, y = Inf, label = label), inherit.aes = FALSE, color = nfy_clr, hjust = -0.05, vjust = 1.1, size = 3.5, label.size = 0, fill = alpha("white", 0.85)) +
  facet_wrap(~rep, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Resto de los factores" = "grey80", "Ventana de cobertura ±25%" = "#F2B680")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Efecto con CCAAT-box − efecto sin CCAAT-box (Hodges-Lehmann)", y = "N° de factores (ReMap)", fill = NULL) +
  theme_pubr(base_size = 16) +
  theme(legend.position = "top")
save_jpg_pdf(p_spec, "tf_especificidad_ccaat", 9, 7.5)

# Figura completa: scatter (A) y efecto vs. cobertura (B) a la izquierda, uno
# sobre otro (los dos con Rep 1 | Rep 2 lado a lado); con - sin CCAAT-box (C)
# a la derecha, a lo alto
figura <- wrap_elements(full = p_scatter) + wrap_elements(full = p_cov) + wrap_elements(full = p_spec) +
  plot_layout(design = "AC
BC", widths = c(12, 8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 26, face = "bold"))
save_jpg_pdf(figura, "nfy_especificidad", 20, 13)

# =========================================================================
# Resumen de percentiles (NF-YA ∩)
# =========================================================================

pct_summary <- bind_rows(
  focal_eff %>% transmute(analisis = "efecto", rep, valor = estimate, n_unidos, n_universo, percentil_universo, rango_universo, n_ventana, percentil_ventana, rango_ventana),
  focal_spec %>% transmute(analisis = "con - sin CCAAT-box", rep, valor = dif_con_sin, n_universo, percentil_universo, rango_universo, n_ventana, percentil_ventana, rango_ventana)
)
print(pct_summary, width = Inf)
write_tsv(pct_summary, file.path(out_dir, "nfy_percentiles.tsv"))

message("Figuras y tablas guardadas en ", out_dir)
