# shape_especificidad_tisular (ver docs/mapping_figuras.csv para el
# numero de figura vigente). Grilla 2x2, en el orden en que Metodos cita
# los paneles (primero especificidad, despues forma):
#
# Panel A: heatmap de especificidad tisular. Filas = promotores EPD de la
#          library con sample_specificity_gini no NA, ordenados por Gini
#          (menor a mayor). Columnas = las mismas unidades sobre las que
#          se calcula el Gini en build_prom_features.R: la columna `sample`
#          de sample_CAGE_activity.tsv, que son 110 terminos de ontologia
#          de FANTOM5 (69 de celulas primarias + 41 de tejidos, cada uno
#          sumando varias muestras individuales). Color = fraccion de la
#          actividad total del promotor en cada termino (TPM normalizado
#          por fila). Columnas ordenadas por clustering jerarquico
#          (distancia 1 - correlacion, enlace promedio); barra superior
#          indica tejido / celula primaria. Cortes de terciles de Gini
#          marcados sobre el eje de filas.
# Panel B: densidad del indice de Gini con los cortes de terciles.
# Panel C: dos promotores de ejemplo de forma (angosto: TATA sin isla CpG;
#          ancho: isla CpG sin TATA), lecturas de CAGE por posicion
#          respecto del TSS de EPD, sumando los 109 archivos de
#          Library_filt_CAGE (tejidos + celulas primarias, la misma senal
#          que se uso para el shape en transcriptional_library/Analysis/
#          scripts/shape_merged.qmd). Sombreado = intervalo entre los
#          cuantiles 0.05 y 0.95 del cluster de CAGE (qLow/qUp con que se
#          calculo interquantile_width en shape_merged.qmd).
# Panel D: densidad del ancho interquantile con los cortes de terciles.
#
# Terciles: quantile(x, c(1/3, 2/3)) sobre todos los promotores de
# prom_df, mismo criterio de cut_number(x, 3) usado en R/09_coocurrencia y
# R/15_summary_features.
#
# Requiere data/processed/prom_df.tsv, data/raw/library.bed y las
# excepciones de R/00_prom_features/analysis_tables_exceptions.R
# (sample_CAGE_activity.tsv, shape_merged.tsv). El panel C lee ~2.9GB de
# Library_filt_CAGE (heavy_data_paths.R) la primera vez; los intermedios
# se cachean en data/processed/gini_heatmap_matrix_v1.tsv y
# data/processed/cage_shape_ejemplos_v1.tsv (borrarlos para recalcular).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/00_prom_features/analysis_tables_exceptions.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "shape_especificidad_tisular"
out_dir <- fig_dir(slug)
base_size <- 16
clr_esp <- "#AD343E"

save_panel <- function(p, name, width = 9, height = 6.75) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = width, height = height, units = "in")
  }
}

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

terciles <- function(x) quantile(x, c(1 / 3, 2 / 3), na.rm = TRUE)

density_terciles <- function(data, var, xlab, fill) {
  ggplot(data, aes(.data[[var]])) +
    geom_density(fill = fill, alpha = 0.6, col = fill, linewidth = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = terciles(data[[var]]), linetype = "dashed", linewidth = 0.5, col = "#333333") +
    labs(x = xlab, y = "Densidad") +
    theme_pubclean(base_size = base_size)
}

# Misma funcion que build_prom_features.R, para verificar el Gini.
gini <- function(x, weights = rep(1, length = length(x))) {
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox] / sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu / nu[n]
  sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

# --- Panel A: heatmap de especificidad tisular -----------------------------

heat_cache <- "data/processed/gini_heatmap_matrix_v1.tsv"
prom_gini <- prom_df %>% filter(!is.na(sample_specificity_gini)) %>% select(name, sample_specificity_gini)

if (file.exists(heat_cache)) {
  heat_wide <- read_tsv(heat_cache, show_col_types = FALSE)
} else {
  # Mismo TPM que build_prom_features.R (counts * 1e6 / libsize por termino).
  heat_wide <- read_tsv(path_sample_cage_activity, show_col_types = FALSE) %>%
    filter(name %in% prom_gini$name) %>%
    mutate(tpm = counts * 1e6 / libsize) %>%
    group_by(name) %>%
    mutate(gini_recalc = gini(tpm), frac = tpm / sum(tpm)) %>%
    ungroup() %>%
    select(name, gini_recalc, sample, frac) %>%
    pivot_wider(names_from = sample, values_from = frac)
  write_tsv(heat_wide, heat_cache)
}

heat_wide <- heat_wide %>%
  inner_join(prom_gini, by = "name") %>%
  arrange(sample_specificity_gini)
stopifnot(nrow(heat_wide) == nrow(prom_gini))
gini_diff <- max(abs(heat_wide$gini_recalc - heat_wide$sample_specificity_gini))
message("Heatmap: ", nrow(heat_wide), " promotores x ", ncol(heat_wide) - 3, " terminos de ontologia; ",
        "max |Gini recalculado - prom_df| = ", signif(gini_diff, 3))
stopifnot(gini_diff < 1e-6)

frac_mat <- heat_wide %>% select(-name, -gini_recalc, -sample_specificity_gini) %>% as.matrix()

ontologia <- bind_rows(
  readxl::read_xls("data/external/FANTOM5/primary_cell_ontology_FANTOM5.xlsx") %>% mutate(origen = "Células primarias"),
  readxl::read_xls(path_fantom5_tissue_ontology) %>% mutate(origen = "Tejidos")
) %>%
  distinct(sample = `Facet ontology term`, origen)
stopifnot(all(colnames(frac_mat) %in% ontologia$sample), !anyDuplicated(ontologia$sample))

# Tejidos y celulas primarias en bloques separados; clustering jerarquico
# dentro de cada bloque.
orden_cols <- ontologia %>%
  filter(sample %in% colnames(frac_mat)) %>%
  group_by(origen) %>%
  group_modify(~ {
    m <- frac_mat[, .x$sample, drop = FALSE]
    tibble(sample = colnames(m)[hclust(as.dist(1 - cor(m)), method = "average")$order])
  }) %>%
  mutate(col = row_number()) %>%
  ungroup() %>%
  mutate(origen = factor(origen, levels = c("Tejidos", "Células primarias")))

heat_long <- heat_wide %>%
  mutate(fila = row_number()) %>%
  select(-name, -gini_recalc, -sample_specificity_gini) %>%
  pivot_longer(-fila, names_to = "sample", values_to = "frac") %>%
  inner_join(orden_cols, by = "sample")

# Cortes de terciles (sobre el mismo vector que el panel B) pasados a
# posicion de fila.
cortes_gini <- terciles(prom_df$sample_specificity_gini)
filas_corte <- map_int(cortes_gini, ~ sum(heat_wide$sample_specificity_gini <= .x)) + 0.5

frac_max <- 0.25
panel_a <- ggplot(heat_long, aes(col, fila, fill = frac)) +
  geom_raster() +
  geom_hline(yintercept = filas_corte, linetype = "dashed", linewidth = 0.5, col = "#333333") +
  scale_fill_gradientn(
    colours = c("white", "#E8B9BD", clr_esp, "#4A1419"), limits = c(0, frac_max),
    trans = "sqrt", oob = scales::squish, breaks = c(0, 0.05, 0.1, 0.25),
    labels = c("0", "0.05", "0.1", "≥0.25"), name = "Fracción de\nla actividad"
  ) +
  scale_y_continuous(breaks = filas_corte, labels = sprintf("Gini %.2f", cortes_gini), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(~origen, scales = "free_x", space = "free_x") +
  labs(x = "Términos de ontología de FANTOM5", y = "Promotores (ordenados por índice de Gini)") +
  theme_pubclean(base_size = base_size) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank(), panel.spacing.x = unit(0.4, "cm"),
    strip.background = element_blank(), strip.text = element_text(size = base_size),
    legend.position = "right", legend.key.height = unit(1, "cm")
  )

save_panel(panel_a, "panel_a_heatmap_especificidad")

# --- Panel B: densidad del indice de Gini ----------------------------------

panel_b <- density_terciles(prom_df, "sample_specificity_gini", "Especificidad tisular (índice de Gini)", clr_esp)
save_panel(panel_b, "panel_b_gini")

# --- Panel C: ejemplos de forma (CAGE por posicion) ------------------------

cortes_ancho <- terciles(prom_df$interquantile_width)
med_angosto <- median(prom_df$interquantile_width[prom_df$interquantile_width <= cortes_ancho[1]], na.rm = TRUE)
med_ancho <- median(prom_df$interquantile_width[prom_df$interquantile_width > cortes_ancho[2]], na.rm = TRUE)

shape <- read_tsv(path_shape_merged, show_col_types = FALSE) %>%
  select(name, cl_start = start, cl_end = end, cl_score = score.x)

# Angosto: candidatos con TATA sin CGI y ancho interquantile a +-1 pb de la
# mediana del tercil angosto; los 8 de mayor senal de CAGE en el cluster
# (se elige uno mas abajo).
# Ancho: fijado a mano por el autor (HDAC2_1: CGI sin TATA, 60 pb, senal
# repartida a lo largo del cluster sin un pico dominante) despues de
# revisar los candidatos con CGI sin TATA cerca de la mediana del tercil
# ancho ordenados por fraccion de lecturas en la posicion pico.
ejemplo_ancho <- "HDAC2_1"

candidatos <- bind_rows(
  prom_df %>% filter(TATA_EPD, !CGI_frag, abs(interquantile_width - med_angosto) <= 1) %>% mutate(clase = "angosto"),
  prom_df %>% filter(name == ejemplo_ancho) %>% mutate(clase = "ancho")
) %>%
  select(seq_id, name, clase, interquantile_width, TATA_EPD, CGI_frag) %>%
  inner_join(shape, by = "name") %>%
  group_by(clase) %>%
  slice_max(cl_score, n = 8, with_ties = FALSE) %>%
  ungroup()
fijo <- candidatos %>% filter(clase == "ancho")
stopifnot(nrow(fijo) == 1, fijo$CGI_frag, !fijo$TATA_EPD, fijo$interquantile_width > cortes_ancho[2])

lib <- read_tsv("data/raw/library.bed", col_names = c("chr", "s0", "e", "seq_id", "sc", "strand"),
                col_select = 1:6, show_col_types = FALSE)

cage_cache <- "data/processed/cage_shape_ejemplos_v3.tsv"
if (file.exists(cage_cache)) {
  cage_ej <- read_tsv(cage_cache, show_col_types = FALSE)
} else {
  # Senal por posicion de los 109 terminos (tejidos + celulas primarias),
  # sumada por posicion, como en shape_merged.qmd. Solo hebra del promotor.
  bed_files <- list.files(path_fantom5_filt_cage, pattern = "\\.bed$", full.names = TRUE)
  cage_ej <- map(bed_files, function(f) {
    data.table::fread(f, header = FALSE, col.names = c("chr", "s0", "e", "seq_id", "score", "strand"))[seq_id %in% candidatos$seq_id]
  }) %>%
    data.table::rbindlist() %>%
    as_tibble() %>%
    semi_join(lib, by = c("seq_id", "strand")) %>%
    mutate(pos = s0 + 1) %>%
    group_by(seq_id, pos) %>%
    summarise(score = sum(score), .groups = "drop")
  write_tsv(cage_ej, cage_cache)
}
stopifnot(all(candidatos$seq_id %in% cage_ej$seq_id))

# Cuantil de CAGEr: primera posicion (orden genomico) donde la suma
# acumulada alcanza la fraccion q del total del cluster.
qpos <- function(pos, sc, q) {
  o <- order(pos)
  pos[o][which(cumsum(sc[o]) / sum(sc) >= q)[1]]
}

# Posicion relativa al TSS de EPD (posicion 237 del fragmento de 252 pb,
# geometria -236/+15), respetando la hebra.
to_rel <- function(pos, s0, e, strand) ifelse(strand == "+", pos - s0, e - pos + 1) - 237

verif <- cage_ej %>%
  inner_join(candidatos, by = "seq_id") %>%
  filter(pos >= cl_start, pos <= cl_end) %>%
  group_by(seq_id) %>%
  summarise(q05 = qpos(pos, score, 0.05), q95 = qpos(pos, score, 0.95), lecturas = sum(score),
            frac_pico = max(score) / sum(score)) %>%
  mutate(ancho_recalc = q95 - q05 + 1) %>%
  inner_join(candidatos, by = "seq_id") %>%
  left_join(lib %>% select(seq_id, s0, e, strand), by = "seq_id") %>%
  mutate(rel_a = to_rel(q05, s0, e, strand), rel_b = to_rel(q95, s0, e, strand))

# Se exige que el ancho recalculado coincida (+-1 pb) con prom_df y que el
# intervalo 5-95 quede al menos 5 pb dentro del fragmento (la senal esta
# recortada a -236/+15, un cluster pegado al borde se ve truncado).
# Angosto: el mas cercano a la mediana de su tercil, desempate por lecturas.
# El ancho fijado tiene que pasar los mismos controles (si no, falla).
validos <- verif %>%
  filter(abs(ancho_recalc - interquantile_width) <= 1, pmin(rel_a, rel_b) >= -231, pmax(rel_a, rel_b) <= 10)
stopifnot(ejemplo_ancho %in% validos$name)
ejemplos <- bind_rows(
  validos %>% filter(clase == "angosto") %>% arrange(abs(interquantile_width - med_angosto), desc(lecturas)) %>% slice(1),
  validos %>% filter(clase == "ancho")
) %>%
  mutate(
    ini = pmin(rel_a, rel_b) - 0.5, fin = pmax(rel_a, rel_b) + 0.5,
    gen = str_remove(name, "_\\d+$"),
    etiqueta = sprintf("%s (%s): ancho interquantile %d pb",
                       if_else(clase == "angosto", "Angosto", "Ancho"), gen, as.integer(interquantile_width)),
    etiqueta = fct_reorder(etiqueta, clase == "ancho")
  )
stopifnot(nrow(ejemplos) == 2)
message("Ejemplos de forma:")
print(ejemplos %>% select(clase, seq_id, name, interquantile_width, ancho_recalc, lecturas) %>% as.data.frame())

senal_ej <- cage_ej %>%
  inner_join(ejemplos %>% select(seq_id, s0, e, strand, etiqueta), by = "seq_id") %>%
  mutate(relpos = to_rel(pos, s0, e, strand)) %>%
  filter(relpos >= -236, relpos <= 15)

# Limite izquierdo: 20 pb antes del intervalo sombreado mas extenso (se
# omite senal suelta lejos del cluster).
x_min <- max(-236, floor((min(ejemplos$ini) - 20) / 10) * 10)

panel_c <- ggplot(senal_ej, aes(relpos, score / 1000)) +
  geom_rect(data = ejemplos, aes(xmin = ini, xmax = fin, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = thesis_clr_bg) +
  geom_col(width = 1, fill = thesis_clr) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, col = "#333333") +
  facet_wrap(~etiqueta, ncol = 1, scales = "free_y") +
  scale_x_continuous(limits = c(x_min, 15.5), breaks = c(seq(-200, 0, by = 20), 15)) +
  labs(x = "Posición relativa al TSS de EPD (pb)", y = "Lecturas de CAGE (miles)") +
  theme_pubclean(base_size = base_size) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, size = base_size))

save_panel(panel_c, "panel_c_ejemplos_forma")

# --- Panel D: densidad del ancho interquantile -----------------------------

panel_d <- density_terciles(prom_df, "interquantile_width", "Ancho interquantile (pb)", thesis_clr)
save_panel(panel_d, "panel_d_ancho")

# --- Figura compuesta ------------------------------------------------------

compuesta <- wrap_elements(full = panel_a) + panel_b + panel_c + panel_d +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

save_panel(compuesta, slug, width = 18, height = 13.5)

message("Figuras guardadas en ", out_dir)
