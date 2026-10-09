# orientacion_promotor (ver docs/mapping_figuras.csv para el numero vigente) - Orientación transcripcional
# del promotor (unidireccional vs. bidireccional), a partir de PRO-cap
# (HEK293, ENCODE). Ver R/00_prom_features/build_promoter_orientation.R
# para el cómputo de la Orientation Index (OI) y la clasificación
# (columna `orientacion` de prom_df.tsv).
#
# Un solo panel con 4 gráficos:
#   1. Histograma de señal en sentido (Rp), para justificar el umbral
#      MIN_SENSE=10 usado para descartar promotores sin señal suficiente.
#   2. Histograma de la Orientation Index (OI), con el umbral OI>=0.9
#      usado para clasificar unidireccional vs. bidireccional.
#   3. Composición: proporción de promotores por clase de orientación.
#   4. Metaplot de señal de PRO-cap (sentido/antisentido) alrededor del
#      TSS, separado en facets por clase - la validación de que la
#      clasificación captura señal real.
#
# A diferencia de promoter_orientation.R (transcriptional_library), NO
# incluye la validación contra los BED pre-clasificados de ENCODE (paso
# 5 de ese script) ni la comparación entre réplicas (sección 8 de ese
# script, tampoco corresponde acá - eso era efecto sobre actividad/ruido,
# no la clasificación en sí).
#
# Requiere data/raw/library.bed, data/processed/prom_df.tsv,
# data/processed/promoter_orientation.tsv (ver build_promoter_
# orientation.R) y los bigWig de PRO-cap (heavy_data_paths.R). Run from
# the TesisDoc repo root.

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(patchwork)
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "orientacion_promotor"
out_dir <- fig_dir(slug)

META_WIN <- 500L
BIN_SIZE <- 10L
MIN_SENSE <- 10L

# Mismos colores que shape_especificidad_tisular (shape_especificidad_tisular.R): azul =
# thesis_clr ("#358AAA"), rojo = "#AD343E".
cls_colors <- c(bidireccional = "#AD343E", unidireccional = "#358AAA")
color_umbral_senal <- "#7FB800" # verde oliva, elegido para combinar con cls_colors

# --- Datos ya calculados -----------------------------------------------

orientation_df <- read_tsv("data/processed/promoter_orientation.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)

# --- Panel 1: histograma de señal en sentido (umbral MIN_SENSE) --------

p_hist <- orientation_df %>%
  mutate(bajo_umbral = Rp_total < MIN_SENSE) %>%
  ggplot(aes(x = log10(Rp_total + 1), fill = bajo_umbral)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.1) +
  geom_vline(xintercept = log10(MIN_SENSE + 1), linetype = "dashed", color = "tomato") +
  annotate("text", x = log10(MIN_SENSE + 1), y = Inf, vjust = 1.5, hjust = -0.1,
           label = paste0("umbral = ", MIN_SENSE), size = 3.5, color = "tomato") +
  scale_fill_manual(values = c("TRUE" = "#BBBBBB", "FALSE" = color_umbral_senal), guide = "none") +
  labs(x = "log10(señal en sentido + 1)", y = "Promotores",
       title = "Umbral de señal mínima") +
  theme_bw(base_size = 13)

# --- Panel 2: histograma de la Orientation Index (umbral OI_THRESH) ----

OI_THRESH <- 0.9

p_oi <- orientation_df %>%
  filter(!is.na(orientacion)) %>%
  ggplot(aes(x = OI, fill = orientacion)) +
  geom_histogram(binwidth = 0.05, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = OI_THRESH, linetype = "dashed", color = "tomato") +
  annotate("text", x = OI_THRESH - 0.01, y = Inf, vjust = 1.5, hjust = 1,
           label = paste0("umbral = ", OI_THRESH), size = 3.5, color = "tomato") +
  scale_fill_manual(values = cls_colors) +
  labs(x = "Orientation Index (OI)", y = "Promotores", fill = NULL,
       title = "Umbral de la Orientation Index") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", legend.text = element_text(size = 13))

# --- Panel 3: composición por clase -------------------------------------

p_comp <- orientation_df %>%
  mutate(orientacion = replace_na(orientacion, "no clasificado")) %>%
  count(orientacion) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = "", y = pct, fill = orientacion)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%s\n%.1f%%\nn=%d", orientacion, pct, n)),
            position = position_stack(vjust = 0.5), size = 3.2, lineheight = 0.9) +
  scale_fill_manual(values = c(cls_colors, "no clasificado" = "#BBBBBB")) +
  labs(x = NULL, y = "% de promotores", title = "Composición") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

# --- Panel 4: metaplot de señal de PRO-cap por clase (facets) -----------

lib <- rtracklayer::import.bed("data/raw/library.bed")
lib_fp <- lib[startsWith(lib$name, "FP")]
plus_idx <- as.character(strand(lib_fp)) == "+"
tss_pos <- ifelse(plus_idx, start(lib_fp) + 236L, end(lib_fp) - 236L)

n_bins <- (2L * META_WIN) %/% BIN_SIZE
bin_mids <- seq(-META_WIN + BIN_SIZE / 2, META_WIN - BIN_SIZE / 2, by = BIN_SIZE)

meta_gr <- GRanges(
  seqnames = seqnames(lib_fp),
  ranges = IRanges(start = pmax(1L, tss_pos - META_WIN), end = tss_pos + META_WIN),
  strand = strand(lib_fp)
)

build_meta_matrix <- function(bw_path) {
  hits <- tryCatch(rtracklayer::import(BigWigFile(bw_path), which = meta_gr), error = function(e) GRanges())
  mat <- matrix(0, nrow = length(lib_fp), ncol = n_bins)
  if (length(hits) == 0L) return(mat)
  hits$score <- abs(score(hits))
  ov <- findOverlaps(meta_gr, hits)
  qi <- queryHits(ov)
  si <- subjectHits(ov)
  h_mid <- (start(hits)[si] + end(hits)[si]) %/% 2L
  rel_pos <- ifelse(plus_idx[qi], h_mid - tss_pos[qi], tss_pos[qi] - h_mid)
  bin_i <- (rel_pos + META_WIN) %/% BIN_SIZE + 1L
  valid <- bin_i >= 1L & bin_i <= n_bins
  if (!any(valid)) return(mat)
  agg <- tapply(score(hits)[si[valid]], paste(qi[valid], bin_i[valid]), sum)
  parts <- strsplit(names(agg), " ")
  row_i <- as.integer(vapply(parts, `[`, character(1), 1))
  col_i <- as.integer(vapply(parts, `[`, character(1), 2))
  mat[(col_i - 1L) * nrow(mat) + row_i] <- as.numeric(agg)
  mat
}

message("Construyendo matrices de metaplot (4 bigWig)...")
mat_plus_r1 <- build_meta_matrix(path_procap_plus_r1)
mat_minus_r1 <- build_meta_matrix(path_procap_minus_r1)
mat_plus_r2 <- build_meta_matrix(path_procap_plus_r2)
mat_minus_r2 <- build_meta_matrix(path_procap_minus_r2)

mat_sense_r1 <- mat_plus_r1; mat_sense_r1[!plus_idx, ] <- mat_minus_r1[!plus_idx, ]
mat_anti_r1 <- mat_minus_r1; mat_anti_r1[!plus_idx, ] <- mat_plus_r1[!plus_idx, ]
mat_sense_r2 <- mat_plus_r2; mat_sense_r2[!plus_idx, ] <- mat_minus_r2[!plus_idx, ]
mat_anti_r2 <- mat_minus_r2; mat_anti_r2[!plus_idx, ] <- mat_plus_r2[!plus_idx, ]

mat_sense_avg <- (mat_sense_r1 + mat_sense_r2) / 2
mat_anti_avg <- (mat_anti_r1 + mat_anti_r2) / 2

orient_by_seq <- orientation_df$orientacion[match(lib_fp$name, orientation_df$seq_id)]
is_bidi <- orient_by_seq == "bidireccional" & !is.na(orient_by_seq)
is_uni <- orient_by_seq == "unidireccional" & !is.na(orient_by_seq)

meta_df <- tibble(
  position = bin_mids,
  sentido_bidi = colMeans(mat_sense_avg[is_bidi, , drop = FALSE]),
  sentido_uni = colMeans(mat_sense_avg[is_uni, , drop = FALSE]),
  antisentido_bidi = colMeans(mat_anti_avg[is_bidi, , drop = FALSE]),
  antisentido_uni = colMeans(mat_anti_avg[is_uni, , drop = FALSE])
)

p_meta <- meta_df %>%
  pivot_longer(-position, names_to = c("hebra", "clase"), names_sep = "_", values_to = "señal") %>%
  mutate(
    clase = factor(clase, levels = c("bidi", "uni"), labels = c("Bidireccional", "Unidireccional")),
    hebra = factor(hebra, levels = c("sentido", "antisentido"), labels = c("Sentido", "Antisentido"))
  ) %>%
  ggplot(aes(x = position, y = señal, color = clase, linetype = hebra)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  facet_wrap(~clase, scales = "free_y") +
  scale_color_manual(values = c(Bidireccional = cls_colors[["bidireccional"]], Unidireccional = cls_colors[["unidireccional"]]), guide = "none") +
  scale_linetype_manual(values = c(Sentido = "solid", Antisentido = "dashed")) +
  labs(x = "Posición relativa al TSS (pb)", y = "Señal media de PRO-cap",
       title = "Metaplot de orientación", color = NULL, linetype = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", legend.text = element_text(size = 13))

# --- Panel combinado -----------------------------------------------------

p_combined <- (p_hist | p_oi | p_comp) / p_meta +
  plot_annotation(title = "Orientación transcripcional del promotor (PRO-cap, HEK293)")

ggsave(file.path(out_dir, "orientacion_promotor.jpg"), p_combined, width = 14, height = 9, units = "in")
message("Figura guardada en ", out_dir)
