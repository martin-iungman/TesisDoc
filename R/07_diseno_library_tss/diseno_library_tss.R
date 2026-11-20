# diseno_library_tss (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Panel A: esquema de estructura de secuencias con adaptadores - manual,
#          NO se genera acá (dibujado con formas nativas de PowerPoint en
#          docs/Fig MyM.pptx, sin datos de por medio). Exportar a mano
#          desde el pptx y guardar en fig_dir(slug) como
#          panel_a_esquema_adaptadores.png; si existe, se pega en la
#          figura compuesta.
# Panel B: cantidad de secuencias diseñadas por tipo (promotor EPD,
#          enhancer, promotor asociado a cáncer WT + mutante), desde la
#          tabla de diseño completa data/raw/seq_data.tsv (incluye las
#          mutantes, que no tienen fila en library.bed), sin duplicados.
# Panel C: heatmap de frecuencia de uso de TSS por CAGE (cultivos primarios,
#          FANTOM5), recortado a la region proximal (-50pb a +15pb del TSS)
#
# Requiere: data/raw/seq_data.tsv y data/raw/library.bed.
# Panel C ademas lee ~2.9GB de datos de FANTOM5 directo de
# transcriptional_library (ver R/00_prom_features/heavy_data_paths.R) -
# tarda ~1-2 minutos; el resultado se cachea en
# data/processed/cage_tss_relpos_v1.tsv y las corridas siguientes lo
# reusan (borrarlo para recalcular). Run from the TesisDoc repo root.

library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(tidyverse)
library(ggpubr)
library(patchwork)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "diseno_library_tss"
out_dir <- fig_dir(slug)
base_size <- 20

save_panel <- function(p, name, width = 9, height = 6.75) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = width, height = height, units = "in")
  }
}

# --- Panel B: cantidad de secuencias por tipo -----------------------------

type_labels <- tibble(
  type = c("promoter", "enhancer", "cancer_wt", "cancer_mut"),
  # WT y mutantes de cáncer van juntos en una sola barra
  tipo_secuencia = c("Promotor EPD", "Enhancer", "Promotor asociado\na cáncer", "Promotor asociado\na cáncer")
)

# Sin duplicados: dupl_id marca secuencias idénticas con otro ID. Se descarta
# la secuencia a la que apunta dupl_id (53 pares EPD-EPD y 36 pares EPD-WT de
# cáncer, donde siempre apunta al WT), así cada secuencia cuenta una sola vez
# y ante la duda queda asignada como promotor EPD.
seq_data <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE)
conteos_tipo <- seq_data %>%
  filter(!seq_id %in% seq_data$dupl_id) %>%
  left_join(type_labels, by = "type") %>%
  count(tipo_secuencia) %>%
  mutate(tipo_secuencia = factor(tipo_secuencia, levels = rev(unique(type_labels$tipo_secuencia))))
stopifnot(!anyNA(conteos_tipo$tipo_secuencia))

panel_b <- conteos_tipo %>%
  ggplot(aes(n, tipo_secuencia, label = n)) +
  geom_col(fill = thesis_clr) +
  geom_label(hjust = 0, nudge_x = max(conteos_tipo$n) * 0.02, fill = thesis_clr_label, size = 6) +
  scale_x_continuous(breaks = seq(0, 20000, by = 5000), expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Cantidad de secuencias", y = NULL) +
  theme_pubclean(base_size = base_size, flip = TRUE) +
  theme(axis.text = element_text(colour = "black"))

save_panel(panel_b, "panel_b_tipo_secuencia")

# --- Panel C: heatmap de uso de TSS por CAGE (FANTOM5 primary cells) -----

relpos_cache <- "data/processed/cage_tss_relpos_v1.tsv"

if (file.exists(relpos_cache)) {
  relpos_df <- read_tsv(relpos_cache, show_col_types = FALSE)
} else {
  lib <- rtracklayer::import.bed("data/raw/library.bed")
  values(lib) <- tibble(seq_id = lib$name)

  # El archivo tiene extension .xlsx pero es formato .xls legado - read_xls
  # es intencional, read_xlsx falla.
  cell_ont <- readxl::read_xls("data/external/FANTOM5/primary_cell_ontology_FANTOM5.xlsx") %>%
    rename(sample_id = `Sample ID`)

  # Solo para listar que muestras existen (nunca se lee su contenido - son
  # ~6.7GB de CTSS crudo, la señal real sale de Library_filt_CAGE mas abajo).
  raw_sample_files <- list.files(path_fantom5_primary_cell_raw, full.names = TRUE, pattern = ".nobarcode.ctss.bed.gz$")
  cell_ont <- tibble(filename = raw_sample_files, sample_id = str_extract(raw_sample_files, "CNhs.{5}")) %>%
    inner_join(cell_ont, by = "sample_id")

  # Señal de CAGE ya pre-filtrada a las posiciones de la library, un archivo
  # por termino de ontologia ("Facet ontology term").
  filt_files <- file.path(path_fantom5_filt_cage, paste0(unique(cell_ont$`Facet ontology term`), ".bed"))
  filt_files <- filt_files[filt_files %in% list.files(path_fantom5_filt_cage, full.names = TRUE)]

  beds <- map(filt_files, import.bed)
  group_bed <- do.call("c", beds)
  beds_sum <- group_bed %>%
    as_tibble() %>%
    group_by(seqnames, start, strand, name) %>%
    summarise(sum_score = sum(score, na.rm = TRUE), .groups = "drop")
  beds_sum_gr <- GRanges(
    seqnames = beds_sum$seqnames, IRanges(start = beds_sum$start, width = 1),
    strand = beds_sum$strand, name = beds_sum$name, score = beds_sum$sum_score
  )

  # Join dirigido (mismo strand) contra la library para saber a que
  # promotor/enhancer pertenece cada posicion con señal de CAGE.
  joined <- beds_sum_gr %>% plyranges::join_overlap_left_directed(lib)

  lib_df <- as_tibble(lib) %>% select(seq_id, lib_start = start, lib_width = width)
  cage_lib <- as_tibble(joined) %>%
    filter(!is.na(seq_id)) %>%
    left_join(lib_df, by = "seq_id") %>%
    mutate(
      # posicion 1-based dentro de la secuencia, ascendente con la
      # coordenada genomica (independiente de strand)
      relpos_genomic = start - lib_start + 1,
      # para strand "-", invertir para que relpos=1 sea siempre el extremo
      # 5' en la direccion de transcripcion, igual que en strand "+"
      relpos = ifelse(strand == "-", lib_width - relpos_genomic + 1, relpos_genomic)
    )

  relpos_df <- cage_lib %>%
    group_by(seq_id) %>%
    mutate(relscore = score / sum(score), sumscore = sum(score)) %>%
    ungroup()

  relpos_df <- relpos_df %>% select(seq_id, relpos, score, relscore, sumscore)
  write_tsv(relpos_df, relpos_cache)
}

# El offset de 237 alinea el TSS anotado (posicion ~237 de 252) en el 0.
# Recortado a la region proximal (-50pb/+15pb) como en la referencia.
panel_c <- relpos_df %>%
  filter(sumscore > 10, str_detect(seq_id, "^FP")) %>%
  arrange(desc(relscore)) %>%
  mutate(seq_id = fct_inorder(seq_id), relpos = relpos - 237) %>%
  filter(relpos > -51) %>%
  ggplot(aes(relpos, seq_id, fill = relscore * 100)) +
  geom_tile(na.rm = TRUE) +
  # magma invertido pero arrancando en blanco, para que el fondo (~0%)
  # no quede amarillo y el TSS resalte
  scale_fill_gradientn(colours = c("white", rev(viridisLite::magma(8, end = 0.9))[-1])) +
  labs(x = "Posición relativa al TSS", y = "Promotor", fill = "Actividad\nrelativa (%)") +
  scale_x_continuous(breaks = c(15, seq(0, -50, by = -15))) +
  theme_pubr(base_size = base_size) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.key.width = unit(1.2, "cm"))

save_panel(panel_c, "panel_c_cage_tss_heatmap")

# --- Figura compuesta ---------------------------------------------------------
# Si el esquema del panel A ya fue exportado, va arriba a todo el ancho;
# si no, la compuesta sale con B y C solos (tags B y C igual).

panel_a_png <- file.path(out_dir, "panel_a_esquema_adaptadores.png")
fila_bc <- (panel_b | panel_c) + plot_layout(widths = c(1.2, 1.5))

if (file.exists(panel_a_png)) {
  panel_a <- wrap_elements(full = grid::rasterGrob(png::readPNG(panel_a_png)))
  compuesta <- panel_a / fila_bc + plot_layout(heights = c(0.6, 1))
  alto <- 11.5
  tags <- list(c("A", "B", "C"))
} else {
  compuesta <- fila_bc
  alto <- 6.75
  tags <- list(c("B", "C"))
}

compuesta <- compuesta +
  plot_annotation(tag_levels = tags) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

save_panel(compuesta, "diseno_library_tss", width = 16, height = alto)

message("Paneles y figura guardados en ", out_dir,
        if (!file.exists(panel_a_png)) ". Panel A (esquema manual) sigue pendiente." else "")
