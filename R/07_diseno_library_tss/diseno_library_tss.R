# diseno_library_tss (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Panel A: esquema de estructura de secuencias con adaptadores - manual,
#          NO se genera acá (dibujado con formas nativas de PowerPoint en
#          docs/Fig MyM.pptx, sin datos de por medio). Exportar a mano
#          desde el pptx y guardar en fig_dir(slug) como
#          panel_a_esquema_adaptadores.png.
# Panel B: cantidad de secuencias por tipo (promotor/enhancer/cancer)
# Panel C: heatmap de frecuencia de uso de TSS por CAGE (cultivos primarios,
#          FANTOM5), recortado a la region proximal (-50pb a +15pb del TSS)
#
# Requiere: data/processed/prom_df.tsv (generado por
# R/00_prom_features/build_prom_features.R) y data/raw/library.bed.
# Panel C ademas lee ~2.9GB de datos de FANTOM5 directo de
# transcriptional_library (ver R/00_prom_features/heavy_data_paths.R) -
# tarda ~1-2 minutos. Run from the TesisDoc repo root.

library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(tidyverse)
library(ggpubr)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "diseno_library_tss"
out_dir <- fig_dir(slug)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)

# --- Panel B: cantidad de secuencias por tipo -----------------------------

type_labels <- tibble(
  type = c("promoter", "enhancer", "cancer_wt"),
  tipo_secuencia = c("Promotor", "Enhancer", "Cancer related\npromoter")
)

panel_b <- prom_df %>%
  count(type) %>%
  left_join(type_labels, by = "type") %>%
  mutate(tipo_secuencia = fct_reorder(tipo_secuencia, n)) %>%
  ggplot(aes(tipo_secuencia, n, label = n)) +
  geom_col(fill = thesis_clr) +
  geom_label(nudge_y = max(prom_df %>% count(type) %>% pull(n)) * 0.06, fill = thesis_clr_label, size = 5) +
  labs(y = "Frecuencia", x = "Tipo de secuencia") +
  theme_pubclean(base_size = 16)

ggsave(file.path(out_dir, "panel_b_tipo_secuencia.jpg"), panel_b, width = 9, height = 6.75, units = "in")

# --- Panel C: heatmap de uso de TSS por CAGE (FANTOM5 primary cells) -----

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

# El offset de 237 alinea el TSS anotado (posicion ~237 de 252) en el 0.
# Recortado a la region proximal (-50pb/+15pb) como en la referencia.
panel_c <- relpos_df %>%
  filter(sumscore > 10, str_detect(seq_id, "^FP")) %>%
  arrange(desc(relscore)) %>%
  mutate(seq_id = fct_inorder(seq_id), relpos = relpos - 237) %>%
  filter(relpos > -51) %>%
  ggplot(aes(relpos, seq_id, fill = relscore * 100)) +
  geom_tile(na.rm = TRUE) +
  scale_fill_gradient(high = "black", low = "white") +
  labs(x = "Posición relativa", y = "Promotor", fill = "Actividad\nrelativa (%)") +
  scale_x_continuous(breaks = c(15, seq(0, -50, by = -15))) +
  theme_pubr() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(file.path(out_dir, "panel_c_cage_tss_heatmap.jpg"), panel_c, width = 9, height = 6.75, units = "in")

message("Paneles B y C guardados en ", out_dir, ". Panel A (esquema manual) sigue pendiente.")
