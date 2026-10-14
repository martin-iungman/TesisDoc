# perfil_nucleotidos_library (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Perfil posicional de frecuencias nucleotidicas de los promotores FP de la
# library, en todos los promotores y con / sin isla CpG (CGI_frag):
# Panel A: frecuencia de A/C/G/T por posicion, -236 a +15
# Panel B: zoom en -50 a +15 (TATA-box en ~-30, iniciador en -1/0)
# Panel C: dinucleotido CpG (CG en p, p+1) observado y esperado por
#   composicion (frec. C en p x frec. G en p+1)
#
# Geometria: fragmento de 252 pb orientado de -236 a +15, TSS de EPD (= 0,
# sin +1 separado) en la posicion TATA_TSS_POS (237) de la secuencia. Se
# verifica la orientacion (pirimidina en -1, purina en 0) antes de graficar.
#
# Requiere: data/processed/prom_df.tsv (generado por
# R/00_prom_features/build_prom_features.R) y la secuencia de la library
# (../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa, solo
# los 252 pb del promotor, sin adaptadores). Run from the TesisDoc repo root.

library(Biostrings)
library(tidyverse)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/functions/tata_scan.R") # TATA_TSS_POS, TATA_FRAG_LEN

select <- dplyr::select
filter <- dplyr::filter

slug <- "perfil_nucleotidos_library"

path_fa <- "../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa"
BASES <- c("A", "C", "G", "T")
GRUPOS <- c("Todos", "Con isla CpG", "Sin isla CpG")
ZOOM <- c(-50, 15)
TATA_BANDA <- c(-31, -26)
col_bases <- c(A = "#7FB800", C = thesis_clr, G = "#D6741F", T = "#B5345A")
col_cpg <- c("Observado" = "#216869", "Esperado (C x G)" = "grey55")

# --- Datos -------------------------------------------------------------------

promoters_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

fa <- readDNAStringSet(path_fa)
fa <- fa[names(fa) %in% promoters_df$seq_id]
stopifnot(all(width(fa) == TATA_FRAG_LEN), all(promoters_df$seq_id %in% names(fa)))
con_n <- vcountPattern("N", fa) > 0
message("Secuencias con N excluidas del perfil: ", sum(con_n), " de ", length(fa))
fa <- fa[!con_n]
seq_mat <- as.matrix(fa)
posiciones <- seq_len(TATA_FRAG_LEN) - TATA_TSS_POS
cgi <- setNames(promoters_df$CGI_frag, promoters_df$seq_id)[names(fa)]
filas_grupo <- list(
  "Todos" = rep(TRUE, nrow(seq_mat)),
  "Con isla CpG" = cgi %in% TRUE,
  "Sin isla CpG" = cgi %in% FALSE
)

# --- Perfil por posicion -----------------------------------------------------

perfil <- imap(filas_grupo, function(filas, grupo) {
  m <- seq_mat[filas, , drop = FALSE]
  bases <- map(set_names(BASES), ~ colMeans(m == .x))
  cpg_obs <- colMeans(m[, -TATA_FRAG_LEN] == "C" & m[, -1] == "G")
  cpg_esp <- bases$C[-TATA_FRAG_LEN] * bases$G[-1]
  bind_rows(
    tibble(base = rep(BASES, each = TATA_FRAG_LEN), posicion = rep(posiciones, 4), frecuencia = unlist(bases)),
    tibble(base = "CpG", posicion = posiciones[-TATA_FRAG_LEN], frecuencia = cpg_obs),
    tibble(base = "CpG_esperado", posicion = posiciones[-TATA_FRAG_LEN], frecuencia = cpg_esp)
  ) %>%
    mutate(grupo = grupo, n = nrow(m))
}) %>%
  list_rbind() %>%
  mutate(grupo = factor(grupo, levels = GRUPOS)) %>%
  select(posicion, grupo, base, frecuencia, n)

# Verificacion de la orientacion: pirimidina en -1 y purina en 0.
orient <- perfil %>%
  filter(grupo == "Todos", posicion %in% c(-1, 0), base %in% BASES) %>%
  group_by(posicion) %>%
  summarise(pir = sum(frecuencia[base %in% c("C", "T")]), pur = sum(frecuencia[base %in% c("A", "G")]))
print(orient)
if (orient$pir[orient$posicion == -1] < 0.7 || orient$pur[orient$posicion == 0] < 0.75) {
  stop("Orientacion inesperada: en -1 deberian predominar pirimidinas y en 0 purinas")
}

# --- Paneles -----------------------------------------------------------------

minus <- function(x) sub("^-", "−", format(x, trim = TRUE))
perfil_bases <- perfil %>% filter(base %in% BASES)
y_max <- max(perfil_bases$frecuencia)
th_perfil <- theme_pubclean(base_size = 14) +
  theme(strip.text = element_text(face = "bold"), legend.position = "right")
marca_tata <- annotate("rect", xmin = TATA_BANDA[1] - 0.5, xmax = TATA_BANDA[2] + 0.5, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.4)
marca_tss <- geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")
x_completo <- scale_x_continuous(breaks = c(-200, -150, -100, -50, 0), labels = minus, expand = expansion(0.01))

panel_a <- ggplot(perfil_bases, aes(posicion, frecuencia, colour = base)) +
  marca_tata +
  marca_tss +
  geom_line(linewidth = 0.5) +
  facet_wrap(~grupo, nrow = 1) +
  scale_colour_manual(values = col_bases, name = "Base") +
  x_completo +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(c(0, 0.03))) +
  labs(x = "Posición respecto del TSS (pb)", y = "Frecuencia") +
  th_perfil

panel_b <- perfil_bases %>%
  filter(between(posicion, ZOOM[1], ZOOM[2])) %>%
  ggplot(aes(posicion, frecuencia, colour = base)) +
  marca_tata +
  marca_tss +
  geom_line(linewidth = 0.6) +
  geom_point(size = 0.8, show.legend = FALSE) +
  facet_wrap(~grupo, nrow = 1) +
  scale_colour_manual(values = col_bases, name = "Base") +
  scale_x_continuous(breaks = seq(-50, 10, 10), labels = minus, expand = expansion(0.02)) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(c(0, 0.03))) +
  labs(x = paste0("Posición respecto del TSS (pb), ", minus(ZOOM[1]), " a +", ZOOM[2]), y = "Frecuencia") +
  th_perfil

panel_c <- perfil %>%
  filter(base %in% c("CpG", "CpG_esperado")) %>%
  mutate(base = factor(ifelse(base == "CpG", "Observado", "Esperado (C x G)"), levels = names(col_cpg))) %>%
  ggplot(aes(posicion, frecuencia, colour = base, linetype = base)) +
  marca_tss +
  geom_line(linewidth = 0.5) +
  facet_wrap(~grupo, nrow = 1) +
  scale_colour_manual(values = col_cpg, name = "CpG") +
  scale_linetype_manual(values = c("solid", "22"), name = "CpG") +
  x_completo +
  scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  labs(x = "Posición respecto del TSS (pb)", y = "Frecuencia CpG") +
  th_perfil

# --- Save individual panels + combined figure -----------------------------

out_dir <- fig_dir(slug)
ggsave(file.path(out_dir, "panel_a_perfil_completo.jpg"), panel_a, width = 14, height = 4.5, units = "in")
ggsave(file.path(out_dir, "panel_b_perfil_zoom.jpg"), panel_b, width = 14, height = 4.5, units = "in")
ggsave(file.path(out_dir, "panel_c_cpg.jpg"), panel_c, width = 14, height = 3.5, units = "in")
write_tsv(perfil, file.path(out_dir, "perfil_nucleotidos.tsv"))

combined <- panel_a / panel_b / panel_c +
  plot_layout(heights = c(1, 1, 0.7), guides = "collect") +
  plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, paste0(slug, ".jpg")), combined, width = 14, height = 11, units = "in")
