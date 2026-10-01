# Orientation Index (OI) de cada promotor FP de la library, a partir de
# datos de PRO-cap (HEK293, ENCODE) - Core & Martins et al. Supp. Fig. 5:
#   OI = 2 x max(Rp, Rm) / (Rp + Rm) - 1
#   Rp, Rm: reads de PRO-cap en sentido/antisentido en los 250pb
#   downstream/upstream del TSS
#   OI in [0, 1]: 0 = bidireccional/divergente, 1 = unidireccional
# Clasificacion: OI >= 0.9 -> unidireccional; OI < 0.9 y señal suficiente
# -> bidireccional; señal insuficiente (Rp_total < 10) -> no clasificado.
#
# Portado de transcriptional_library/Analysis/scripts/promoter_orientation.R
# (pasos 1-4 de ese script; se omite la validacion contra los BED
# pre-clasificados de ENCODE - paso 5 - y todo lo de la seccion 8 de ese
# script, que no corresponde a esta figura). Agrega la columna
# `orientacion` (unidireccional/bidireccional/NA) a prom_df.tsv.
#
# Requiere data/raw/library.bed (promotores FP de la library) y los
# bigWig de PRO-cap referenciados en heavy_data_paths.R (no copiados,
# se leen directo de transcriptional_library). Corre sobre TODOS los
# promotores FP de la library, no solo los presentes en prom_df.tsv, asi
# que el join posterior deja NA para los que no tengan promoter FP
# correspondiente (no deberia pasar, pero por robustez).
#
# Run from the TesisDoc repo root.

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
source("R/00_prom_features/heavy_data_paths.R")

OI_WIN <- 250L
PSEUDOCOUNT <- 1
OI_THRESH <- 0.9
MIN_SENSE <- 10L

# --- 1. Promotores FP de la library ----------------------------------------

lib <- rtracklayer::import.bed("data/raw/library.bed")
lib_fp <- lib[startsWith(lib$name, "FP")]
stopifnot(all(width(lib_fp) == 252L))
message("Promotores FP: ", length(lib_fp))

plus_idx <- as.character(strand(lib_fp)) == "+"

# TSS (1-based, ventana BED de 252pb): + strand start+236, - strand end-236
tss_pos <- ifelse(plus_idx, start(lib_fp) + 236L, end(lib_fp) - 236L)

# --- 2. Ventanas de la OI: sentido (downstream) vs. antisentido (upstream) --

oi_down_gr <- GRanges(
  seqnames = seqnames(lib_fp),
  ranges = IRanges(
    start = ifelse(plus_idx, tss_pos, tss_pos - OI_WIN + 1L),
    end   = ifelse(plus_idx, tss_pos + OI_WIN - 1L, tss_pos)
  ),
  strand = strand(lib_fp)
)
oi_up_gr <- GRanges(
  seqnames = seqnames(lib_fp),
  ranges = IRanges(
    start = ifelse(plus_idx, tss_pos - OI_WIN, tss_pos + 1L),
    end   = ifelse(plus_idx, tss_pos - 1L, tss_pos + OI_WIN)
  ),
  strand = strand(lib_fp)
)

# --- 3. Extraer señal de los bigWig -----------------------------------------

bw_region_sum <- function(bw_path, regions) {
  hits <- tryCatch(
    rtracklayer::import(BigWigFile(bw_path), which = regions),
    error = function(e) GRanges()
  )
  result <- rep(0, length(regions))
  if (length(hits) == 0L) return(result)
  hits$score <- abs(score(hits))
  ov <- findOverlaps(regions, hits)
  s <- tapply(score(hits)[subjectHits(ov)], queryHits(ov), sum)
  result[as.integer(names(s))] <- as.numeric(s)
  result
}

message("Extrayendo señal de PRO-cap en las ventanas de OI...")
pd_r1 <- bw_region_sum(path_procap_plus_r1, oi_down_gr)
md_r1 <- bw_region_sum(path_procap_minus_r1, oi_down_gr)
pu_r1 <- bw_region_sum(path_procap_plus_r1, oi_up_gr)
mu_r1 <- bw_region_sum(path_procap_minus_r1, oi_up_gr)
pd_r2 <- bw_region_sum(path_procap_plus_r2, oi_down_gr)
md_r2 <- bw_region_sum(path_procap_minus_r2, oi_down_gr)
pu_r2 <- bw_region_sum(path_procap_plus_r2, oi_up_gr)
mu_r2 <- bw_region_sum(path_procap_minus_r2, oi_up_gr)

rp_r1 <- ifelse(plus_idx, pd_r1, md_r1)
rm_r1 <- ifelse(plus_idx, mu_r1, pu_r1)
rp_r2 <- ifelse(plus_idx, pd_r2, md_r2)
rm_r2 <- ifelse(plus_idx, mu_r2, pu_r2)

# --- 4. Calcular OI (pooleando las 2 replicas tecnicas) y clasificar ------

oi_fun <- function(rp, rm, eps = PSEUDOCOUNT) 2 * pmax(rp + eps, rm + eps) / (rp + rm + 2 * eps) - 1
oi_pooled <- oi_fun(rp_r1 + rp_r2, rm_r1 + rm_r2)
rp_total <- rp_r1 + rp_r2
rm_total <- rm_r1 + rm_r2
total_signal <- rp_total + rm_total + 0 # (ya son sumas de reps; sin pseudocount)

orientation_class <- case_when(
  rp_total < MIN_SENSE ~ NA_character_,
  oi_pooled >= OI_THRESH ~ "unidireccional",
  TRUE ~ "bidireccional"
)

orientation_df <- tibble(
  seq_id = lib_fp$name,
  OI = oi_pooled,
  Rp_total = rp_total,
  Rm_total = rm_total,
  total_signal = total_signal,
  orientacion = orientation_class
)

message("─── Clasificación de orientación ───────────────────")
print(table(orientation_df$orientacion, useNA = "always"))
message("──────────────────────────────────────────────────────")

write_tsv(orientation_df, "data/processed/promoter_orientation.tsv")
message("Guardado data/processed/promoter_orientation.tsv")

# --- 5. Agregar `orientacion` a prom_df.tsv ---------------------------------
# PENDING (igual que hek_tpm/mean_dnase en su momento): integrar este
# calculo directamente a build_prom_features.R en vez de un join posterior.

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
if ("orientacion" %in% names(prom_df)) prom_df <- select(prom_df, -orientacion)
prom_df <- left_join(prom_df, orientation_df %>% select(seq_id, orientacion), by = "seq_id")
write_tsv(prom_df, "data/processed/prom_df.tsv")
message("Columna `orientacion` agregada a data/processed/prom_df.tsv")
