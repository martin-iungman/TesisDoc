# G4 observados in vitro (G4-seq, Marsico et al. 2019, GSE110582) en cada
# fragmento FP de la library, separados por hebra respecto del promotor.
# Dos condiciones: K+ (GSM3003539, umbral 25 del archivo) y K+ + PDS
# (GSM3003540, umbral 35). Los llamados (OQ) estan en hg19 y separados por
# hebra del genoma (archivos plus / minus); la columna 4 es el puntaje de
# G4-seq (diferencia de mismatches entre la lectura en K+ (o K+ + PDS) y
# la de Li+). Se pasa la library de hg38 a hg19 con liftOver (se descartan
# los fragmentos que no quedan contiguos con 252 pb) y un OQ cuenta para el
# fragmento si su centro cae dentro de los 252 pb.
#
# Hebra: verificado empiricamente que los OQ del archivo "plus" tienen las G
# en la hebra - del genoma (en la secuencia orientada de un promotor en +
# son ricos en C: 45% C vs. 26% G con K+) y coinciden con el regex canonico
# de la hebra opuesta. Entonces un OQ "plus" es G4 en la hebra codificante
# de un promotor en - y en la hebra molde de uno en +. La asignacion se
# chequea al final del script (G vs. C de cada OQ en la secuencia orientada).
#
# Salidas:
#   data/processed/g4seq_lib_v1.tsv: por fragmento x condicion, numero de
#     OQ por hebra, puntaje maximo y posicion (centro respecto del TSS) del
#     OQ mas cercano al TSS; liftover_ok = el fragmento paso a hg19.
#   data/processed/g4seq_hits_lib_v1.tsv: un OQ por fila (para histogramas).
# Requiere data/raw/library.bed y data/external/G4seq/ (ver
# data/external/download_scripts/README.md).
# Run from the TesisDoc repo root.

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)

dir_g4seq <- "data/external/G4seq"
ARCHIVOS <- tribble(
  ~condicion, ~hebra_genoma, ~archivo,
  "K+", "+", "GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz",
  "K+", "-", "GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz",
  "K+ + PDS", "+", "GSM3003540_Homo_all_w15_th-1_plus.hits.max.PDS.w50.35.bed.gz",
  "K+ + PDS", "-", "GSM3003540_Homo_all_w15_th-1_minus.hits.max.PDS.w50.35.bed.gz"
)
TSS_OFFSET <- 236L # TSS = start + 236 (+) o end - 236 (-), coordenadas 1-based

lib <- rtracklayer::import.bed("data/raw/library.bed")
lib_fp <- lib[startsWith(lib$name, "FP")]
stopifnot(all(width(lib_fp) == 252L))

# liftOver hg38 -> hg19 (import.chain no lee .gz)
chain_tmp <- tempfile(fileext = ".chain")
writeLines(readLines(gzfile(file.path(dir_g4seq, "hg38ToHg19.over.chain.gz"))), chain_tmp)
lifted <- rtracklayer::liftOver(lib_fp, rtracklayer::import.chain(chain_tmp)) %>% GenomicRanges::reduce(min.gapwidth = 1L)
ok <- lengths(lifted) == 1L
ok[ok] <- width(unlist(lifted[ok])) == 252L
lib19 <- unlist(lifted[ok])
lib19$seq_id <- lib_fp$name[ok]
message("Fragmentos que pasan a hg19 contiguos con 252 pb: ", sum(ok), " / ", length(lib_fp))

hits <- pmap_dfr(ARCHIVOS, function(condicion, hebra_genoma, archivo) {
  read_tsv(file.path(dir_g4seq, archivo), col_names = c("chr", "start0", "end", "score"),
           col_types = "ciid", progress = FALSE) %>%
    mutate(condicion = condicion, hebra_genoma = hebra_genoma)
}) %>%
  mutate(centro = (start0 + 1 + end) / 2)

# centro del OQ dentro del fragmento
hits_gr <- GRanges(hits$chr, IRanges(floor(hits$centro), width = 1))
ov <- findOverlaps(hits_gr, lib19, ignore.strand = TRUE)
frag_strand <- as.character(strand(lib19))
hits_lib <- hits[queryHits(ov), ] %>%
  mutate(
    seq_id = lib19$seq_id[subjectHits(ov)],
    frag_strand = frag_strand[subjectHits(ov)],
    tss = if_else(frag_strand == "+", start(lib19)[subjectHits(ov)] + TSS_OFFSET, end(lib19)[subjectHits(ov)] - TSS_OFFSET),
    pos = if_else(frag_strand == "+", centro - tss, tss - centro),
    largo = end - start0,
    hebra = if_else(hebra_genoma != frag_strand, "Codificante", "Molde") # ver encabezado
  ) %>%
  select(seq_id, condicion, hebra, pos, largo, score, chr_hg19 = chr, start0_hg19 = start0, end_hg19 = end)

# chequeo de hebra: en la secuencia orientada del fragmento, los OQ de la
# hebra codificante deben ser ricos en G y los de la molde ricos en C
fa <- Biostrings::readDNAStringSet("../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa")
chk <- hits_lib %>%
  filter(seq_id %in% names(fa)) %>%
  mutate(
    s = substring(as.character(fa[seq_id]), pmax(1, round(pos + 237 - largo / 2)), pmin(252, round(pos + 237 + largo / 2))),
    G = str_count(s, "G") / nchar(s), C = str_count(s, "C") / nchar(s)
  ) %>%
  group_by(condicion, hebra) %>%
  summarise(n = n(), G = mean(G), C = mean(C), .groups = "drop")
print(chk)
stopifnot(with(chk, all(G[hebra == "Codificante"] > C[hebra == "Codificante"]), all(C[hebra == "Molde"] > G[hebra == "Molde"])))

cerca <- function(p) if (length(p)) p[which.min(abs(p))] else NA_real_
por_frag <- expand_grid(seq_id = lib_fp$name, condicion = unique(ARCHIVOS$condicion)) %>%
  left_join(
    hits_lib %>%
      group_by(seq_id, condicion) %>%
      summarise(
        n_oq_cod = sum(hebra == "Codificante"), n_oq_molde = sum(hebra == "Molde"),
        score_max_cod = suppressWarnings(max(score[hebra == "Codificante"])),
        score_max_molde = suppressWarnings(max(score[hebra == "Molde"])),
        pos_oq_cod = cerca(pos[hebra == "Codificante"]), pos_oq_molde = cerca(pos[hebra == "Molde"]),
        .groups = "drop"
      ),
    by = c("seq_id", "condicion")
  ) %>%
  mutate(
    liftover_ok = seq_id %in% lib19$seq_id,
    across(c(n_oq_cod, n_oq_molde), ~ if_else(liftover_ok, replace_na(.x, 0L), NA_integer_)),
    across(c(score_max_cod, score_max_molde), ~ if_else(is.finite(.x), .x, NA_real_))
  )

print(por_frag %>% filter(liftover_ok) %>% group_by(condicion) %>%
  summarise(n = n(), pct_cod = 100 * mean(n_oq_cod > 0), pct_molde = 100 * mean(n_oq_molde > 0)))
write_tsv(por_frag, "data/processed/g4seq_lib_v1.tsv")
write_tsv(hits_lib, "data/processed/g4seq_hits_lib_v1.tsv")
message("Guardado data/processed/g4seq_lib_v1.tsv y g4seq_hits_lib_v1.tsv")
