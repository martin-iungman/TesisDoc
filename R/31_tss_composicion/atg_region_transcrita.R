# Control exploratorio (sin figuras ni tablas): presencia de ATG en la
# region transcrita incluida en el fragmento (TSS 0..+15).
# Busca ATG en la hebra del promotor con inicio entre 0 y +13 (asi queda
# completo dentro del fragmento de 252 pb, TSS de EPD en la posicion 237).
#   1. Promotores FP de la library, secuencia ensayada
#      (../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa).
#   2. Todo EPDnew hg38, secuencia genomica -236..+15 (misma ventana que la
#      library), con CGI_frag calculado igual que en build_prom_features.R.
# Resultado (oct 2026): ningun promotor FP tiene ATG en 0..+13, ni en la
# secuencia ensayada ni en la genomica -> fue un criterio de seleccion de la
# library, no una mutacion. En EPD completo ~11% tiene al menos uno.
# Requiere: data/processed/prom_df.tsv, data/raw/library.bed,
# data/external/EPD/human38_epdnew.bed. Solo imprime en consola.
# Run from the TesisDoc repo root.

suppressMessages({
  library(tidyverse)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
source("R/functions/tata_scan.R") # TATA_TSS_POS, TATA_FRAG_LEN

select <- dplyr::select
filter <- dplyr::filter

path_fa <- "../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa"
TSS_POS <- TATA_TSS_POS # 237
POS_ATG <- 0:13

# n_atg por secuencia (ATG con inicio en POS_ATG) y conteo por posicion
atg_scan <- function(dna) {
  m <- as.matrix(dna)
  L <- ncol(m)
  tri <- paste0(m[, 1:(L - 2)], m[, 2:(L - 1)], m[, 3:L])
  dim(tri) <- c(nrow(m), L - 2)
  win <- tri[, TSS_POS + POS_ATG, drop = FALSE] == "ATG"
  list(n_atg = rowSums(win), win = win)
}

resumen <- function(x) {
  tibble(
    N = nrow(x), con_ATG = sum(x$n_atg > 0), pct = round(100 * mean(x$n_atg > 0), 2),
    mas_de_uno = sum(x$n_atg > 1), pct_mas1 = round(100 * mean(x$n_atg > 1), 2)
  )
}

por_grupo <- function(d, extra = list()) {
  grupos <- c(
    list(Todos = d, `Con isla CpG` = filter(d, cgi %in% TRUE), `Sin isla CpG` = filter(d, cgi %in% FALSE)),
    extra
  )
  map(grupos, resumen) %>% bind_rows(.id = "grupo") %>% as.data.frame()
}

por_posicion <- function(win, cgi) {
  tibble(
    pos = POS_ATG, todos = colSums(win),
    con_CGI = colSums(win[cgi %in% TRUE, , drop = FALSE]),
    sin_CGI = colSums(win[cgi %in% FALSE, , drop = FALSE])
  ) %>% as.data.frame()
}

# --- 1. Library (secuencia ensayada) ---------------------------------------

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fa <- readDNAStringSet(path_fa)
fa <- fa[prom_df$seq_id]
stopifnot(all(width(fa) == TATA_FRAG_LEN))

lib_scan <- atg_scan(fa)
d_lib <- tibble(seq_id = names(fa), cgi = prom_df$CGI_frag, n_atg = lib_scan$n_atg)
cat("== Library FP (secuencia ensayada) ==\n")
print(por_grupo(d_lib))
print(por_posicion(lib_scan$win, d_lib$cgi))

# --- 2. EPDnew completo (secuencia genomica) -------------------------------

epd <- import.bed("data/external/EPD/human38_epdnew.bed")
epd <- epd[seqnames(epd) %in% standardChromosomes(Hsapiens)]
frag <- IRanges::promoters(epd, upstream = TSS_POS - 1, downstream = TATA_FRAG_LEN - TSS_POS + 1)
dna <- getSeq(Hsapiens, frag)
stopifnot(all(width(dna) == TATA_FRAG_LEN))

nf <- oligonucleotideFrequency(dna, width = 1)
cg <- vcountPattern("CG", dna)
oe <- ifelse(nf[, "C"] * nf[, "G"] > 0, cg * TATA_FRAG_LEN / (nf[, "C"] * nf[, "G"]), NA_real_)
lib <- import.bed("data/raw/library.bed")
lib_names <- str_remove(lib$name[startsWith(lib$name, "FP")], "^FP.{6}_")

epd_scan <- atg_scan(dna)
d_epd <- tibble(
  name = epd$name,
  cgi = (nf[, "C"] + nf[, "G"]) / TATA_FRAG_LEN >= 0.5 & !is.na(oe) & oe > 0.6,
  n_atg = epd_scan$n_atg,
  en_lib = name %in% lib_names
)
cat("\n== EPDnew hg38 (secuencia genomica -236..+15) ==\n")
print(por_grupo(d_epd, list(
  `FP en library` = filter(d_epd, en_lib), `EPD fuera de library` = filter(d_epd, !en_lib)
)))
print(por_posicion(epd_scan$win, d_epd$cgi))
cat("\nATG por promotor:\n")
print(table(d_epd$n_atg))
