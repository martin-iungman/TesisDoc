# Metilacion de CpG por WGBS en HEK293 y HEK293T sobre cada promotor de la
# library. Complementa build_metilacion_hek293.R (RRBS de ENCODE), que casi no cubre
# promotores pobres en CpG ni LTR, justo el grupo de silentes sin marcas.
#
# Fuente: GEO, solo muestras control (archivos procesados por los autores):
#   HEK293:  GSE243785 (control con vehiculo, rep 1 y 2; Bismark .cov, hg19)
#            GSE233830 (TNRC18 WT, rep 1 a 3; cytosine report de DRAGEN, hg19)
#   HEK293T: GSE270034 (WT, rep 1 y 3; CpG con >= 5 lecturas, hg38)
# Las ventanas de la library se llevan a hg19 con la cadena de UCSC hg38ToHg19 (se
# descartan las que cambian de largo en mas de 10%), en lugar de levantar
# decenas de millones de CpG.
#
# Por promotor y muestra: lecturas metiladas / totales sumando todos los CpG de la
# ventana, y numero de CpG con >= 5 lecturas. Ventanas: el fragmento de la library
# y 1 kb centrada en el.
#
# Descarga ~1.8 GB a un directorio temporal, que se borra al final (o a WGBS_DIR,
# si esta definido, y ahi se conserva). No se incorpora a prom_df: lo leen
# directamente los scripts que lo usan (silentes_ltr y el exploratorio de
# silentes). Salida:
# data/processed/metilacion_hek293_wgbs_v1.tsv. Requiere data/raw/library.bed.
# Run from the TesisDoc repo root.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
})

options(timeout = 7200)
MIN_LECTURAS_CPG <- 5
geo <- "https://ftp.ncbi.nlm.nih.gov/geo/samples"
muestras <- tribble(
  ~muestra, ~serie, ~linea, ~genoma, ~formato, ~ruta,
  "GSM7795652", "GSE243785", "HEK293", "hg19", "cov", "GSM7795nnn/GSM7795652/suppl/GSM7795652_WGBS_control_rep1.cov.gz",
  "GSM7795653", "GSE243785", "HEK293", "hg19", "cov", "GSM7795nnn/GSM7795653/suppl/GSM7795653_WGBS_control_rep2.cov.gz",
  "GSM7437599", "GSE233830", "HEK293", "hg19", "cpg_report", "GSM7437nnn/GSM7437599/suppl/GSM7437599_22144D-01-01.CpG_report.txt.gz",
  "GSM7437600", "GSE233830", "HEK293", "hg19", "cpg_report", "GSM7437nnn/GSM7437600/suppl/GSM7437600_22144D-01-02.CpG_report.txt.gz",
  "GSM7437601", "GSE233830", "HEK293", "hg19", "cpg_report", "GSM7437nnn/GSM7437601/suppl/GSM7437601_22144D-01-03.CpG_report.txt.gz",
  "GSM8333112", "GSE270034", "HEK293T", "hg38", "mc5x", "GSM8333nnn/GSM8333112/suppl/GSM8333112_HEK293T_WT_CG_mC_5x.txt.gz",
  "GSM8543460", "GSE270034", "HEK293T", "hg38", "mc5x", "GSM8543nnn/GSM8543460/suppl/GSM8543460_HEK293T_WT_CG_mC_5x_CG.rep3.txt.gz"
)

dir_dl <- Sys.getenv("WGBS_DIR", unset = tempfile("wgbs_hek293_"))
dir.create(dir_dl, showWarnings = FALSE, recursive = TRUE)
conservar <- nzchar(Sys.getenv("WGBS_DIR"))
bajar <- function(url, destino) {
  if (file.exists(destino)) {
    return(destino)
  }
  # curl retoma la descarga donde se corto (el FTP de NCBI corta archivos grandes)
  parcial <- paste0(destino, ".part")
  ok <- FALSE
  for (intento in 1:10) {
    estado <- system2("curl", c("-s", "-L", "--fail", "-C", "-", "--retry", "5", "--retry-delay", "10", "-o", shQuote(parcial), shQuote(url)))
    if (estado == 0) {
      ok <- TRUE
      break
    }
    Sys.sleep(15)
  }
  stopifnot(ok)
  file.rename(parcial, destino)
  destino
}

# --- Ventanas en hg38 y hg19 ----------------------------------------------------

lib <- read_tsv("data/raw/library.bed", col_names = FALSE, show_col_types = FALSE, col_select = 1:4) %>%
  distinct(X4, .keep_all = TRUE)
frag <- GRanges(lib$X1, IRanges(lib$X2 + 1, lib$X3), seq_id = lib$X4, ventana = "Fragmento")
ventanas38 <- c(frag, resize(frag, 1000, fix = "center") %>% `mcols<-`(value = DataFrame(seq_id = frag$seq_id, ventana = "1 kb")))

cadena_gz <- bajar("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz", file.path(dir_dl, "hg38ToHg19.over.chain.gz"))
cadena_txt <- file.path(dir_dl, "hg38ToHg19.over.chain")
writeLines(readLines(gzfile(cadena_gz)), cadena_txt)
cadena <- import.chain(cadena_txt)
lev <- liftOver(ventanas38, cadena)
largo <- sum(width(lev))
unica <- elementNROWS(GenomicRanges::reduce(lev)) == 1
ok <- unica & abs(largo - width(ventanas38)) <= 0.1 * width(ventanas38)
ventanas19 <- unlist(range(lev[ok]))
mcols(ventanas19) <- mcols(ventanas38)[ok, ]
cat("Ventanas llevadas a hg19:", sum(ok), "de", length(ventanas38), "\n")

# --- Lectura de cada formato ------------------------------------------------------

leer <- function(f, formato) {
  if (formato == "cov") {
    d <- fread(f, select = c(1, 2, 5, 6), col.names = c("chr", "pos", "met", "nomet"))
  } else if (formato == "cpg_report") {
    d <- fread(f, select = c(1, 2, 4, 5), col.names = c("chr", "pos", "met", "nomet"))
  } else {
    d <- fread(f, select = c(1, 2, 6, 7), col.names = c("chr", "pos", "total", "met"))
    d[, nomet := total - met][, total := NULL]
  }
  d[met + nomet > 0]
}

resumir <- function(d, ventanas) {
  cpg <- GRanges(d$chr, IRanges(d$pos, width = 1))
  ov <- findOverlaps(ventanas, cpg)
  dt <- data.table(i = queryHits(ov), met = d$met[subjectHits(ov)], tot = d$met[subjectHits(ov)] + d$nomet[subjectHits(ov)])
  r <- dt[, .(lecturas = sum(tot), metilacion = 100 * sum(met) / sum(tot), n_cpg5 = sum(tot >= MIN_LECTURAS_CPG)), by = i]
  as_tibble(r) %>%
    mutate(seq_id = ventanas$seq_id[i], ventana = ventanas$ventana[i]) %>%
    select(-i)
}

res <- pmap(muestras, function(muestra, serie, linea, genoma, formato, ruta) {
  cat(muestra, serie, linea, format(Sys.time(), "%H:%M"), "\n")
  f <- bajar(file.path(geo, ruta), file.path(dir_dl, basename(ruta)))
  d <- leer(f, formato)
  cat("  CpG con lecturas:", nrow(d), "| lecturas medias:", round(mean(d$met + d$nomet), 1), "\n")
  r <- resumir(d, if (genoma == "hg19") ventanas19 else ventanas38)
  rm(d)
  gc()
  if (!conservar) unlink(f)
  mutate(r, muestra = muestra, serie = serie, linea = linea)
}) %>%
  list_rbind() %>%
  mutate(metilacion = round(metilacion, 2))

write_tsv(res, "data/processed/metilacion_hek293_wgbs_v1.tsv")
if (!conservar) unlink(dir_dl, recursive = TRUE)
cat("Listo\n")
