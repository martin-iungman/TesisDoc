# Builds the shared promoter-feature table (prom_df) consumed by the
# figure scripts in R/06_composicion_library, R/07_motivos_promotor,
# R/08_conservacion_evolucion, R/09_coocurrencia, R/16_tf_contexto_endogeno.
# Run from the TesisDoc repo root. Output: data/processed/prom_df.tsv

library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(tidyverse)

# biomaRt/AnnotationDbi (loaded below) define S4 generics that shadow
# dplyr's select/filter/rename; pin them to the dplyr versions.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# --- Library import ---------------------------------------------------

lib <- rtracklayer::import.bed("data/raw/library.bed")
dna <- getSeq(Hsapiens, lib)
values(lib) <- tibble(seq_id = lib$name)

seq_data <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE)
prom_df <- seq_data %>% filter(seq_id %in% lib$seq_id)

# --- Ensembl gene IDs (needs internet access) ---------------------------
# Cached after the first successful query since Ensembl's API is
# occasionally down and this mapping barely changes between runs.

ids_cache_path <- "data/external/ensembl_gene_ids.tsv"
if (file.exists(ids_cache_path)) {
  ids_df <- read_tsv(ids_cache_path, show_col_types = FALSE)
} else {
  library(biomaRt)
  ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  ids_df <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
    filters = "hgnc_symbol",
    values = prom_df$gene_sym,
    mart = ensembl
  )
  ids_df <- ids_df %>%
    filter(chromosome_name %in% c(1:22, "X", "MT")) %>%
    distinct() %>%
    group_by(hgnc_symbol) %>%
    filter(ensembl_gene_id == sample(ensembl_gene_id, 1)) %>%
    select(-chromosome_name) %>%
    ungroup()
  write_tsv(ids_df, ids_cache_path)
}
prom_df <- left_join(prom_df, ids_df, by = c("gene_sym" = "hgnc_symbol"))

# --- Chromosome -----------------------------------------------------------

prom_df <- lib %>%
  as_tibble() %>%
  select(seq_id, seqnames) %>%
  left_join(prom_df, .) %>%
  rename(chromosome = seqnames)

# --- G+C content ------------------------------------------------------

G_C <- Biostrings::oligonucleotideFrequency(dna, width = 1, as.prob = TRUE)[, c("C", "G")] %>% rowSums()
lib$g_c <- round(G_C, 3)
prom_df <- elementMetadata(lib) %>% as_tibble() %>% dplyr::select(seq_id, g_c) %>% left_join(prom_df, .)

# --- CpG content --------------------------------------------------------

CpG <- Biostrings::oligonucleotideFrequency(dna, width = 2)[, c("CG", "GC")]
lib$CpG <- round(rowSums(CpG) / width(lib), 3)
prom_df <- elementMetadata(lib) %>% as_tibble() %>% select(seq_id, CpG) %>% left_join(prom_df, .)

# --- Observed vs expected CpG -------------------------------------------

prom_df <- Biostrings::oligonucleotideFrequency(dna, width = 1)[, c("C", "G")] %>%
  as_tibble() %>%
  mutate(seq_id = lib$seq_id, exp_cpg = (C * G) / (width(dna)^2)) %>%
  select(-c(C, G)) %>%
  left_join(prom_df, ., by = "seq_id")

# --- CpG islands (UCSC) ---------------------------------------------------
# Source: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz
# Already downloaded into data/external/cpgIslandExt.hg38.bed

cgi_df <- read_tsv(
  "data/external/cpgIslandExt.hg38.bed",
  col_names = c("seqname", "start", "end", "cgi_name", "length", "cpgNum", "gcNum", "perCpG", "perGC", "obsExp"),
  show_col_types = FALSE
)
cgi_lib <- GRanges(cgi_df)
message(sum(countOverlaps(cgi_lib, lib) > 1), " CGI overlap more than one promoter")
cgi_lib <- plyranges::join_overlap_left(cgi_lib, lib)
cgi_lib_inner <- plyranges::join_overlap_intersect(lib, cgi_lib)
cgi_lib_inner$width_cgi_overlap <- width(cgi_lib_inner)
# plyranges suffixes the seq_id column shared by both sides of the join
# (seq_id.x from lib, seq_id.y from cgi_lib); keep the lib one.
prom_df <- prom_df %>%
  left_join(as_tibble(cgi_lib_inner) %>% select(seq_id = seq_id.x, width_cgi_overlap)) %>%
  mutate(width_cgi_overlap = replace_na(width_cgi_overlap, 0)) %>%
  group_by(seq_id) %>%
  mutate(width_cgi_overlap = sum(width_cgi_overlap)) %>%
  distinct() %>%
  mutate(CGI = width_cgi_overlap > 100) %>%
  ungroup()

# --- EPD motifs (TATA, INR, CCAAT, GC-box) -------------------------------
# Only defined for EPD promoters (type == "promoter").
# Source: https://epd.epfl.ch/EPDnew_select.php subsets, already in data/external/EPD/

tata_epd <- import.bed("data/external/EPD/human38_epd_TATA.bed")
prom_df <- prom_df %>% mutate(TATA_EPD = ifelse(type == "promoter", ifelse(name %in% tata_epd$name, TRUE, FALSE), NA))

inr_epd <- import.bed("data/external/EPD/human38_epd_INR.bed")
prom_df <- prom_df %>% mutate(INR_EPD = ifelse(type == "promoter", ifelse(name %in% inr_epd$name, TRUE, FALSE), NA))

ccaat_epd <- import.bed("data/external/EPD/human38_epd_CCAAT.bed")
prom_df <- prom_df %>% mutate(CCAAT_EPD = ifelse(type == "promoter", ifelse(name %in% ccaat_epd$name, TRUE, FALSE), NA))

gc_epd <- import.bed("data/external/EPD/human38_epd_GCbox.bed")
prom_df <- prom_df %>% mutate(GCbox_EPD = ifelse(type == "promoter", ifelse(name %in% gc_epd$name, TRUE, FALSE), NA))

# --- TSS motif ------------------------------------------------------------
# Only defined for 252bp sequences (promoters with a canonical TSS window).

tss_seq <- lib %>%
  subset(width == 252) %>%
  getSeq(Hsapiens, .) %>%
  as.character() %>%
  str_extract(".{18}$") %>%
  str_extract("^.{5}")
tss_df <- tibble(
  seq_id = (lib %>% subset(width == 252))$seq_id,
  TSS_seq = tss_seq,
  INR_strong_TSS = str_detect(TSS_seq, "(T|C)CA(G|C)(T|A)"),
  TCT_TSS = str_detect(TSS_seq, ".(T|C)(C|G)(T|C)(C|T)"),
  CA_TSS = str_detect(TSS_seq, ".CA.."),
  CG_TSS = str_detect(TSS_seq, ".CG.."),
  TA_TSS = str_detect(TSS_seq, ".TA.."),
  TG_TSS = str_detect(TSS_seq, ".TG.."),
  GC_TSS = str_detect(TSS_seq, ".GC.."),
  PyPu_TSS = str_detect(TSS_seq, ".(C|T)(A|G).."),
  other_TSS = INR_strong_TSS + TCT_TSS + CA_TSS + CG_TSS + TA_TSS + TG_TSS == 0
)
prom_df <- left_join(prom_df, tss_df, by = "seq_id")

# --- PhyloP conservation score -------------------------------------------
# hg38.phyloP100way_summary_lib.bed is small (1MB), copied to data/external/.
# The per-bp version is 280MB; see R/00_prom_features/heavy_data_paths.R.

source("R/00_prom_features/heavy_data_paths.R")

phylo100 <- import.bed("data/external/hg38.phyloP100way_summary_lib.bed") %>% sort()
lib <- sort(lib)
lib$phylo100 <- as.numeric(phylo100$name) %>% round(4)
# lib accumulates columns (g_c, CpG...) added by earlier sections; only
# bring in the new one to avoid re-joining duplicates.
prom_df <- left_join(prom_df, as_tibble(values(lib)) %>% select(seq_id, phylo100), by = "seq_id")

phylo1pb <- rtracklayer::import.bed(path_phylo_1bp)
phylo1pb$score <- phylo1pb$name
phylo1pb$name <- NULL
phylo1pb <- plyranges::join_overlap_inner(phylo1pb, lib)
phylo1pb_prom <- subset(phylo1pb, str_detect(seq_id, "^FP.{6}_"))
lib_df <- as_tibble(lib)
names(lib_df) <- paste0(names(lib_df), "_seq")
values(phylo1pb_prom) <- left_join(as_tibble(values(phylo1pb_prom)), lib_df, by = c("seq_id" = "seq_id_seq"))
strand(phylo1pb_prom) <- phylo1pb_prom$strand_seq
phylo1pb_prom$seq <- getSeq(Hsapiens, phylo1pb_prom)
phylo1pb_prom_df <- distinct(as_tibble(phylo1pb_prom)) %>%
  mutate(
    relpos = ifelse(strand_seq == "-", abs(start_seq - start), abs(start - end_seq)),
    score = as.numeric(score)
  )

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(phylo1pb_prom_df, "data/processed/phylo100_1pb_lib.tsv")
message("Wrote data/processed/phylo100_1pb_lib.tsv with ", nrow(phylo1pb_prom_df), " rows")

range_sum_phylo <- phylo1pb_prom_df %>%
  mutate(range = ifelse((relpos - 16) <= 50, "close", ifelse((relpos - 16) <= 150, "intermediate", "far"))) %>%
  group_by(seq_id, range) %>%
  summarise(mean_range_phylo100 = mean(score, na.rm = TRUE), .groups = "drop")
prom_df <- range_sum_phylo %>%
  pivot_wider(values_from = mean_range_phylo100, names_from = range, names_prefix = "phylop100_") %>%
  left_join(prom_df, ., by = "seq_id")

# --- Human-mouse functional turnover (Young et al. 2015) -----------------

epd_hg19 <- import.bed("data/external/EPD/EPD_hg19.bed")
young <- read_tsv(
  "data/external/Sup1_mouse_human_Young2015.txt",
  skip = 2,
  col_names = c("seqnames", "start", "end", "human_pos", "mouse_pos", "turnover"),
  show_col_types = FALSE
)
young$strand <- str_extract(young$human_pos, ".$")
hg19_young_lib <- epd_hg19 %>%
  subset(name %in% prom_df$name) %>%
  plyranges::join_overlap_left_directed(GRanges(young))
na_proms <- subset(hg19_young_lib, is.na(turnover))
extra_proms <- epd_hg19 %>%
  subset(name %in% na_proms$name) %>%
  promoters(upstream = 235, downstream = 16) %>%
  plyranges::join_overlap_left_directed(GRanges(young)) %>%
  subset(!is.na(turnover)) %>%
  values() %>%
  as_tibble() %>%
  select(-human_pos, -mouse_pos) %>%
  distinct() %>%
  add_count(name) %>%
  filter(n == 1) %>%
  select(-n, -score)
turnover_df <- values(hg19_young_lib) %>%
  as_tibble() %>%
  filter(!is.na(turnover)) %>%
  select(-human_pos, -mouse_pos, -score) %>%
  bind_rows(extra_proms)

prom_df <- prom_df %>% left_join(turnover_df, by = "name")

# --- Write output -----------------------------------------------------

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(prom_df, "data/processed/prom_df.tsv")
message("Wrote data/processed/prom_df.tsv with ", nrow(prom_df), " rows and ", ncol(prom_df), " columns")
