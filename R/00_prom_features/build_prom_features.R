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

# --- Repeated elements (RepeatMasker) -------------------------------------
# strand is ignored here (matches the original pipeline): RepeatMasker's
# "cons_strand" column uses +/C, not +/-/*, so GRanges() leaves strand="*"
# for all rows and the join below is not strand-aware.

repeatmasker <- read_delim(path_repeatmasker, col_names = TRUE, delim = " ", show_col_types = FALSE)
repeatmasker_gr <- GRanges(repeatmasker) %>%
  plyranges::join_overlap_intersect(lib, .) %>%
  subset(width > 10)

dna_repeat <- getSeq(Hsapiens, repeatmasker_gr)
G_C_repeat <- Biostrings::oligonucleotideFrequency(dna_repeat, width = 1, as.prob = TRUE)[, c("C", "G")] %>% rowSums()
repeatmasker_gr$G_C_repeat <- round(G_C_repeat, 3)
repeatmasker_df <- repeatmasker_gr %>%
  as_tibble() %>%
  mutate(
    repeat_family = ifelse(repeat_class %in% c("Low_complexity", "Simple_repeat"), "LCR", "TE"),
    repeat_superclass = str_remove(repeat_class, "/.+$")
  ) %>%
  rename(repeat_overlap = width)

LCR_df <- repeatmasker_df %>%
  filter(repeat_family == "LCR") %>%
  select(seq_id, repeat_overlap, G_C_repeat) %>%
  group_by(seq_id) %>%
  summarise(LCR_overlap = sum(repeat_overlap), G_C_LCR = sum(G_C_repeat * repeat_overlap) / LCR_overlap)
TE_df <- repeatmasker_df %>%
  filter(repeat_family == "TE") %>%
  select(seq_id, repeat_overlap, G_C_repeat, repeat_superclass) %>%
  group_by(seq_id) %>%
  summarise(TE_overlap = sum(repeat_overlap), G_C_TE = sum(G_C_repeat * repeat_overlap) / TE_overlap, TE_superclass = unique(repeat_superclass)[[1]])

prom_df <- prom_df %>%
  left_join(LCR_df, by = "seq_id") %>%
  left_join(TE_df, by = "seq_id") %>%
  mutate(
    across(c(LCR_overlap, TE_overlap, G_C_TE, G_C_LCR), ~ replace_na(.x, 0)),
    TE_superclass = TE_superclass %>% fct_lump_min(100)
  )
rm(repeatmasker, repeatmasker_gr, repeatmasker_df, dna_repeat)

# --- Enhancers (FANTOM5, distance-window overlap counts) ------------------

enhancers <- import.bed("data/external/F5.hg38.enhancers.bed")
lib10kb <- promoters(lib, upstream = 5000, downstream = 5000)
lib50kb <- promoters(lib, upstream = 25000, downstream = 25000)
lib100kb <- promoters(lib, upstream = 50000, downstream = 50000)

lib$enh10kb <- countOverlaps(lib10kb, enhancers)
lib$enh50kb <- countOverlaps(lib50kb, enhancers)
lib$enh100kb <- countOverlaps(lib100kb, enhancers)
prom_df <- left_join(prom_df, values(lib) %>% as_tibble() %>% select(seq_id, enh10kb, enh50kb, enh100kb), by = "seq_id")

# --- Chromatin accessibility (DNase, ENCODE) -------------------------------

dnase <- rtracklayer::import.bw(path_dnase, as = "GRanges")
lib2 <- lib
strand(lib2) <- "*"
# restrict the genome-wide DNase track down to the library's chromosomes
seqlevels(dnase, pruning.mode = "coarse") <- seqlevels(lib2)
ol <- plyranges::find_overlaps(dnase, lib2)
long_dnase <- as_tibble(ol) %>%
  select(-strand) %>%
  left_join(as_tibble(lib) %>% select(seq_id, strand), by = "seq_id") %>%
  filter(str_detect(seq_id, "^FP")) %>%
  group_by(seq_id) %>%
  mutate(range_order = (ifelse(strand == "-", row_number(-start), row_number(start)) - 1) * 25)
wide_dnase <- long_dnase %>%
  select(seq_id, range_order, score) %>%
  group_by(seq_id) %>%
  mutate(mean_dnase = mean(score)) %>%
  pivot_wider(names_from = range_order, values_from = score, names_prefix = "dnase_") %>%
  select(-dnase_275)
prom_df <- left_join(prom_df, wide_dnase, by = "seq_id")
rm(dnase, lib2, ol, long_dnase, wide_dnase)

# --- Cis-regulatory modules (CRM), tissue specificity, shape, HEK293 -----
# EXCEPTION: reads transcriptional_library/Analysis/Tables directly - see
# R/00_prom_features/analysis_tables_exceptions.R for why and what's pending.

source("R/00_prom_features/analysis_tables_exceptions.R")

crm <- read_tsv(path_library_remap_crm,
  col_names = c(
    "seqnames", "start", "end", "seq_id", "score", "strand", "seqnames_peak", "start_peak",
    "end_peak", "peak", "N_TF_CRM", "strand_peak", "thickstart_peak", "thickend_peak", "rgb"
  ),
  show_col_types = FALSE
)
prom_df <- crm %>%
  select(seq_id, N_TF_CRM) %>%
  group_by(seq_id) %>%
  summarise(N_TF_CRM = sum(N_TF_CRM)) %>%
  right_join(prom_df, by = "seq_id") %>%
  mutate(N_TF_CRM = replace_na(N_TF_CRM, 0))
rm(crm)

# Tissue specificity: classification (discrete) + Gini index (numeric).
# sample_specific = >10TPM in some sample, active (>=1TPM) in <10 samples
# group_enrichment = 50% of total activity accumulated in <=7 samples
# sample_enhanced = some sample >=4x the mean of active samples
# all_samples_detected = none of the above, active in >=80% of samples
# mixed = everything else
sample_activity <- read_tsv(path_sample_cage_activity, show_col_types = FALSE)
sample_activity_lib <- sample_activity %>% filter(name %in% prom_df$name)
n_samples <- unique(sample_activity_lib$sample) %>% length()
obj <- sample_activity_lib %>%
  mutate(tpm = counts * 1e6 / libsize) %>%
  group_by(name) %>%
  mutate(active = tpm > 1, n_detected = sum(active), sample_enhancement = tpm / (sum(tpm) / n_detected)) %>%
  ungroup()
enrichment <- obj %>%
  group_by(name) %>%
  arrange(name, desc(tpm)) %>%
  mutate(i = row_number(), group_enrichment = cumsum(tpm) / (sum(tpm) - cumsum(tpm))) %>%
  filter(active, group_enrichment > 1, i <= (0.1 * n_samples)) %>%
  filter(i == min(i)) %>%
  mutate(sample_specificity_class = "group_enrichment")
enhanced <- obj %>% filter(sample_enhancement >= 4, active == TRUE)
class_long <- obj %>%
  mutate(sample_specificity_class = ifelse(name %in% (enrichment$name[enrichment$sample_specificity_class == "group_enrichment"]), "group_enrichment",
    ifelse(n_detected == 0, "non_detected",
      ifelse(name %in% enhanced$name, "sample_enhanced",
        ifelse(n_detected > 0.8 * n_samples, "all_samples_detected", "mixed")
      )
    )
  ))
class_short <- class_long %>%
  select(name, sample_specificity_class, n_detected) %>%
  distinct()
prom_df <- left_join(prom_df, class_short %>% rename(active_fantom5_samples = n_detected), by = "name")

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
gini_values <- obj %>%
  filter(n_detected > 0) %>%
  group_by(name) %>%
  summarise(sample_specificity_gini = gini(tpm), sum_tpm_fantom5 = sum(tpm))
prom_df <- left_join(prom_df, gini_values, by = "name")
rm(sample_activity, sample_activity_lib, obj, enrichment, enhanced, class_long, class_short, gini_values)

# Promoter shape (CAGEr interquantile width): Broad if >25bp else Narrow.
# shape_merged's "name" holds short gene-symbol-like names (e.g.
# "KLHL17_1"), matching prom_df's "name" column - NOT prom_df's seq_id.
shape <- read_tsv(path_shape_merged, show_col_types = FALSE) %>%
  group_by(name) %>%
  filter(dominant_ctss.score == max(dominant_ctss.score)) %>%
  ungroup() %>%
  mutate(shape_class = ifelse(interquantile_width > 25, "Broad", "Narrow")) %>%
  select(name, shape_class, interquantile_width) %>%
  distinct(name, .keep_all = TRUE)
prom_df <- left_join(prom_df, shape, by = "name")
rm(shape)

# Endogenous HEK293 activity (CAGEr).
hek_df <- read_tsv(path_hek_cager, show_col_types = FALSE) %>%
  group_by(name) %>%
  summarise(hek_tpm = sum(score.y, na.rm = TRUE), hek_interq_width = sum(interquantile_width))
prom_df <- left_join(prom_df, hek_df, by = c("seq_id" = "name"))
rm(hek_df)

# --- Write output -----------------------------------------------------

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(prom_df, "data/processed/prom_df.tsv")
message("Wrote data/processed/prom_df.tsv with ", nrow(prom_df), " rows and ", ncol(prom_df), " columns")
