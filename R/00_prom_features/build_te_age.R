# Edad del elemento transponible (TE) que se superpone con cada promotor de
# la library, a partir de la divergencia de RepeatMasker (perc_div =
# milliDiv/10: % de sustituciones de la copia respecto del consenso de su
# subfamilia; mas divergencia = copia mas antigua).
#
# Misma regla de seleccion que build_prom_features.R: TE = todo lo que no
# sea Low_complexity/Simple_repeat, superposicion > 10 pb con el fragmento
# de 252 pb. Por promotor se resume el TE dominante (mayor superposicion)
# y el promedio de perc_div ponderado por superposicion.
#
# Intermedio separado (no se toca prom_df.tsv): lee el RepeatMasker crudo
# (514MB, ver heavy_data_paths.R) una sola vez, igual que
# build_remap_tf_hits.R.
#
# Output: data/processed/te_age.tsv (seq_id, TE_name, TE_class,
# TE_superclass, TE_perc_div_dom, TE_overlap_dom, TE_perc_div_wmean,
# TE_n, TE_subfam_med_div, TE_subfam_n). Run from the TesisDoc repo root.

library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
source("R/00_prom_features/heavy_data_paths.R")

lib <- rtracklayer::import.bed("data/raw/library.bed")
values(lib) <- tibble(seq_id = lib$name)

rm_raw <- data.table::fread(path_repeatmasker, sep = " ", fill = TRUE, select = c("perc_div", "seqname", "start", "end", "matching_repeat", "repeat_class"))
rm_raw <- rm_raw[!repeat_class %in% c("Low_complexity", "Simple_repeat")]
# Edad por subfamilia: mediana genomica de perc_div de cada subfamilia
# (promedia el ruido de cada copia; no corrige el sesgo CpG pero es
# consistente entre copias de la misma subfamilia).
subfam_med <- rm_raw[, .(TE_subfam_med_div = median(perc_div), TE_subfam_n = .N), by = matching_repeat]
rm_raw <- rm_raw[seqname %in% seqlevels(lib)]
rm_gr <- GRanges(as.data.frame(rm_raw))

ov <- plyranges::join_overlap_intersect(lib, rm_gr) %>%
  subset(width > 10) %>%
  as_tibble() %>%
  transmute(
    seq_id, overlap = width, perc_div,
    TE_name = matching_repeat, TE_class = repeat_class,
    TE_superclass = str_remove(repeat_class, "/.+$")
  )

te_age <- ov %>%
  group_by(seq_id) %>%
  summarise(
    dom = which.max(overlap),
    TE_name = TE_name[dom], TE_class = TE_class[dom], TE_superclass = TE_superclass[dom],
    TE_perc_div_dom = perc_div[dom], TE_overlap_dom = overlap[dom],
    TE_perc_div_wmean = sum(perc_div * overlap) / sum(overlap),
    TE_n = n(),
    .groups = "drop"
  ) %>%
  select(-dom) %>%
  left_join(as_tibble(subfam_med), by = c("TE_name" = "matching_repeat"))

write_tsv(te_age, "data/processed/te_age.tsv")
message("Wrote data/processed/te_age.tsv with ", nrow(te_age), " promoters")
print(summary(te_age$TE_perc_div_dom))
print(count(te_age, TE_superclass, sort = TRUE))
