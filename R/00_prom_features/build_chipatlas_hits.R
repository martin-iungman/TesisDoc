# Shared intermediate: per-(seq_id, feature) ChIP-Atlas peak counts
# (n_samples), aggregated once via streaming awk over the raw peaks
# file (path_chipatlas_peaks_raw, 5.5GB, 157M rows - one row per raw
# peak) instead of loading it into R.
#
# Each raw line is "<n> <seq_id>\t<feature>". The leading "<n> " is NOT
# the sample count despite looking like one - verified against the old
# transcriptional_library/Analysis/Tables/allPeaks_light.hg38.50_lib_final.tsv
# (e.g. ALYREF-WT/"RNA polymerase II" has n_samples=838 there; in the
# raw file every one of its 838 rows has "1 " as the leading number,
# while counting raw LINES per (seq_id, feature) gives exactly 838).
# ported from transcriptional_library/Analysis/scripts/final_github.R
# ("ChIP-Atlas" section) and Histone_chipatlas.qmd, both of which read
# this file with read_tsv(col_names=c("seq_id","feature")) and then
# count(seq_id, feature) - that silently produces n_samples=1 for every
# row, because the un-stripped leading number becomes part of seq_id
# and blocks the aggregation from grouping correctly. Fixed here by
# stripping the leading number before counting.
#
# Run from the TesisDoc repo root. Output:
# data/external/allPeaks_chipatlas_counted.tsv (seq_id, feature, n_samples)
# Tarda ~10-15 minutos (157M lineas).

source("R/00_prom_features/heavy_data_paths.R")

# No header row - read with read_tsv(col_names = c("seq_id", "feature",
# "n_samples")) downstream (R/22_chipatlas_histonas/chipatlas_histonas.R).
awk_script <- '{
  seqid = $1
  sub(/^[ \t]*[0-9]+ /, "", seqid)
  count[seqid "\t" $2]++
}
END {
  for (k in count) print k "\t" count[k]
}'

out_path <- "data/external/allPeaks_chipatlas_counted.tsv"
dir.create("data/external", showWarnings = FALSE, recursive = TRUE)
status <- system2("awk", c("-F", "'\\t'", shQuote(awk_script), shQuote(path_chipatlas_peaks_raw)), stdout = out_path)
if (status != 0) stop("awk aggregation failed with status ", status)
message("Wrote ", out_path)
