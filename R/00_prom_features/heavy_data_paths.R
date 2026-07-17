# Paths to large external/intermediate files NOT copied into TesisDoc/data
# (too big to duplicate on disk - see data/external/download_scripts/ for
# how each one is produced). Scripts source() this file and read from
# transcriptional_library for now. Once a file is actually copied locally,
# flip its path to the commented TesisDoc-native one below.
#
# Assumes TesisDoc and transcriptional_library are sibling folders
# (both directly under .../Labo/).

translib <- "../transcriptional_library"

# RepeatMasker hg38 annotation (514MB).
# Source: data/external/download_scripts/download_data.sh ("Repeated Elements" section)
# path_repeatmasker <- "data/external/repeatMasker_hg38.txt"
path_repeatmasker <- file.path(translib, "External_data/repeatMasker_hg38.txt")

# ENCODE DNase-seq bigWig (630MB). ENCFF165GHP.
# Source: ENCODE (https://www.encodeproject.org/files/ENCFF165GHP/) - not
# listed in download_data.sh, downloaded separately; confirm cell line with author.
# path_dnase <- "data/external/ENCFF165GHP_DNase.bigWig"
path_dnase <- file.path(translib, "External_data/ENCFF165GHP_DNase.bigWig")

# phyloP100way score, per-bp summary over the 1bp-resolution library (280MB).
# Source: data/external/download_scripts/phyloP_lib.sh (multiBigwigSummary
# over hg38.phyloP100way.bw, itself downloaded in download_data.sh)
# path_phylo_1bp <- "data/external/hg38.phyloP100way_summary_lib_1bp.bed"
path_phylo_1bp <- file.path(translib, "External_data/hg38.phyloP100way_summary_lib_1bp.bed")
