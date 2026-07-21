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

# FANTOM5 primary-cell raw CTSS files (6.7GB) - only used via list.files()
# to enumerate which samples exist, never read in full; safe to reference
# directly without copying.
path_fantom5_primary_cell_raw <- file.path(translib, "External_data/FANTOM5/hg38_primary_cell")

# FANTOM5 CAGE signal pre-filtered to library positions, per ontology term
# (2.9GB, 109 files). Source: transcriptional_library/External_data/
# filter_lib_fantom5.sh (not yet copied here - ask if it needs bringing in).
# path_fantom5_filt_cage <- "data/external/FANTOM5/Library_filt_CAGE"
path_fantom5_filt_cage <- file.path(translib, "External_data/FANTOM5/Library_filt_CAGE")

# ReMap non-redundant peaks (5.2GB raw), used for the CRM ("any TF peak
# overlap") feature - see build_prom_features.R.
# path_remap_nr <- "data/external/remap2022_nr_macs2_hg38_v1_0.bed"
path_remap_nr <- file.path(translib, "External_data/remap2022_nr_macs2_hg38_v1_0.bed")

# Flow cytometry validation, individual per-cell EGFP fluorescence tables
# (519MB, one .cells.csv per sample/promoter). Source: FlowJo export, no
# regeneration script (raw instrument data).
# path_citometry_stable_validation <- "data/external/Citometry_Stable_validation"
path_citometry_stable_validation <- file.path(translib, "Experimental_data/Citometry/Stable_validation/Sample Group - 1/Tables")
