# EXCEPTION to the "always regenerate from raw data" rule for TesisDoc -
# these read directly from transcriptional_library/Analysis/Tables (a
# precomputed intermediate of the OLD pipeline, not raw/reference data).
# Explicitly authorized by the author (2026-07-20) for R7
# (summary_features_secuencia) because no reasonably-reproducible raw
# recomputation was found within scope. PENDING: replace each of these
# with a from-raw-data build once the original generating pipeline is
# located/reconstructed.
#
# Assumes TesisDoc and transcriptional_library are sibling folders.

translib_analysis <- "../transcriptional_library/Analysis/Tables"

# Cis-regulatory module (CRM) counts: non-redundant ReMap ChIP-seq TF
# peaks per promoter. The raw source (External_data/remap2022_nr_macs2_
# hg38_v1_0.bed, 5.2GB, ~one row per individual TF:celltype peak) would
# need a bedtools-style merge + per-region TF aggregation we don't have
# reconstructed. PENDING: find/rebuild that aggregation script.
path_library_remap_crm <- file.path(translib_analysis, "library_remap_CRM.bed")

# FANTOM5 per-sample CAGE activity (136MB), used for tissue-specificity
# classification (sample_specificity_class) and the Gini index
# (sample_specificity_gini). PENDING: locate/rebuild the raw FANTOM5 ->
# per-sample-per-promoter counts aggregation.
path_sample_cage_activity <- file.path(translib_analysis, "sample_CAGE_activity.tsv")

# CAGEr-derived promoter shape (interquantile width) per promoter (1.8MB).
# PENDING: locate/rebuild the raw CAGEr shape-calling pipeline.
path_shape_merged <- file.path(translib_analysis, "shape_merged.tsv")

# HEK293 endogenous CAGE activity (3MB). PENDING: locate/rebuild the raw
# HEK293 CAGE processing (CAGEr) pipeline.
path_hek_cager <- file.path(translib_analysis, "HEK_CAGEr_df.tsv")
