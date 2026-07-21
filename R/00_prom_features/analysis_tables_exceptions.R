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

# Cis-regulatory module (CRM) overlap: RESOLVED 2026-07-21, now built as a
# 0/1 "any ReMap TF peak overlap" flag directly from the raw
# External_data/remap2022_nr_macs2_hg38_v1_0.bed (5.2GB) in
# build_prom_features.R - no longer reads library_remap_CRM.bed. Note this
# is a simplification, not a port: it drops ReMap's CRM (merged-region,
# TF-count) aggregation since N_TF_CRM was only ever used downstream as
# N_TF_CRM == 0.

# FANTOM5 per-sample CAGE activity (136MB), used for tissue-specificity
# classification (sample_specificity_class) and the Gini index
# (sample_specificity_gini). PENDING: locate/rebuild the raw FANTOM5 ->
# per-sample-per-promoter counts aggregation.
path_sample_cage_activity <- file.path(translib_analysis, "sample_CAGE_activity.tsv")

# CAGEr-derived promoter shape (interquantile width) per promoter (1.8MB).
# PENDING: locate/rebuild the raw CAGEr shape-calling pipeline.
path_shape_merged <- file.path(translib_analysis, "shape_merged.tsv")

# HEK293 endogenous CAGE activity: RESOLVED 2026-07-21, now built from
# raw FANTOM5 data in build_prom_features.R (ported from
# transcriptional_library/Analysis/scripts/hek_cage.qmd via CAGEr) -
# no longer reads from Analysis/Tables.
