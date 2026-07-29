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

# Same underlying data as path_sample_cage_activity above (per-sample
# FANTOM5 CAGE counts per promoter, columns name/counts/sample/libsize),
# split into tissue (38MB) and primary_cell (83MB) instead of one
# combined file - used for R4.2/R4.3/R4.4 (transcriptional_library/
# Analysis/scripts/final_github.R "Tissue specificity" section: mean/max
# endogenous activity across samples, tissue-specificity boxplots).
# Extends the same authorized exception, not a new one. PENDING: same
# raw-data rebuild as path_sample_cage_activity.
path_tissue_cage_activity <- file.path(translib_analysis, "tissue_CAGE_activity.tsv")
path_primary_cell_cage_activity <- file.path(translib_analysis, "primary_cell_CAGE_activity.tsv")

# CAGEr-derived promoter shape (interquantile width) per promoter (1.8MB).
# PENDING: locate/rebuild the raw CAGEr shape-calling pipeline.
path_shape_merged <- file.path(translib_analysis, "shape_merged.tsv")

# HEK293 endogenous CAGE activity: RESOLVED 2026-07-21, now built from
# raw FANTOM5 data in build_prom_features.R (ported from
# transcriptional_library/Analysis/scripts/hek_cage.qmd via CAGEr) -
# no longer reads from Analysis/Tables.

# PUFFIN (Dudnyk et al. 2024) prediction scores for our own library
# sequences (teoTSS_score and related motif/effect columns), used by R4.5
# (puffin_selectivity_predscore). Produced by running the external
# puffin.py model + weights (not present in this repo) over ~24000
# promoters x 5 random upstream-flank replicates - see
# R/25_puffin/puffin_prediction_model.R for the ported (NOT RUN)
# pipeline, and transcriptional_library/Analysis/scripts/
# puffin_processing.qmd for the original. The per-sequence raw prediction
# CSVs that pipeline writes no longer exist on disk; only this aggregated
# summary survived. Authorized 2026-07-29 for R4.5 because rerunning the
# actual model is out of scope. PENDING: rerun puffin_prediction_model.R
# from raw library sequences if puffin.py/weights become available again.
path_puffin_pred <- file.path(translib_analysis, "Dudnyk_puffin_prediction_summ.tsv")
