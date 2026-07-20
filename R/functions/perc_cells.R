# Relative sampling-effort weights per sorting gate (1-6), from the
# % of cells collected in each gate during FACS sorting. Shared between
# R/01_activity_stats (normalization) and R/05_correccion_sesgo (the
# sample_effort plot), so both use the exact same numbers.
perc_cells <- tibble::tibble(sample = factor(1:6), perc = c(14.59, 8.87, 6.04, 3.09, 1.41, 1.03))
perc_cells$perc <- perc_cells$perc / sum(perc_cells$perc)
perc_cells$rel_factor <- perc_cells$perc / min(perc_cells$perc)
