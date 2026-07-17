# Renames figures/*_<slug>/ folders to match the current numero_actual
# in docs/mapping_figuras.csv. Run this after editing numero_actual for
# any row (e.g. after re-checking docs/Fig MyM.pptx for a renumbering -
# see docs/rutina_resync_figuras.md).
#
# Does NOT touch R/<numero>_<slug>/ script folders - those are numbered
# by pipeline order, not figure number, and never change.
# Does NOT decide whether the numbering changed - that requires actually
# reading the thesis content; this script only makes the folders match
# whatever is already in the CSV.

library(tidyverse)
source("R/functions/fig_paths.R")

mapping <- read_csv("docs/mapping_figuras.csv", show_col_types = FALSE)
existing <- list.dirs("figures", recursive = FALSE, full.names = FALSE)

for (slug in mapping$slug) {
  current <- existing[str_detect(existing, paste0("^([A-Za-z0-9]+_)?", slug, "$"))]
  if (length(current) == 0) next
  if (length(current) > 1) {
    warning("Multiple folders match slug '", slug, "': ", paste(current, collapse = ", "), " - skipping")
    next
  }
  target <- paste0(fig_number_prefix(slug), "_", slug)
  if (current != target) {
    message(current, " -> ", target)
    file.rename(file.path("figures", current), file.path("figures", target))
  }
}
