# Resolves a figure's current numbered output folder from
# docs/mapping_figuras.csv, so scripts never hardcode a figure number
# that can change when the thesis is reorganized (only the slug is
# hardcoded in the calling script).

# Returns "M10", "R8", etc. for an unambiguous numero_actual, or
# "sinnum" when the CSV marks it unresolved (contains "?" or "/").
fig_number_prefix <- function(slug, mapping_path = "docs/mapping_figuras.csv") {
  mapping <- readr::read_csv(mapping_path, show_col_types = FALSE)
  row <- mapping[mapping$slug == slug, ]
  if (nrow(row) == 0) stop("slug '", slug, "' not found in ", mapping_path)
  m <- regmatches(row$numero_actual, regexpr("^Fig\\.\\s*[MR][0-9]+$", row$numero_actual))
  if (length(m) == 0) return("sinnum")
  sub("^Fig\\.\\s*", "", m)
}

# Returns (and by default creates) figures/<numero>_<slug>/ for the
# given slug. Scripts should always write their output through this,
# e.g. ggsave(file.path(fig_dir("gating_facs"), "panel_a.jpg"), ...).
fig_dir <- function(slug, create = TRUE, mapping_path = "docs/mapping_figuras.csv") {
  prefix <- fig_number_prefix(slug, mapping_path)
  dir <- file.path("figures", paste0(prefix, "_", slug))
  if (create) dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}
