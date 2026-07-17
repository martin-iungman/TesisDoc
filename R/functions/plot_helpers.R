# Shared ggplot helpers/theme for thesis figures (Spanish labels, ggpubr style).

library(ggplot2)
library(ggpubr)

thesis_clr <- "#358AAA"
thesis_clr_bg <- "#dfedf6"
thesis_clr_label <- "#AADAD4"

# Barplot of counts against a total: a full-height background bar (the
# total), a foreground bar (the count) and a value label on top. Used for
# "N of total promoters have this feature" figures (EPD motifs, CGI, TSS
# motifs, conservation classes...).
count_bar_labeled <- function(data, x, n, total, base_size = 16,
                               label_nudge_frac = 0.07, coord_flip = FALSE) {
  x <- rlang::enquo(x)
  n <- rlang::enquo(n)
  total <- rlang::enquo(total)
  # `total` may be a scalar from the caller's env (e.g. n_promoters) or a
  # column in `data` (e.g. total_n); eval_tidy resolves either.
  total_val <- max(rlang::eval_tidy(total, data), na.rm = TRUE)
  p <- ggplot(data, aes(x = !!x, label = !!n)) +
    geom_col(aes(y = !!total), fill = thesis_clr_bg) +
    geom_col(aes(y = !!n), fill = thesis_clr) +
    geom_label(aes(y = !!n), nudge_y = total_val * label_nudge_frac, fill = thesis_clr_label, size = 4) +
    theme_pubclean(base_size = base_size)
  if (coord_flip) p <- p + coord_flip()
  p
}
