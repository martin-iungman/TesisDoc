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

# Bins de 100 promotores segun rango de actividad media, por replica (se
# descarta el resto de la division para que los bins queden parejos). Usado
# para los scatters de union de TF/motivo vs. actividad (R/16_tf_contexto_
# endogeno, R/17_remap_activity_gsea).
add_mean_sw_bins <- function(df) {
  excess <- df %>% group_by(rep) %>% summarise(excess = n() %% 100, .groups = "drop")
  df %>%
    group_split(rep) %>%
    purrr::map2(excess$excess, ~ .x %>%
      mutate(mean_rank_sw = row_number(mean)) %>%
      filter(mean_rank_sw > .y) %>%
      mutate(mean_rank_sw = row_number(mean), mean_sw = ceiling(mean_rank_sw / 100)) %>%
      ungroup()) %>%
    purrr::list_rbind()
}

# Proporcion de promotores con `feature`==TRUE por bin de actividad
# (mean_sw, ver add_mean_sw_bins), con smooth loess por replica.
tf_scatter <- function(df, feature, titulo) {
  df %>%
    group_by(rep, mean_sw) %>%
    summarise(prop = sum(.data[[feature]]) / n(), .groups = "drop") %>%
    ggplot(aes(mean_sw, prop)) +
    geom_point(col = thesis_clr) +
    geom_smooth(col = "#216869") +
    facet_wrap(~rep) +
    ylim(0, 1) +
    labs(x = "Actividad media (bins de 100, ordenados por rango)", y = "Proporción de promotores", title = titulo) +
    theme_bw(base_size = 16)
}
