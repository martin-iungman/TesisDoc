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

# ROC curve (por umbral de rango) de una feature booleana contra
# var_rank_sw (rango de varianza dentro de cada bin de actividad media -
# ver add_mean_sw_bins), para testear si una feature de secuencia predice
# ruido inusualmente alto/bajo. Ported from transcriptional_library/
# Analysis/scripts/noise_analysis.qmd.
roc_curve <- function(data, feature) {
  purrr::map(1:99, ~ data %>%
    count(gr = var_rank_sw > .x, across(all_of(feature))) %>%
    mutate(threshold = .x)) %>%
    purrr::list_rbind() %>%
    pivot_wider(names_from = c(gr, all_of(feature)), values_from = n) %>%
    mutate(TPR = TRUE_TRUE / (TRUE_TRUE + FALSE_TRUE), FPR = TRUE_FALSE / (FALSE_FALSE + TRUE_FALSE)) %>%
    arrange(FPR)
}

# Area bajo la curva (regla del trapecio) de un roc_curve().
noise_auc <- function(roc) {
  total <- 0
  for (i in seq_len(nrow(roc) - 1)) {
    deltax <- roc$FPR[i + 1] - roc$FPR[i]
    deltay <- roc$TPR[i + 1] - roc$TPR[i]
    total <- sum(total, deltax * roc$TPR[i] + deltax * deltay / 2, na.rm = TRUE)
  }
  total
}

# Scatter de ruido (var_rank_sw) vs. bin de actividad media, coloreado por
# una feature booleana, con smooth por nivel (paleta separada para puntos
# vs. curvas via ggnewscale, igual que el original).
noise_scatter <- function(data, feature, titulo, leyenda) {
  data %>%
    ggplot(aes(mean_sw, var_rank_sw, col = .data[[feature]])) +
    geom_point(size = 0.2, alpha = 0.5) +
    facet_wrap(~rep) +
    labs(x = "Bin de actividad media", y = "Rango de varianza (dentro del bin)", col = leyenda, title = titulo) +
    scale_color_manual(values = c("grey", thesis_clr), guide = guide_legend(override.aes = list(size = 3))) +
    ggnewscale::new_scale_color() +
    geom_smooth(aes(col = .data[[feature]]), se = FALSE, linewidth = 1.2, show.legend = FALSE) +
    scale_color_manual(values = c("#525252", "#136869")) +
    theme_pubr(base_size = 16)
}

# Curva ROC (una linea por replica, con AUC en el titulo) de una feature
# booleana prediciendo var_rank_sw alto.
noise_roc <- function(data, feature, titulo) {
  data %>%
    group_split(rep) %>%
    purrr::map(~ roc_curve(.x, feature) %>% mutate(rep = unique(.x$rep), AUC = noise_auc(.))) %>%
    purrr::list_rbind() %>%
    ggplot(aes(FPR, TPR, col = rep)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", col = "grey50") +
    geom_text(
      data = . %>% distinct(AUC, rep) %>% mutate(AUC = paste("AUC:", round(AUC, 3))) %>%
        bind_cols(tibble(TPR = c(0.9, 0.8), FPR = 0.02)),
      aes(label = AUC), hjust = "left", show.legend = FALSE
    ) +
    ggtitle(titulo) +
    scale_color_manual(values = c("#0D2C54", "#AD343E")) +
    theme_pubr(base_size = 16) +
    labs(col = "")
}
