# especificidad_tisular_fantom (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Actividad endogena maxima (FANTOM5, todas las muestras) vs. actividad
# del reportero, separando promotores housekeeping (baja especificidad
# tisular) vs. tissue-specific (alta), segun tercil de sample_
# specificity_gini: tissue_specificity_boxplot.jpg (boxplots por deciles
# de actividad DENTRO de cada tercil) y tissue_specificity_residuals.jpg
# (residuos de un modelo lineal log10(tpm)~rank(actividad) por tercil,
# como medida de que tan bien la actividad del reportero predice la
# actividad endogena en cada grupo). Cada figura se repite restringida a
# promotores activos en HEK293 (hek_tpm > 0, sufijo _hek_active), con los
# MISMOS cortes de tercil que el conjunto completo (recalcularlos dentro
# del subconjunto bajaba el corte de tejido-especifico de gini 0.56 a 0.37
# y cambiaba la definicion del grupo). En
# cada boxplot se anota rho de Spearman por grupo con su IC 95% (bootstrap
# con la N propia de cada grupo) y la N. La
# figura completa (especificidad_tisular_fantom.jpg/.pdf, 2x2) junta las
# cuatro: todos los promotores a la izquierda, HEK293 a la derecha, con
# ambas replicas. Los boxplots tienen bigotes en los percentiles 10-90.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv y data/processed/fantom_endo_activity_summary.tsv (ver
# R/00_prom_features/build_fantom_endo_activity_summary.R - EXCEPCION,
# ver R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# ("Tissue specificity" section).

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "especificidad_tisular_fantom"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

set.seed(1)

add_terciles <- function(df) {
  df %>%
    group_by(rep) %>%
    mutate(specificity_tercile = cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) %>%
    ungroup()
}

# IC 95% de rho de Spearman por bootstrap (remuestreo con reposicion, con
# la N propia del grupo).
rho_boot_ci <- function(x, y, n_iter = 1000) {
  r <- replicate(n_iter, {
    i <- sample(length(x), replace = TRUE)
    cor(x[i], y[i], method = "spearman")
  })
  list(lo = unname(quantile(r, 0.025)), hi = unname(quantile(r, 0.975)))
}

data_all <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  mutate(max_tpm = replace_na(max_tpm, 0)) %>%
  filter(!is.na(sample_specificity_gini)) %>%
  add_terciles()

patron_colors <- c("Housekeeping" = "#1B8C8E", "Tejido-específico" = "#0D2C54")
label_patron <- function(tercile) ifelse(tercile == 1, "Housekeeping", "Tejido-específico")

save_both <- function(plot, name, width = 9, height = 6.75) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}

# Una faceta por replica solo si hay mas de una (la figura completa usa una).
facet_reps <- function(data, ...) {
  if (n_distinct(data$rep) > 1) facet_wrap(~rep, ...) else NULL
}

# --- Boxplot: actividad endogena maxima por bin de actividad, dentro de cada tercil ---
# Bigotes en p10-p90 (los de 1.5*IQR abarcaban 4 ordenes de magnitud y
# aplastaban las cajas).

plot_box <- function(data) {
  box_stats <- data %>%
    filter(specificity_tercile != 2) %>%
    group_by(rep, specificity_tercile) %>%
    mutate(mean_bin = row_number(mean) %>% cut_number(10) %>% as.numeric()) %>%
    group_by(rep, specificity_tercile, mean_bin) %>%
    summarise(
      ymin = quantile(max_tpm, 0.1), lower = quantile(max_tpm, 0.25), middle = median(max_tpm),
      upper = quantile(max_tpm, 0.75), ymax = quantile(max_tpm, 0.9),
      .groups = "drop"
    ) %>%
    mutate(patron = label_patron(specificity_tercile))

  # Spearman entre actividad del reportero y TPM maximo, por grupo (sobre
  # promotores individuales, no bins), arriba a la izquierda con su color.
  rho_df <- data %>%
    filter(specificity_tercile != 2) %>%
    group_by(rep, specificity_tercile) %>%
    summarise(
      rho = cor(mean, max_tpm, method = "spearman"),
      ci = list(rho_boot_ci(mean, max_tpm)),
      n = n(),
      .groups = "drop"
    ) %>%
    unnest_wider(ci) %>%
    mutate(
      patron = label_patron(specificity_tercile),
      label = sprintf("rho == \"%.2f [%.2f, %.2f]; n = %d\"", rho, lo, hi, n),
      vjust = ifelse(specificity_tercile == 1, 1.5, 3)
    )

  box_stats %>%
    ggplot(aes(factor(mean_bin), fill = patron)) +
    geom_boxplot(
      aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
      stat = "identity", alpha = 0.6, width = 0.75, position = position_dodge(0.8),
      key_glyph = draw_key_rect
    ) +
    geom_text(
      data = rho_df, aes(x = 0.6, y = Inf, label = label, colour = patron, vjust = vjust),
      hjust = 0, size = 4, parse = TRUE, inherit.aes = FALSE, show.legend = FALSE
    ) +
    scale_colour_manual(values = patron_colors, guide = "none") +
    facet_reps(data, scales = "free_x") +
    scale_y_log10() +
    labs(
      x = "Deciles de actividad del promotor reportero (dentro de cada grupo)",
      y = "TPM (valor máximo a través de FANTOM5)", fill = "Patrón de expresión"
    ) +
    scale_fill_manual(values = patron_colors) +
    theme_bw(base_size = 15) +
    theme(legend.position = "top")
}

# --- Residuos de actividad endogena ~ rango de actividad, por tercil -----

plot_residuals <- function(data) {
  lm_data <- data %>%
    filter(specificity_tercile != 2) %>%
    group_by(rep, specificity_tercile) %>%
    mutate(mean_rank = dense_rank(mean), log_tpm = log10(max_tpm + 0.1), grp = paste(rep, specificity_tercile, sep = "_")) %>%
    ungroup()

  lm_l <- lm_data %>% split(.$grp) %>% map(~ lm(log_tpm ~ mean_rank, data = .x))

  residuals_df <- lm_l %>%
    imap(~ tibble(residuals = residuals(.x), grp = .y) %>% separate(grp, into = c("rep", "specificity_tercile"), sep = "_")) %>%
    list_rbind() %>%
    mutate(patron = label_patron(specificity_tercile))

  residuals_df %>%
    ggplot(aes(residuals, fill = patron)) +
    geom_density(alpha = 0.6, key_glyph = draw_key_rect) +
    facet_reps(data) +
    scale_fill_manual(values = patron_colors) +
    labs(x = "Residuos (log10 TPM ~ rango de actividad)", y = "Densidad", fill = "Patrón de expresión") +
    theme_bw(base_size = 15) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
}

subsets <- list(
  todos = list(data = data_all, suffix = "", title = "Todos los promotores"),
  hek_active = list(
    data = data_all %>% filter(!is.na(hek_tpm), hek_tpm > 0), suffix = "_hek_active",
    title = "Promotores activos en HEK293"
  )
)

panels <- map(subsets, function(s) {
  df <- s$data
  save_both(plot_box(df), paste0("tissue_specificity_boxplot", s$suffix))
  save_both(plot_residuals(df), paste0("tissue_specificity_residuals", s$suffix))
  # Mismos limites en ambas columnas para poder comparar todos vs. HEK293
  list(
    box = plot_box(df) + ggtitle(s$title) + coord_cartesian(ylim = c(1, 5000)),
    res = plot_residuals(df) + coord_cartesian(xlim = c(-2, 3), ylim = c(0, 1.05))
  )
})

# Figura completa: columna izquierda todos los promotores,
# derecha activos en HEK293; boxplot arriba y sus residuos abajo.
figura <- panels$todos$box + panels$todos$res + panels$hek_active$box + panels$hek_active$res +
  plot_layout(ncol = 2, byrow = FALSE, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top", plot.tag = element_text(size = 24, face = "bold"))
save_both(figura, slug, 20, 12)

message("Figuras guardadas en ", out_dir)
