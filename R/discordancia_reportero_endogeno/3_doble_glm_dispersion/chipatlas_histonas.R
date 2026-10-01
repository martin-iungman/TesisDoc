# discordancia_reportero_endogeno / 3_doble_glm_dispersion / histonas
#
# Analogo a remap_tf.R (misma carpeta) pero con la union de cada marca
# de histona (ChIP-Atlas) como feature. Doble-GLM (glmmTMB con
# dispformula) que ajusta simultaneamente la MEDIA y la DISPERSION de
# dif_signed. Sin GSEA (a diferencia del panel de TFs) - mismo criterio
# que R3.4/1_original_sin_corregir: no hay una ontologia funcional
# sensata para marcas de histona.
#
# Requiere: data/external/allPeaks_chipatlas_counted.tsv, data/processed/
# activity_stats_highconf.tsv, data/processed/prom_df.tsv y
# data/processed/fantom_endo_activity_summary.tsv (EXCEPCION, ver
# R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.

library(tidyverse)
library(glmmTMB)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

out_dir <- "figures/discordancia_reportero_endogeno/3_doble_glm_dispersion"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance()

chipatlas <- read_tsv("data/external/allPeaks_chipatlas_counted.tsv", col_names = c("seq_id", "feature", "n_samples"), show_col_types = FALSE)
explist <- read_tsv(path_chipatlas_explist, col_names = FALSE, show_col_types = FALSE) %>%
  select(X3, X4) %>%
  distinct(X4, .keep_all = TRUE)
chipatlas <- chipatlas %>%
  left_join(explist, by = c("feature" = "X4")) %>%
  rename(group = X3)

histones <- chipatlas %>% filter(group == "Histone")

binary_hist <- histones %>%
  distinct(seq_id, feature) %>%
  mutate(value = TRUE, feature = paste0("Hist_", feature)) %>%
  pivot_wider(names_from = feature, values_from = value, values_fill = FALSE)

hist_data <- data %>%
  select(seq_id, rep, avg_rank, dif_signed) %>%
  left_join(binary_hist, by = "seq_id") %>%
  mutate(across(starts_with("Hist_"), ~ replace_na(.x, FALSE)))

nHist <- hist_data %>%
  group_by(rep) %>%
  summarise(across(starts_with("Hist_"), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(starts_with("Hist_")) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  distinct(name)

hist_data <- hist_data %>% mutate(rep = as.factor(rep))

fit_double_glm <- function(var) {
  df <- hist_data
  df$feat_test <- df[[var]]
  fit <- tryCatch(
    glmmTMB(dif_signed ~ avg_rank + rep + feat_test, dispformula = ~ avg_rank + feat_test, data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(tibble(variable = str_remove(var, "^Hist_"), mean_estimate = NA_real_, mean_pval = NA_real_, disp_estimate = NA_real_, disp_pval = NA_real_))
  }
  s <- summary(fit)
  cc_mean <- s$coefficients$cond
  cc_disp <- s$coefficients$disp
  tibble(
    variable = str_remove(var, "^Hist_"),
    mean_estimate = cc_mean["feat_testTRUE", "Estimate"], mean_pval = cc_mean["feat_testTRUE", "Pr(>|z|)"],
    disp_estimate = cc_disp["feat_testTRUE", "Estimate"], disp_pval = cc_disp["feat_testTRUE", "Pr(>|z|)"]
  )
}

message("Ajustando doble-GLM para ", length(nHist$name), " marcas de histona...")
resultados <- map(nHist$name, fit_double_glm) %>%
  list_rbind() %>%
  mutate(mean_padj = p.adjust(mean_pval, method = "BH"), disp_padj = p.adjust(disp_pval, method = "BH"))
write_tsv(resultados, file.path(out_dir, "chipatlas_histonas_doble_glm.tsv"))

sig_colors <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")
add_sig <- function(df, padj_col) {
  df %>% mutate(sig = case_when(
    .data[[padj_col]] < 0.001 ~ "FDR < 0.001", .data[[padj_col]] < 0.01 ~ "FDR < 0.01",
    .data[[padj_col]] < 0.05 ~ "FDR < 0.05", TRUE ~ "ns"
  ) %>% factor(levels = names(sig_colors)))
}

p_mean <- resultados %>% add_sig("mean_padj") %>% mutate(variable = fct_reorder(variable, mean_estimate)) %>%
  ggplot(aes(x = mean_estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = mean_estimate, xmax = mean_estimate), height = 0) +
  geom_point(size = 3) +
  scale_color_manual(values = sig_colors) +
  labs(title = "Doble-GLM: efecto de marcas de histona sobre la MEDIA de la discordancia", x = "Estimado (β) — media", y = NULL, color = NULL) +
  theme_bw(base_size = 14) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "chipatlas_histonas_media.jpg"), p_mean, width = 9, height = 6.75, units = "in")

p_disp <- resultados %>% add_sig("disp_padj") %>% mutate(variable = fct_reorder(variable, disp_estimate)) %>%
  ggplot(aes(x = disp_estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  scale_color_manual(values = sig_colors) +
  labs(title = "Doble-GLM: efecto de marcas de histona sobre la DISPERSIÓN de la discordancia", x = "Estimado (β, escala log) — dispersión", y = NULL, color = NULL) +
  theme_bw(base_size = 14) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "chipatlas_histonas_dispersion.jpg"), p_disp, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
