# summary_features_secuencia = Fig. R7 (ver docs/mapping_figuras.csv)
# Resumen del efecto de features de secuencia (seq) y de contexto
# endogeno/genomico (endo) sobre la actividad transcripcional media:
# para cada feature booleana, test de Wilcoxon (mean ~ feature, pareado
# por replica) con correccion BH: las que dan significativas se grafican
# como barras con el tamano de efecto (estimate) e IC 95%.
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y
# R/01_activity_stats - incluye las features de
# R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(coin)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")

slug <- "summary_features_secuencia"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

# --- Tabla de features booleanas (todas en espanol) -----------------------
# build_tidy_features()/feature_groups en R/functions/plot_helpers.R
# (compartido con R3.2, el mismo resumen pero para ruido en vez de actividad).

# Islas CpG definidas por la composicion del fragmento (CGI_frag, default
# de build_tidy_features()).
tidy_data <- build_tidy_features(data, keep = "mean")
vbles_split <- map(c("seq", "endo"), ~ feature_groups$feature[feature_groups$group == .x])

# --- Wilcoxon: efecto de cada feature sobre la actividad media ----------
# wilcox_effect_summary() en R/functions/plot_helpers.R (compartido con
# R5.3, el mismo patron pero para las categorias de promotor alternativo).

wilcox_df <- wilcox_effect_summary(tidy_data)

# saved for R/09_coocurrencia (coocurrencia_motivos), which needs the significant/consistent
# feature list without recomputing every Wilcoxon test.
write_tsv(wilcox_df, "data/processed/activity_summary.tsv")

# --- Graficos resumen: seq y endo ------------------------------------------

plot_summary <- function(df, features) {
  df %>%
    filter(feature %in% features, pval_corr < 0.05) %>%
    mutate(act = ifelse(sign(estimate) == 1, "Actividad alta", "Actividad baja")) %>%
    group_by(feature) %>%
    mutate(n = length(unique(act))) %>%
    filter(n == 1) %>%
    arrange(estimate, desc(rep)) %>%
    pivot_wider(names_from = val, values_from = estimate) %>%
    filter(pval_corr < 0.05) %>%
    ggplot(aes(x = estimate, y = fct_inorder(feature), group = fct_inorder(rep))) +
    geom_col(orientation = "y", position = "dodge", aes(fill = act, alpha = rep)) +
    scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7)) +
    theme_pubr() +
    theme(text = element_text(size = 20), legend.position = "top") +
    labs(fill = "Efecto", x = "Efecto sobre la actividad", y = "Features", alpha = "") +
    geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5) +
    xlim(c(-0.3, 0.55)) +
    scale_fill_manual(values = c("Actividad alta" = "#D6741F", "Actividad baja" = "#7FB800")) +
    guides(alpha = "none")
}

panel_seq <- plot_summary(wilcox_df, vbles_split[[1]])
ggsave(file.path(out_dir, "summary_seq.jpg"), panel_seq, width = 13.5, height = 6.75, units = "in")

panel_endo <- plot_summary(wilcox_df, vbles_split[[2]])
ggsave(file.path(out_dir, "summary_endo.jpg"), panel_endo, width = 13.5, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
