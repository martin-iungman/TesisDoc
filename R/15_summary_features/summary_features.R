# Dos figuras (ver docs/mapping_figuras.csv):
#   - summary_features_secuencia: panel de secuencia (summary_seq.jpg)
#   - summary_features_endogeno: panel de contexto endogeno
#     (summary_endo.jpg), figura propia
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

out_dir_seq <- fig_dir("summary_features_secuencia")
out_dir_endo <- fig_dir("summary_features_endogeno")

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

plot_summary <- function(df, features, y_label) {
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
    labs(fill = "Efecto", x = "Efecto sobre la actividad", y = y_label, alpha = "") +
    geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5) +
    # -0.35: el IC de "Tracto poli(dA:dT) proximal" llega a -0.30 (Rep 2)
    xlim(c(-0.35, 0.55)) +
    scale_fill_manual(values = c("Actividad alta" = "#D6741F", "Actividad baja" = "#7FB800")) +
    guides(alpha = "none")
}

save_jpg_pdf <- function(plot, dir, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}

panel_seq <- plot_summary(wilcox_df, vbles_split[[1]], "Elementos del promotor")
save_jpg_pdf(panel_seq, out_dir_seq, "summary_seq", 13.5, 6.75)

panel_endo <- plot_summary(wilcox_df, vbles_split[[2]], "Características del promotor")
save_jpg_pdf(panel_endo, out_dir_endo, "summary_endo", 13.5, 6.75)
# el panel A de summary_features_endogeno se guarda como objeto para que
# R/40_contexto_ajustado/contexto_ajustado.R arme la figura completa (A + B)
saveRDS(panel_endo, "data/processed/panel_summary_endo.rds")

message("Figuras guardadas en ", out_dir_seq, " y ", out_dir_endo)
