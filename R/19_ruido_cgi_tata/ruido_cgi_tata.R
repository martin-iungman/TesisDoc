# ruido_cgi y ruido_tata (ver docs/mapping_figuras.csv para el numero de
# figura vigente) - antes una sola figura (ruido_cgi_tata).
# Ruido transcripcional (rango de varianza dentro de cada bin de
# actividad media, var_rank_sw) segun presencia de isla CpG / TATA-box:
# scatter de ruido vs. bin de actividad + curva ROC (con AUC por
# replica) de la feature prediciendo ruido alto.
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y
# R/01_activity_stats). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "## Noise", bloques CGI y TATA) y noise_analysis.qmd
# (funciones roc_curve/auc, aca roc_curve/noise_auc en
# R/functions/plot_helpers.R).

library(tidyverse)
library(ggpubr)
library(ggnewscale)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

out_dir_cgi <- fig_dir("ruido_cgi")
out_dir_tata <- fig_dir("ruido_tata")

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data <- add_noise_rank(data)

save_jpg_pdf <- function(plot, dir, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}
tag_theme <- theme(plot.tag = element_text(size = 24, face = "bold"))

# --- CGI ---------------------------------------------------------------
# Islas CpG definidas por la composicion del fragmento (CGI_frag: G+C >= 50%
# y CpG o/e > 0.6), no por la anotacion genomica UCSC (CGI).

scatter_cgi <- noise_scatter(data, "CGI_frag", "Isla CpG", NULL)
roc_cgi <- noise_roc(data, "CGI_frag", "Curva ROC - Isla CpG")
save_jpg_pdf(scatter_cgi, out_dir_cgi, "noise_scatter_CGI", 9, 6.75)
save_jpg_pdf(roc_cgi, out_dir_cgi, "ROC_CGI", 7, 6.75)

# Figura completa: scatter (A) y curva ROC (B) lado a lado
figura_cgi <- scatter_cgi + roc_cgi + plot_layout(widths = c(9, 7)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_jpg_pdf(figura_cgi, out_dir_cgi, "ruido_cgi", 16, 6.75)

# --- TATA-box ------------------------------------------------------------

scatter_tata <- noise_scatter(data, "TATA_EPD", "TATA-box", NULL)
roc_tata <- noise_roc(data, "TATA_EPD", "Curva ROC - TATA-box")
roc_tata_expr <- noise_roc_by_expr(data, "TATA_EPD", "Curva ROC - TATA-box")
save_jpg_pdf(scatter_tata, out_dir_tata, "noise_scatter_TATA", 9, 6.75)
save_jpg_pdf(roc_tata, out_dir_tata, "ROC_TATA", 7, 6.75)
save_jpg_pdf(roc_tata_expr, out_dir_tata, "ROC_TATA_expr", 5, 5)

# Figura completa: scatter (A), curva ROC (B) y curva ROC por nivel de
# actividad (C) en una fila. Todo scatter de ruido vs. media va a 9 x 6.75 in
# dentro de la figura, igual que en ruido_cgi.
figura_tata <- scatter_tata + roc_tata + roc_tata_expr +
  plot_layout(widths = c(9, 7, 7)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_jpg_pdf(figura_tata, out_dir_tata, "ruido_tata", 23, 6.75)

message("Figuras guardadas en ", out_dir_cgi, " y ", out_dir_tata)
