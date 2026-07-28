# ruido_cgi_tata (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
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
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "ruido_cgi_tata"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data <- add_mean_sw_bins(data) %>%
  group_by(rep, mean_sw) %>%
  mutate(var_rank_sw = row_number(var)) %>%
  ungroup()

# --- CGI ---------------------------------------------------------------

ggsave(
  file.path(out_dir, "noise_scatter_CGI.jpg"),
  noise_scatter(data, "CGI", "Ruido vs. actividad, por isla CpG", "Isla CpG"),
  width = 9, height = 6.75, units = "in"
)
ggsave(
  file.path(out_dir, "ROC_CGI.jpg"),
  noise_roc(data, "CGI", "Curva ROC - Isla CpG"),
  width = 7, height = 6.75, units = "in"
)

# --- TATA-box ------------------------------------------------------------

ggsave(
  file.path(out_dir, "noise_scatter_TATA.jpg"),
  noise_scatter(data, "TATA_EPD", "Ruido vs. actividad, por TATA-box", "TATA-box"),
  width = 9, height = 6.75, units = "in"
)
ggsave(
  file.path(out_dir, "ROC_TATA.jpg"),
  noise_roc(data, "TATA_EPD", "Curva ROC - TATA-box"),
  width = 7, height = 6.75, units = "in"
)
ggsave(
  file.path(out_dir, "ROC_TATA_expr.jpg"),
  noise_roc_by_expr(data, "TATA_EPD", "Curva ROC - TATA-box"),
  width = 5, height = 5, units = "in"
)

message("Figuras guardadas en ", out_dir)
