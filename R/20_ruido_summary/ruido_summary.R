# ruido_summary (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Resumen del efecto de cada feature booleana del promotor sobre el
# ruido transcripcional: AUC (con IC bootstrap) de la feature
# prediciendo var_rank_sw alto (rango de varianza dentro de cada bin de
# actividad media - ver R/19_ruido_cgi_tata), un solo panel combinado
# (sin separar seq/endo, a diferencia de R7).
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y
# R/01_activity_stats). Tarda unos minutos: ~32 features x 2 replicas x
# 2000 iteraciones de bootstrap para el IC de cada AUC. Run from the
# TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (bloque de resumen de la seccion de ruido, plot_auc2). build_tidy_
# features()/feature_groups (identicas a las de R/15_summary_features,
# el resumen de actividad) viven en R/functions/plot_helpers.R.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(pROC)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "ruido_summary"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data <- add_mean_sw_bins(data) %>%
  group_by(rep, mean_sw) %>%
  mutate(var_rank_sw = row_number(var)) %>%
  ungroup()

tidy_data <- build_tidy_features(data, keep = "var_rank_sw")

# --- AUC (con IC bootstrap) de cada feature prediciendo ruido alto -------

vbles <- tidy_data %>% select(-rep, -var_rank_sw) %>% names()
repname <- unique(tidy_data$rep)

auc_df <- map(repname, function(r) {
  df <- tidy_data %>% filter(rep == r)
  map(vbles, function(f) {
    roc_obj <- pROC::roc(df[[f]], df[["var_rank_sw"]], ci = TRUE, boot.n = 2000, direction = ">", quiet = TRUE)
    tibble(feature = f, rep = r, AUC = as.numeric(roc_obj$ci[2]), ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]))
  }) %>% list_rbind()
}) %>% list_rbind()

write_tsv(auc_df, file.path(out_dir, "ruido_summary_auc.tsv"))

# --- Grafico resumen --------------------------------------------------
# Features con IC que no cruza 0.5 Y mismo sentido (ruido alto/bajo) en
# ambas replicas (mismo criterio de consistencia que R7, no estaba en
# el plot_auc2 original).

panel <- auc_df %>%
  mutate(noise = ifelse(AUC > 0.5, "Ruido alto", "Ruido bajo")) %>%
  filter((ci2.5 > 0.5 & ci97.5 > 0.5) | (ci2.5 < 0.5 & ci97.5 < 0.5)) %>%
  group_by(feature) %>%
  mutate(n_dir = length(unique(noise))) %>%
  filter(n_dir == 1, n() == 2) %>%
  ungroup() %>%
  arrange(desc(rep), AUC) %>%
  mutate(feature = fct_inorder(feature)) %>%
  ggplot(aes(x = AUC - 0.5, y = feature, group = fct_inorder(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = noise, alpha = rep)) +
  geom_errorbarh(aes(xmax = ci2.5 - 0.5, xmin = ci97.5 - 0.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5) +
  scale_x_continuous(labels = function(x) x + 0.5) +
  scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7)) +
  scale_fill_manual(values = c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800")) +
  theme_pubr() +
  theme(text = element_text(size = 20), legend.position = "top") +
  labs(fill = "Efecto", x = "AUC (efecto sobre el ruido)", y = "Features", alpha = "")

ggsave(file.path(out_dir, "summary_noise.jpg"), panel, width = 13.5, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
