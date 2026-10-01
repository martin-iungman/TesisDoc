# discordancia_reportero_endogeno / 2_dif_signed_bland_altman /
# element_core_promoter_activity_wilcoxon_en
#
# EXPLORATORIO - reconstruccion en ingles del analisis Wilcoxon simple
# (efecto de cada elemento sobre la actividad CRUDA del reportero, sin
# el modelo de discordancia/avg_rank) que existia en una version
# anterior de element_core_promoter.R (misma carpeta) y que ya no esta
# en ese archivo - fue reemplazado por el enfoque LMM sobre dif_signed.
# El .jpg viejo (element_core_promoter_activity_wilcoxon.jpg) seguia en
# disco pero el codigo que lo genero se habia perdido; reconstruido a
# pedido del autor, restringido a los 4 elementos "downstream" del TSS
# (>+15) segun la nomenclatura de ElemeNT: MTE (~+18 a +27), bridge
# (~+17 a +22), DPE (~+28 a +32) y PB/Pause Button (~+30 a +40).
#
# Usa wilcox_effect_summary() (R/functions/plot_helpers.R, compartido
# con R7/summary_features.R y R5.3/promalt_summary.R) en vez de un
# calculo Wilcoxon propio.
#
# Fuente de los elementos: herramienta ElemeNT (Sloutskin et al.) sobre
# EPDnew, External_data/CORE_EPDnew_Aug2023.xlsx hoja "human"
# (path_element_epd en heavy_data_paths.R) - cada celda es "no" o una
# lista "posicion,score;..." de sitios encontrados, colapsada a
# presencia/ausencia (T si hay al menos un sitio).
#
# Requiere: data/processed/activity_stats_highconf.tsv,
# data/processed/prom_df.tsv y path_element_epd (heavy_data_paths.R -
# EPDnew, no copiado a TesisDoc). Run from the TesisDoc repo root.

library(tidyverse)
library(readxl)
library(coin)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

downstream_elements <- c("MTE", "bridge", "DPE", "PB")

element <- read_excel(path_element_epd, sheet = "human") %>%
  select(name, all_of(downstream_elements)) %>%
  mutate(across(all_of(downstream_elements), ~ .x != "no"))

tidy_data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(element, by = "name") %>%
  select(rep, mean, all_of(downstream_elements)) %>%
  mutate(
    rep = as.factor(rep),
    across(all_of(downstream_elements), ~ factor(.x, levels = c("TRUE", "FALSE")))
  )

wilcox_df <- wilcox_effect_summary(tidy_data)

p <- wilcox_df %>%
  pivot_wider(names_from = val, values_from = estimate) %>%
  mutate(
    sig = ifelse(pval_corr < 0.05, "FDR < 0.05", "ns"),
    feature = fct_reorder(feature, estimate)
  ) %>%
  ggplot(aes(x = estimate, y = feature, group = fct_inorder(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = sig, alpha = rep)) +
  scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.6)) +
  geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1) +
  scale_fill_manual(values = c("FDR < 0.05" = "#2980b9", "ns" = "grey70")) +
  labs(
    title = "MTE, DPE, bridge and PB: effect on reporter activity",
    x = "Effect on activity (Wilcoxon)", y = NULL, fill = NULL, alpha = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "element_core_promoter_activity_wilcoxon_en.jpg"), p, width = 9, height = 5, units = "in")

message("Figure saved to ", out_dir)
