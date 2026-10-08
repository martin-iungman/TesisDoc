# EXPLORATORIO - actividad media segun union del TF (ReMap, cualquier
# linea celular), estratificada por presencia del motivo correspondiente:
# NFYA x CCAAT-box y SP1 x GC-box. Violines TF ausente/presente dentro de
# cada estrato del motivo, por replica, con el efecto de Wilcoxon
# (Hodges-Lehmann, presente - ausente) por estrato. Generaliza
# sp2_estratificado_ccaat.R (misma carpeta).
#
# Requiere: remap_tf_hits.tsv, activity_stats_highconf.tsv, prom_df.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

out_dir <- fig_dir("tf_nfya_sp_contexto_endogeno")

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

stratified_tf_plot <- function(tf, motif_col, motif_label, file_stem) {
  d <- data %>%
    mutate(
      TF = factor(ifelse(seq_id %in% remap_hits$seq_id[remap_hits$TF == tf], "Presencia", "Ausencia"), levels = c("Ausencia", "Presencia")),
      motivo = factor(ifelse(.data[[motif_col]], paste("Con", motif_label), paste("Sin", motif_label)),
        levels = paste(c("Sin", "Con"), motif_label)
      )
    )

  effect <- d %>%
    group_by(rep, motivo) %>%
    group_modify(~ {
      w <- wilcox.test(mean ~ TF, data = .x, conf.int = TRUE)
      tibble(n_tf = sum(.x$TF == "Presencia"), n = nrow(.x), estimate = -unname(w$estimate), p = w$p.value)
    }) %>%
    ungroup()
  print(effect)
  write_tsv(effect, file.path(out_dir, paste0(file_stem, ".tsv")))

  lab <- effect %>% mutate(label = sprintf("efecto = %.2f\nn(%s) = %d / %d", estimate, tf, n_tf, n))

  p <- ggplot(d, aes(TF, mean, fill = TF)) +
    geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
    stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), geom = "pointrange", size = 0.4) +
    geom_text(data = lab, aes(x = 1.5, y = 6.6, label = label), inherit.aes = FALSE, size = 4) +
    facet_grid(rep ~ motivo) +
    scale_fill_manual(values = c("Ausencia" = "#14AFB2", "Presencia" = "#216869")) +
    coord_cartesian(ylim = c(0.5, 7.2)) +
    labs(x = paste0("Unión de ", tf, " (ReMap)"), y = "Actividad media") +
    theme_pubclean(base_size = 16) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, paste0(file_stem, ".jpg")), p, width = 9, height = 7, units = "in")
}

stratified_tf_plot("NFYA", "CCAAT_EPD", "CCAAT-box", "nfya_estratificado_ccaat")
stratified_tf_plot("SP1", "GCbox_EPD", "GC-box", "sp1_estratificado_gcbox")

message("Figuras guardadas en ", out_dir)
