# EXPLORATORIO - efecto de la union de SP2 (ReMap, cualquier linea) sobre
# la actividad media, estratificado por presencia de CCAAT-box. Violines
# SP2 ausente/presente dentro de cada estrato de CCAAT, por replica, con
# el efecto de Wilcoxon (Hodges-Lehmann) dentro de cada estrato.
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

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  mutate(
    SP2 = factor(ifelse(seq_id %in% remap_hits$seq_id[remap_hits$TF == "SP2"], "Presencia", "Ausencia"), levels = c("Ausencia", "Presencia")),
    CCAAT = factor(ifelse(CCAAT_EPD, "Con CCAAT-box", "Sin CCAAT-box"), levels = c("Sin CCAAT-box", "Con CCAAT-box"))
  )

# efecto (Hodges-Lehmann, SP2 presente - ausente) y n por estrato
effect <- data %>%
  group_by(rep, CCAAT) %>%
  group_modify(~ {
    w <- wilcox.test(mean ~ SP2, data = .x, conf.int = TRUE)
    tibble(
      n_sp2 = sum(.x$SP2 == "Presencia"), n = nrow(.x),
      estimate = -unname(w$estimate), p = w$p.value
    )
  }) %>%
  ungroup()
print(effect)
write_tsv(effect, file.path(out_dir, "sp2_estratificado_ccaat.tsv"))

lab <- effect %>%
  mutate(label = sprintf("efecto = %.2f\nn(SP2) = %d / %d", estimate, n_sp2, n))

p <- ggplot(data, aes(SP2, mean, fill = SP2)) +
  geom_violin(alpha = 0.7, draw_quantiles = 0.5) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), geom = "pointrange", size = 0.4) +
  geom_text(data = lab, aes(x = 1.5, y = 6.6, label = label), inherit.aes = FALSE, size = 4) +
  facet_grid(rep ~ CCAAT) +
  scale_fill_manual(values = c("Ausencia" = "#14AFB2", "Presencia" = "#216869")) +
  coord_cartesian(ylim = c(0.5, 7.2)) +
  labs(x = "Unión de SP2 (ReMap)", y = "Actividad media") +
  theme_pubclean(base_size = 16) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "sp2_estratificado_ccaat.jpg"), p, width = 9, height = 7, units = "in")

message("Figura guardada en ", out_dir)
