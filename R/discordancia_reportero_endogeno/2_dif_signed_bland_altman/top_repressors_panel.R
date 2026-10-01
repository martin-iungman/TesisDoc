# Panel: discordance (binned line + violin, same design as
# tata_discordance_example.R / repressive_mechanisms_example.R) for the
# 6 most interesting, well-characterized repressors among the 84 TFs
# concordant (significant, same direction) between the two max_tpm +
# HEK293-active-promoters screens (all ReMap ChIP vs HEK293/HEK293T-only
# ChIP; see remap_tf_maxtpm_hek_active.R / remap_tf_maxtpm_hek_chip_hek_
# active.R, same folder). All 6 show sobreestima (consistent with
# reviewer hypothesis B - repression lost when the promoter is taken out
# of its native context).
#
# Uses ReMap ChIP from any cell line + max_tpm endogenous, HEK293-active
# promoters only (the "all-ChIP" analysis) for the actual visualization -
# the effects also replicate with HEK293/HEK293T-restricted ChIP (see the
# concordance table, remap_tf_concordant_maxtpm_hek_active.tsv).
#
# Exploratory - not part of the numbered pipeline.

library(tidyverse)
library(ggpubr)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)
remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)

marcadores <- c("SETDB1", "NCOR1", "L3MBTL2", "TBL1X", "ZEB1", "ASXL1")

markers_bin <- remap_hits %>%
  filter(TF %in% marcadores) %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  add_rank_discordance()

both <- data %>%
  select(seq_id, rep, avg_rank, dif_signed) %>%
  left_join(markers_bin, by = "seq_id") %>%
  mutate(across(all_of(marcadores), ~ replace_na(.x, FALSE))) %>%
  pivot_longer(all_of(marcadores), names_to = "marker", values_to = "bound") %>%
  mutate(
    marker = factor(marker, levels = marcadores),
    bound = ifelse(bound, "Present", "Absent")
  )

binned <- both %>%
  group_by(marker, rep) %>%
  mutate(bin = ntile(avg_rank, 100)) %>%
  group_by(marker, bin, bound) %>%
  summarise(mean_dif = mean(dif_signed), se_dif = sd(dif_signed) / sqrt(n()), n = n(), .groups = "drop")

p_top <- binned %>%
  ggplot(aes(bin, mean_dif, color = bound, fill = bound)) +
  geom_ribbon(aes(ymin = mean_dif - 1.96 * se_dif, ymax = mean_dif + 1.96 * se_dif), alpha = 0.2, color = NA) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~marker, ncol = 3) +
  scale_color_manual(values = c("Present" = "#c0392b", "Absent" = "grey50")) +
  scale_fill_manual(values = c("Present" = "#c0392b", "Absent" = "grey50")) +
  labs(x = "Activity bin (avg_rank, 100 bins)", y = "Mean disagreement (reporter - endogenous)", color = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

p_bottom <- both %>%
  ggplot(aes(bound, dif_signed, fill = bound)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  facet_wrap(~marker, ncol = 3) +
  scale_fill_manual(values = c("Present" = "#c0392b", "Absent" = "grey70")) +
  labs(x = NULL, y = "Disagreement (reporter - endogenous)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p <- ggarrange(p_top, p_bottom, ncol = 1, heights = c(1, 1.1))
ggsave(file.path(out_dir, "top_repressors_panel.jpg"), p, width = 12, height = 10, units = "in")

message("Saved to ", out_dir)
