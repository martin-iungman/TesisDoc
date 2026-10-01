# Panel: discordance (binned line + violin, same design as
# tata_discordance_example.R / repressive_mechanisms_example.R) for the
# 6 most interesting, well-characterized repressors among the 391
# sobreestima TFs in the "all ReMap ChIP peaks" screen only
# (remap_tf_maxtpm_hek_active.R, same folder: max_tpm endogenous,
# HEK293-active promoters, ChIP from any cell line) - NOT requiring
# concordance with the HEK293/HEK293T-ChIP-restricted screen (that
# smaller/underpowered analysis excludes many real repressors, like
# RNF2/EZH2/JARID2 earlier, purely for lack of HEK293-specific ChIP
# coverage, not lack of effect). Companion to top_repressors_panel.R
# (the 6-repressor panel restricted to TFs concordant in both screens).
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

marcadores <- c("BCL6", "MNT", "KDM5A", "MGA", "E2F8", "MBD2")

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
ggsave(file.path(out_dir, "top_repressors_panel_allpeaks.jpg"), p, width = 12, height = 10, units = "in")

message("Saved to ", out_dir)
