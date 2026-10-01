# Worked example, same design as tata_discordance_example.R /
# repressive_mechanisms_example.R (same folder): visualize discordance
# for HIC1 and ZBTB14 (the two annotated repressors in the sobreestima
# set from the fully HEK293-matched screen, session 2026-09-11).
# HEK293 endogenous definition only (hek_tpm, HEK293-active promoters) -
# note both FLIP to subestima against max_tpm (any FANTOM5 sample), same
# pattern as TATA-box, not shown here per author request.
#
# Uses ReMap ChIP from any cell line (remap_tf_hits.tsv) as the predictor
# (not the HEK293/HEK293T-restricted ChIP).
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

markers <- remap_hits %>%
  filter(TF %in% c("HIC1", "ZBTB14")) %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

base <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data_hek <- base %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  group_by(rep) %>%
  mutate(
    rr = dense_rank(mean) / n(),
    re = row_number(hek_tpm) / n(),
    avg_rank = (rr + re) / 2,
    dif_signed = rr - re
  ) %>%
  ungroup()

both <- data_hek %>%
  select(seq_id, rep, avg_rank, dif_signed) %>%
  left_join(markers, by = "seq_id") %>%
  mutate(across(c(HIC1, ZBTB14), ~ replace_na(.x, FALSE))) %>%
  pivot_longer(c(HIC1, ZBTB14), names_to = "marker", values_to = "bound") %>%
  mutate(bound = ifelse(bound, "Present", "Absent"))

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
  facet_wrap(~marker) +
  scale_color_manual(values = c("Present" = "#c0392b", "Absent" = "grey50")) +
  scale_fill_manual(values = c("Present" = "#c0392b", "Absent" = "grey50")) +
  labs(x = "Activity bin (avg_rank, 100 bins)", y = "Mean disagreement (reporter - endogenous, vs. HEK293)", color = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

p_bottom <- both %>%
  ggplot(aes(bound, dif_signed, fill = bound)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  facet_wrap(~marker) +
  scale_fill_manual(values = c("Present" = "#c0392b", "Absent" = "grey70")) +
  labs(x = NULL, y = "Disagreement (reporter - endogenous, vs. HEK293)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p <- ggarrange(p_top, p_bottom, ncol = 1, heights = c(1, 1.1))
ggsave(file.path(out_dir, "hic1_zbtb14_example.jpg"), p, width = 9, height = 8, units = "in")

message("Saved to ", out_dir)
