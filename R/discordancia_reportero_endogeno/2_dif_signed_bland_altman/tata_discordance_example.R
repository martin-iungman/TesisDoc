# Worked example: how to visualize the discordance (dif_signed) of a
# single feature (TATA-box), instead of just its LMM point estimate.
# Two panels: (top) binned dif_signed vs avg_rank, TATA present/absent as
# two lines with 95% CI ribbon - shows whether the gap holds across the
# activity range; (bottom) violin of dif_signed, TATA present vs absent.
# Shown for BOTH endogenous definitions (max_tpm, hek_tpm - HEK293-active
# only) side by side, since TATA is the feature whose sign flips between
# them (session 2026-09-11) - makes that dependence visually obvious
# instead of only a table entry.
#
# Exploratory - not part of the numbered pipeline.

library(tidyverse)
library(ggpubr)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data_max <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance() %>%
  mutate(endo_def = "vs. max activity (any FANTOM5 sample)")

data_hek <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  group_by(rep) %>%
  mutate(
    rr = dense_rank(mean) / n(),
    re = row_number(hek_tpm) / n(),
    avg_rank = (rr + re) / 2,
    dif_signed = rr - re
  ) %>%
  ungroup() %>%
  mutate(endo_def = "vs. HEK293 activity (active only)")

both <- bind_rows(
  data_max %>% select(seq_id, rep, avg_rank, dif_signed, `TATA-box` = TATA_EPD, endo_def),
  data_hek %>% select(seq_id, rep, avg_rank, dif_signed, `TATA-box` = TATA_EPD, endo_def)
) %>%
  mutate(
    `TATA-box` = ifelse(`TATA-box`, "TATA-box present", "TATA-box absent"),
    endo_def = factor(endo_def, levels = c("vs. max activity (any FANTOM5 sample)", "vs. HEK293 activity (active only)"))
  )

# --- Top: binned dif_signed vs avg_rank ------------------------------------

binned <- both %>%
  group_by(endo_def, rep) %>%
  mutate(bin = ntile(avg_rank, 100)) %>%
  group_by(endo_def, bin, `TATA-box`) %>%
  summarise(mean_dif = mean(dif_signed), se_dif = sd(dif_signed) / sqrt(n()), n = n(), .groups = "drop")

p_top <- binned %>%
  ggplot(aes(bin, mean_dif, color = `TATA-box`, fill = `TATA-box`)) +
  geom_ribbon(aes(ymin = mean_dif - 1.96 * se_dif, ymax = mean_dif + 1.96 * se_dif), alpha = 0.2, color = NA) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~endo_def) +
  scale_color_manual(values = c("TATA-box present" = "#c0392b", "TATA-box absent" = "grey50")) +
  scale_fill_manual(values = c("TATA-box present" = "#c0392b", "TATA-box absent" = "grey50")) +
  labs(x = "Activity bin (avg_rank, 100 bins)", y = "Mean disagreement\n(reporter - endogenous)", color = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# --- Bottom: violin ----------------------------------------------------

p_bottom <- both %>%
  ggplot(aes(`TATA-box`, dif_signed, fill = `TATA-box`)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  facet_wrap(~endo_def) +
  scale_fill_manual(values = c("TATA-box present" = "#c0392b", "TATA-box absent" = "grey70")) +
  labs(x = NULL, y = "Disagreement\n(reporter - endogenous)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p <- ggarrange(p_top, p_bottom, ncol = 1, heights = c(1, 1.1))
ggsave(file.path(out_dir, "tata_discordance_example.jpg"), p, width = 10, height = 9, units = "in")

message("Saved to ", out_dir)
