# Paper-ready version of top_repressors_panel.R (same folder): only the
# top (binned line) panel, no title, editable vector text via cairo_pdf.
# Same 6 concordant repressors (SETDB1, NCOR1, L3MBTL2, TBL1X, ZEB1,
# ASXL1), same data (max_tpm + all-ChIP + HEK293-active promoters).
#
# Output goes to ../transcriptional_library/Plots/Review/.

library(tidyverse)
library(ggpubr)
source("R/functions/plot_helpers.R")

plots_review_dir <- "../transcriptional_library/Plots/Review"

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

dir.create(plots_review_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plots_review_dir, "top_repressors_panel_top_english.jpg"), p_top, width = 12, height = 6, units = "in")
ggsave(file.path(plots_review_dir, "top_repressors_panel_top_english.pdf"), p_top, width = 12, height = 6, units = "in", device = cairo_pdf)

message("Figures saved to ", plots_review_dir)
