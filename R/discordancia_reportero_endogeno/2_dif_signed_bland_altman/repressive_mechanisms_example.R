# Worked example, same idea as tata_discordance_example.R (same folder):
# visualize discordance (dif_signed) for specific features rather than
# just their LMM point estimate, comparing two distinct repressive
# mechanisms the reviewer's hypothesis B could operate through:
#   - Polycomb (PRC1): RNF2
#   - KRAB-ZNF heterochromatin: TRIM28 (KAP1, the universal KRAB-ZNF
#     corepressor) and H3K9me3 (the mark it deposits, via SETDB1)
# under both endogenous definitions (max_tpm, hek_tpm HEK293-active only).
# All three show SUBESTIMA under max_tpm (session 2026-09-11 remap_tf.R /
# chipatlas_histonas.R runs) - i.e. none confirm hypothesis B as stated;
# this checks whether that holds under hek_tpm too, same way TATA-box's
# sign flipped.
#
# Exploratory - not part of the numbered pipeline.

library(tidyverse)
library(ggpubr)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

# --- Markers: RNF2/TRIM28 (ReMap) + H3K9me3 (ChIP-Atlas) -------------------

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
tf_binary <- remap_hits %>%
  filter(TF %in% c("RNF2", "TRIM28")) %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

chipatlas <- read_tsv("data/external/allPeaks_chipatlas_counted.tsv", col_names = c("seq_id", "feature", "n_samples"), show_col_types = FALSE)
explist <- read_tsv(path_chipatlas_explist, col_names = FALSE, show_col_types = FALSE) %>%
  select(X3, X4) %>%
  distinct(X4, .keep_all = TRUE)
h3k9me3_binary <- chipatlas %>%
  left_join(explist, by = c("feature" = "X4")) %>%
  rename(group = X3) %>%
  filter(group == "Histone", feature == "H3K9me3") %>%
  distinct(seq_id) %>%
  mutate(H3K9me3 = TRUE)

markers <- tf_binary %>%
  full_join(h3k9me3_binary, by = "seq_id") %>%
  mutate(across(c(RNF2, TRIM28, H3K9me3), ~ replace_na(.x, FALSE)))

# --- Base discordance data, both endogenous definitions --------------------

base <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data_max <- base %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance() %>%
  mutate(endo_def = "vs. max activity (any FANTOM5 sample)")

data_hek <- base %>%
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
  data_max %>% select(seq_id, rep, avg_rank, dif_signed, endo_def),
  data_hek %>% select(seq_id, rep, avg_rank, dif_signed, endo_def)
) %>%
  left_join(markers, by = "seq_id") %>%
  mutate(across(c(RNF2, TRIM28, H3K9me3), ~ replace_na(.x, FALSE))) %>%
  pivot_longer(c(RNF2, TRIM28, H3K9me3), names_to = "marker", values_to = "bound") %>%
  mutate(
    marker = factor(marker, levels = c("RNF2", "TRIM28", "H3K9me3"),
                     labels = c("RNF2 (Polycomb/PRC1)", "TRIM28 (KRAB-ZNF)", "H3K9me3")),
    bound = ifelse(bound, "Present", "Absent"),
    endo_def = factor(endo_def, levels = c("vs. max activity (any FANTOM5 sample)", "vs. HEK293 activity (active only)"))
  )

# --- Top: binned dif_signed vs avg_rank ------------------------------------

binned <- both %>%
  group_by(marker, endo_def, rep) %>%
  mutate(bin = ntile(avg_rank, 100)) %>%
  group_by(marker, endo_def, bin, bound) %>%
  summarise(mean_dif = mean(dif_signed), se_dif = sd(dif_signed) / sqrt(n()), n = n(), .groups = "drop")

p_top <- binned %>%
  ggplot(aes(bin, mean_dif, color = bound, fill = bound)) +
  geom_ribbon(aes(ymin = mean_dif - 1.96 * se_dif, ymax = mean_dif + 1.96 * se_dif), alpha = 0.2, color = NA) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_grid(marker ~ endo_def) +
  scale_color_manual(values = c("Present" = "#c0392b", "Absent" = "grey50")) +
  scale_fill_manual(values = c("Present" = "#c0392b", "Absent" = "grey50")) +
  labs(x = "Activity bin (avg_rank, 100 bins)", y = "Mean disagreement (reporter - endogenous)", color = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# --- Bottom: violin ----------------------------------------------------

p_bottom <- both %>%
  ggplot(aes(bound, dif_signed, fill = bound)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  facet_grid(marker ~ endo_def) +
  scale_fill_manual(values = c("Present" = "#c0392b", "Absent" = "grey70")) +
  labs(x = NULL, y = "Disagreement (reporter - endogenous)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p <- ggarrange(p_top, p_bottom, ncol = 1, heights = c(1, 1.1))
ggsave(file.path(out_dir, "repressive_mechanisms_example.jpg"), p, width = 10, height = 13, units = "in")

message("Saved to ", out_dir)
