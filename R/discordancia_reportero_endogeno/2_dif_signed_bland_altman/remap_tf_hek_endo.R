# Variant of remap_tf.R (same folder) using HEK293 CAGE activity (hek_tpm)
# as the endogenous side instead of max_tpm across all FANTOM5 samples -
# same motivation/design as prom_df_features_hek_endo.R (session
# 2026-09-11). HEK293-active promoters only (hek_tpm > 0) - the "clean"
# variant (no floor/tie block at hek_tpm==0), matching what we've been
# using for the feature-level comparison. No GSEA/ORA here, just the LMM
# screen (add separately if needed).
#
# Does not touch remap_tf.R / add_rank_discordance() (still max_tpm-based).

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  group_by(rep) %>%
  mutate(
    rr = dense_rank(mean) / n(),
    re = row_number(hek_tpm) / n(),
    avg_rank = (rr + re) / 2,
    dif_signed = rr - re
  ) %>%
  ungroup()

n_promoters <- data %>% distinct(seq_id) %>% nrow()
message("n = ", n_promoters, " HEK293-active promoters")

binary_df <- remap_hits %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, dif_signed, avg_rank) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, dif_signed, avg_rank), ~ replace_na(.x, FALSE)))

# Same >100-in-both-reps, non-histone filter as remap_tf.R
nTF <- data_TF %>%
  group_by(rep) %>%
  summarise(across(where(is.logical), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(-rep) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  filter(!str_detect(name, "^H\\d")) %>%
  distinct(name)

message(nrow(nTF), " TFs pass the >100-promoters-in-both-reps filter (HEK293-active set)")

data_TF <- data_TF %>%
  mutate(
    rep = as.factor(rep),
    across(all_of(nTF$name), ~ factor(.x, levels = c("TRUE", "FALSE")))
  )

resultados_lmm <- map(nTF$name, function(var) {
  formula <- as.formula(paste("dif_signed ~ avg_rank +", paste0("`", var, "`"), "+ (1 | rep)"))
  fit <- lmerTest::lmer(formula, data = data_TF, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = var)
}) %>%
  list_rbind() %>%
  mutate(
    conf.low_true = -conf.high,
    conf.high_true = -conf.low,
    estimate = -estimate,
    conf.low = conf.low_true,
    conf.high = conf.high_true
  ) %>%
  select(-conf.low_true, -conf.high_true) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(resultados_lmm, file.path(out_dir, "remap_tf_rank_dif_lmm_hek_endo_excl_inactive.tsv"))

n_sub <- sum(resultados_lmm$p_adj < 0.05 & resultados_lmm$estimate < 0)
n_sobre <- sum(resultados_lmm$p_adj < 0.05 & resultados_lmm$estimate > 0)
message(n_sub, " subestima / ", n_sobre, " sobreestima (FDR<0.05) of ", nrow(resultados_lmm), " TFs tested")

message("Saved to ", out_dir)
