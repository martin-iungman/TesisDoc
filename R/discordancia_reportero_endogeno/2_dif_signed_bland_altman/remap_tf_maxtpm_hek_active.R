# max_tpm as the endogenous definition (re = rank(max_tpm), same as the
# canonical remap_tf.R), but restricted to HEK293-active promoters
# (hek_tpm > 0) - same population filter as rank_dif_lmm_hek_active.R /
# prom_df_features_hek_active_english.R, applied here to the full TF
# screen. ChIP-seq binding from ANY cell line (remap_tf_hits.tsv, same
# as the canonical remap_tf.R) - this is the "max_tpm + all ReMap ChIP +
# HEK293-active only" combination.
#
# For the HEK293/HEK293T-ChIP-restricted counterpart, see
# remap_tf_maxtpm_hek_chip_hek_active.R (same folder).

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  add_rank_discordance()

n_promoters <- data %>% distinct(seq_id) %>% nrow()
message("n = ", n_promoters, " HEK293-active promoters (endo = max_tpm)")

binary_df <- remap_hits %>%
  distinct(seq_id, TF) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = FALSE)

data_TF <- data %>%
  select(seq_id, rep, dif_signed, avg_rank) %>%
  left_join(binary_df, by = "seq_id") %>%
  mutate(across(-c(seq_id, rep, dif_signed, avg_rank), ~ replace_na(.x, FALSE)))

nTF <- data_TF %>%
  group_by(rep) %>%
  summarise(across(where(is.logical), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(-rep) %>%
  filter(value > 100) %>%
  add_count(name) %>%
  filter(n == 2) %>%
  filter(!str_detect(name, "^H\\d")) %>%
  distinct(name)

message(nrow(nTF), " TFs pass the >100-promoters-in-both-reps filter")

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

write_tsv(resultados_lmm, file.path(out_dir, "remap_tf_rank_dif_lmm_maxtpm_hek_active.tsv"))

n_sub <- sum(resultados_lmm$p_adj < 0.05 & resultados_lmm$estimate < 0)
n_sobre <- sum(resultados_lmm$p_adj < 0.05 & resultados_lmm$estimate > 0)
message(n_sub, " subestima / ", n_sobre, " sobreestima (FDR<0.05) of ", nrow(resultados_lmm), " TFs tested")

message("\nTFs sobreestima (FDR<0.05):")
resultados_lmm %>% filter(p_adj < 0.05, estimate > 0) %>% arrange(desc(estimate)) %>%
  select(variable, estimate, p_adj) %>% print(n = 100)

message("Saved to ", out_dir)
