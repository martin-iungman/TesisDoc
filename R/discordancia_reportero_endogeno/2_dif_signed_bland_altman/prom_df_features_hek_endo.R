# Variant of prom_df_features.R (same folder) using HEK293 CAGE activity
# (hek_tpm) as the endogenous side instead of max TPM across all FANTOM5
# samples (max_tpm). Motivation (session 2026-09-11): re = rank(max_tpm)
# is mechanically inflated for promoters detected across many samples
# (breadth of expression), which is plausibly why CpG islands/high
# conservation/broad promoters/high G+C showed up as "subestima" in the
# max_tpm design - not because the reporter is missing real context, but
# because the endogenous "ceiling" itself is a max over ~1800 draws.
# hek_tpm is a single, matched-context measurement (same cell line as the
# reporter assay), so it doesn't have that inflation - but ~46% of
# promoters show hek_tpm==0 (not detected in HEK293 CAGE), so we run BOTH:
# (1) including them (re_hek has a large tied block at 0 - same kind of
# floor risk as "No detectado (FANTOM5)" in the max_tpm design) and
# (2) restricted to hek_tpm>0 (smaller n, no floor block, but drops
# promoters not naturally active in HEK293).
#
# `Alta actividad en HEK293` is EXCLUDED from the feature list here: it's
# a boolean derived directly from hek_tpm, which is now half of avg_rank
# and the other term of dif_signed - testing it would be tautological
# (same reason "alta actividad del reportero" was dropped from the
# original v1 sensitivity control, see prom_df_features.R header).
#
# Does not touch add_rank_discordance() (still max_tpm-based, used
# elsewhere) - defines rr/re/avg_rank/dif_signed locally instead.

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

data_all <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

add_rank_discordance_hek <- function(df) {
  df %>%
    group_by(rep) %>%
    mutate(
      rr = dense_rank(mean) / n(),
      re = row_number(hek_tpm) / n(),
      avg_rank = (rr + re) / 2,
      dif_signed = rr - re
    ) %>%
    ungroup()
}

run_lmm <- function(data, label) {
  n_prom <- data %>% distinct(seq_id) %>% nrow()
  message(label, ": n = ", n_prom, " promotores")

  tidy_data <- build_tidy_features(data, keep = c("seq_id", "dif_signed", "avg_rank"))
  vars_dicotomicas <- tidy_data %>%
    select(-seq_id, -rep, -dif_signed, -avg_rank) %>%
    select(-`Alta actividad en HEK293`) %>% # tautological here, see header
    names()

  resultados_lmm <- map(vars_dicotomicas, function(var) {
    formula <- as.formula(paste("dif_signed ~ avg_rank +", paste0("`", var, "`"), "+ (1 | rep)"))
    fit <- lmerTest::lmer(formula, data = tidy_data, REML = FALSE)
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

  write_tsv(resultados_lmm, file.path(out_dir, paste0("rank_dif_lmm_hek_endo_", label, ".tsv")))

  p <- resultados_lmm %>%
    mutate(
      sig = case_when(
        p_adj < 0.001 ~ "FDR < 0.001",
        p_adj < 0.01 ~ "FDR < 0.01",
        p_adj < 0.05 ~ "FDR < 0.05",
        TRUE ~ "ns"
      ),
      sig = factor(sig, levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns")),
      variable = fct_reorder(variable, estimate)
    ) %>%
    ggplot(aes(x = estimate, y = variable, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, alpha = 0.7) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70"
    )) +
    labs(
      title = paste0("Efecto de features sobre discordancia (reportero - endógena en HEK293) — ", label),
      subtitle = paste0("re = rank(hek_tpm) en vez de max_tpm. n = ", n_prom, " promotores"),
      x = "Estimado (β) con IC 95% — positivo: el reportero sobreestima", y = NULL, color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, paste0("rank_dif_lmm_hek_endo_", label, ".jpg")), p, width = 9, height = 6.75, units = "in")

  resultados_lmm
}

# --- (1) including HEK293-inactive promoters (hek_tpm == 0) ---------------
data_incl <- data_all %>% mutate(hek_tpm = replace_na(hek_tpm, 0)) %>% add_rank_discordance_hek()
res_incl <- run_lmm(data_incl, "incl_inactive")

# --- (2) excluding HEK293-inactive promoters (hek_tpm > 0 only) -----------
data_excl <- data_all %>% filter(!is.na(hek_tpm), hek_tpm > 0) %>% add_rank_discordance_hek()
res_excl <- run_lmm(data_excl, "excl_inactive")

message("Figuras guardadas en ", out_dir)
