# Fully matched-context variant of remap_tf_hek_endo.R (same folder):
# ChIP-seq binding restricted to peaks with HEK293/HEK293T among the
# contributing cell lines (data/processed/remap_tf_hits_hek293.tsv,
# session 2026-09-11), combined with dif_signed defined against HEK293
# CAGE activity (hek_tpm, HEK293-active promoters only) - both the
# predictor (TF binding) and the endogenous side of the outcome now come
# from the same cell line, instead of pooling ChIP-seq across every cell
# line ReMap has data for.
#
# Same >100-promoters-in-both-reps filter as remap_tf.R/remap_tf_hek_endo.R,
# applied to the HEK293/HEK293T-restricted hits.

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

remap_hits_hek <- read_tsv("data/processed/remap_tf_hits_hek293.tsv", show_col_types = FALSE)
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

binary_df <- remap_hits_hek %>%
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

message(nrow(nTF), " TFs pass the >100-promoters-in-both-reps filter (HEK293/HEK293T ChIP, HEK293-active promoters)")

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

write_tsv(resultados_lmm, file.path(out_dir, "remap_tf_rank_dif_lmm_hek_chip_hek_endo.tsv"))

n_sub <- sum(resultados_lmm$p_adj < 0.05 & resultados_lmm$estimate < 0)
n_sobre <- sum(resultados_lmm$p_adj < 0.05 & resultados_lmm$estimate > 0)
message(n_sub, " subestima / ", n_sobre, " sobreestima (FDR<0.05) of ", nrow(resultados_lmm), " TFs tested")

# --- Repressors of interest -------------------------------------------

repressors <- c("RNF2", "EZH2", "SUZ12", "CBX2", "CBX4", "CBX7", "JARID2", "TRIM28", "SETDB1")
message("\nRepressor coefficients (matched HEK293/HEK293T ChIP + HEK293 endogenous):")
resultados_lmm %>% filter(variable %in% repressors) %>%
  select(variable, estimate, p_adj) %>%
  print(n = 20)

# --- Main dot-plot, same style as remap_tf_hek_endo_plot.R -----------------

sig_colors <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")
n_tf <- nrow(resultados_lmm)

p <- resultados_lmm %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "FDR < 0.001",
      p_adj < 0.01 ~ "FDR < 0.01",
      p_adj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "ns"
    ) %>% factor(levels = names(sig_colors)),
    variable = fct_reorder(variable, estimate)
  ) %>%
  ggplot(aes(x = estimate, y = variable, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "Efecto de la unión de cada TF (ReMap, ChIP-seq en HEK293/HEK293T)\nsobre la discordancia (con signo) vs. actividad endógena en HEK293",
    subtitle = "Discordancia (reportero - endo en HEK293) ~ TF + nivel de actividad + (1|réplica) — ChIP y endógeno, misma línea celular",
    x = "Estimado (β) — positivo: el reportero sobreestima",
    y = paste0("Factores de transcripción (ChIP-seq en HEK293/HEK293T, n=", n_tf, ")"), color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(file.path(out_dir, "remap_rank_dif_lmm_hek_chip_hek_endo.jpg"), p, width = 9, height = 6.75, units = "in")

message("Saved to ", out_dir)
