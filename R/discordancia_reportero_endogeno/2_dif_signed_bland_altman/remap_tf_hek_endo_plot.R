# Main dot-plot for the hek_tpm-based TF screen (remap_tf_hek_endo.R,
# same folder) - same style as remap_tf.R's panel_all, but using
# dif_signed defined against HEK293 CAGE activity (HEK293-active
# promoters only) instead of max_tpm across all FANTOM5 samples.
# The >100-promoters-in-both-reps filter is already applied upstream in
# remap_tf_hek_endo.R (954 of 986 TFs pass it on this smaller, HEK293-
# active subset) - this script only plots the already-filtered results.

library(tidyverse)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
resultados_lmm <- read_tsv(file.path(out_dir, "remap_tf_rank_dif_lmm_hek_endo_excl_inactive.tsv"), show_col_types = FALSE)

n_tf <- nrow(resultados_lmm)
sig_colors <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")

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
    title = "Efecto de la unión de cada TF (ReMap) sobre la discordancia (con signo)\nentre actividad del reportero y endógena en HEK293",
    subtitle = "Discordancia (reportero - endo en HEK293) ~ TF + nivel de actividad + (1|réplica) — solo promotores activos en HEK293",
    x = "Estimado (β) — positivo: el reportero sobreestima",
    y = paste0("Factores de transcripción (ReMap, n=", n_tf, ")"), color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

ggsave(file.path(out_dir, "remap_rank_dif_lmm_hek_endo_excl_inactive.jpg"), p, width = 9, height = 6.75, units = "in")

message("Saved to ", out_dir)
