# English main dot-plot for remap_tf_maxtpm_hek_active.R (same folder):
# max_tpm endogenous, ChIP from any cell line, HEK293-active promoters
# only. Same style as remap_tf_english.R / remap_tf_hek_endo_plot.R.
# Does not refit the model - reads the already-computed TSV.
#
# Output goes to ../transcriptional_library/Plots/Review/.

library(tidyverse)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
plots_review_dir <- "../transcriptional_library/Plots/Review"

resultados_lmm <- read_tsv(file.path(out_dir, "remap_tf_rank_dif_lmm_maxtpm_hek_active.tsv"), show_col_types = FALSE)
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
    x = "Estimate (β) — positive: reporter overestimates",
    y = paste0("Transcription factors (ReMap, any cell line, n=", n_tf, ")"), color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

dir.create(plots_review_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plots_review_dir, "remap_rank_dif_lmm_maxtpm_hek_active_english.jpg"), p, width = 9, height = 6.75, units = "in")
ggsave(file.path(plots_review_dir, "remap_rank_dif_lmm_maxtpm_hek_active_english.pdf"), p, width = 9, height = 6.75, units = "in", device = cairo_pdf)

message("Figures saved to ", plots_review_dir)
