# English version of remap_rank_dif_lmm (see remap_tf.R, same folder, for
# the Spanish original and the full model/GSEA/ORA) - all 986 TFs, no
# HEK293-activity filter. Relabels the already-computed results (TF gene
# symbols need no translation) with English axis/legend text, matching
# prom_df_features_english.R / prom_df_features_hek_active_english.R.
# Does not refit the model - reads remap_tf_rank_dif_lmm.tsv.
#
# Output goes to ../transcriptional_library/Plots/Review/ (English-labelled
# figures live there, not in this repo's figures/ tree).

library(tidyverse)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
plots_review_dir <- "../transcriptional_library/Plots/Review"

resultados_lmm <- read_tsv(file.path(out_dir, "remap_tf_rank_dif_lmm.tsv"), show_col_types = FALSE)

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
    y = "Transcription factors (ReMap, n=986)", color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

dir.create(plots_review_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plots_review_dir, "remap_rank_dif_lmm_english.jpg"), p, width = 9, height = 6.75, units = "in")
ggsave(file.path(plots_review_dir, "remap_rank_dif_lmm_english.pdf"), p, width = 9, height = 6.75, units = "in", device = cairo_pdf)

message("Figures saved to ", plots_review_dir)
