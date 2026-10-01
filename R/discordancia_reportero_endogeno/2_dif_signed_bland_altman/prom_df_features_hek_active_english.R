# English version of rank_dif_lmm_hek_active (see prom_df_features.R,
# same folder, for the Spanish original and the full model). Relabels the
# already-computed results with the English feature names used throughout
# the rest of the paper (transcriptional_library/Analysis/Tables/
# tidy_names.tsv), for a figure meant to go alongside the paper panels.
# Does not refit the model - reads rank_dif_lmm_hek_active.tsv.
#
# Output goes to ../transcriptional_library/Plots/Review/ (English-labelled
# figures live there, not in this repo's figures/ tree).

library(tidyverse)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
plots_review_dir <- "../transcriptional_library/Plots/Review"
resultados_lmm <- read_tsv(file.path(out_dir, "rank_dif_lmm_hek_active.tsv"), show_col_types = FALSE)

# Spanish (this repo) -> English (transcriptional_library/Analysis/Tables/tidy_names.tsv)
label_en <- c(
  "No detectado (FANTOM5)" = "Non detected activity (FANTOM5)",
  "LTR" = "LTR",
  "Sin actividad en ratón" = "No activity in mouse",
  "Insertado en humanos" = "Human inserted",
  "Retrotransposón" = "Retrotransposon",
  "Alta accesibilidad de cromatina (DNase-seq)" = "High chromatin accesibility (DNase-seq)",
  "CCAAT" = "CCAAT",
  "LINE" = "LINE",
  "Elementos transponibles" = "Transposable elements",
  "Sin módulo cis-regulatorio" = "None Cis Regulatory Module",
  "Baja especificidad tisular" = "Low tissue specificity",
  "SINE" = "SINE",
  "GC-box" = "GC-box",
  "TA en TSS" = "TA at TSS",
  "Islas CpG" = "CpG islands",
  "TG en TSS" = "TG at TSS",
  "Promotores angostos" = "Narrow promoters",
  "CG en TSS" = "CG at TSS",
  "CA en TSS" = "CA at TSS",
  "TSS no canónico" = "Non canonical TSS",
  "TSS fuerte" = "Strong TSS",
  "Promotores anchos" = "Broad promoters",
  "TCT" = "TCT",
  "Alta conservación (-150 a -235)" = "High conservation (-150 to -235)",
  "Alta especificidad tisular" = "High tissue specificity",
  "Alto contenido G+C" = "High G+C content",
  "Repeticiones de baja complejidad" = "Low Complexity Repeats",
  "Alta conservación (-50 a -150)" = "High conservation (-50 to -150)",
  "Alta conservación (16 a -50pb)" = "High conservation (16 to -50pb)",
  "TATA-box" = "TATA-box",
  "Enhancers a 50kb" = "Enhancers at 50kb window",
  "Alta actividad en HEK293" = "High activity in HEK293",
  "Promotor unidireccional" = "Unidirectional promoter"
)

stopifnot(all(resultados_lmm$variable %in% names(label_en)))

# Excluded per author request: with max_tpm as the endogenous side (this
# model), "Alta actividad en HEK293" is a real, non-tautological feature -
# but it's excluded from this figure for consistency with the hek_tpm-
# endogenous English figures, where it IS tautological and always dropped.
resultados_lmm <- resultados_lmm %>%
  filter(variable != "Alta actividad en HEK293") %>%
  mutate(variable = recode(variable, !!!label_en))

n_promoters <- 6692 # see prom_df_features.R run log (6692 of 12284 active in HEK293)

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
    x = "Effect (β), 95% CI — positive: reporter overestimates", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

dir.create(plots_review_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plots_review_dir, "rank_dif_lmm_hek_active_english.jpg"), p, width = 9, height = 6.75, units = "in")
ggsave(file.path(plots_review_dir, "rank_dif_lmm_hek_active_english.pdf"), p, width = 9, height = 6.75, units = "in", device = cairo_pdf)

message("Figures saved to ", plots_review_dir)
