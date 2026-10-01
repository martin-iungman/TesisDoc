# English version of rank_dif_lmm_hek_endo_excl_inactive (see
# prom_df_features_hek_endo.R, same folder, for the model). HEK293-active
# promoters only (hek_tpm > 0), further excluding:
#   - "Non detected activity (FANTOM5)" (floor/tie artifact, see session
#     2026-09-11 discussion)
#   - "Retrotransposon" (excluded per author request, as in the other
#     English figures)
#   - any feature with fewer than 100 promoters carrying it (TRUE), among
#     the HEK293-active set - too small/unstable a group to report
#
# Does not refit the model - reads rank_dif_lmm_hek_endo_excl_inactive.tsv
# and recomputes group sizes from the same HEK293-active data to apply the
# n>=100 filter.
#
# Output goes to ../transcriptional_library/Plots/Review/.

library(tidyverse)
source("R/functions/plot_helpers.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
plots_review_dir <- "../transcriptional_library/Plots/Review"

resultados_lmm <- read_tsv(file.path(out_dir, "rank_dif_lmm_hek_endo_excl_inactive.tsv"), show_col_types = FALSE)

# --- Recompute per-feature group sizes on the same HEK293-active set ------

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data_excl <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  filter(!is.na(hek_tpm), hek_tpm > 0)

tidy_data <- build_tidy_features(data_excl, keep = "seq_id")
n_true <- tidy_data %>%
  distinct(seq_id, .keep_all = TRUE) %>%
  select(-seq_id, -rep) %>%
  summarise(across(everything(), ~ sum(as.character(.x) == "TRUE"))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_true")

small_features <- n_true %>% filter(n_true < 100) %>% pull(variable)
message("Excluded for n<100: ", paste(small_features, collapse = ", "))

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
  "Enhancers a 50kb" = "Enhancers at 50kb window"
  # "Alta actividad en HEK293" intentionally omitted - not in this model
  # (tautological with the hek_tpm-based outcome, see prom_df_features_hek_endo.R)
)

excluded_features <- c("No detectado (FANTOM5)", "Retrotransposón", small_features)

resultados_lmm <- resultados_lmm %>%
  filter(!variable %in% excluded_features) %>%
  mutate(variable = recode(variable, !!!label_en))

n_promoters <- data_excl %>% distinct(seq_id) %>% nrow()

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
ggsave(file.path(plots_review_dir, "rank_dif_lmm_hek_endo_excl_inactive_english.jpg"), p, width = 9, height = 6.75, units = "in")
ggsave(file.path(plots_review_dir, "rank_dif_lmm_hek_endo_excl_inactive_english.pdf"), p, width = 9, height = 6.75, units = "in", device = cairo_pdf)

message("Figures saved to ", plots_review_dir)
