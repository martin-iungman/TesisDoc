# discordancia_reportero_endogeno / 2_dif_signed_bland_altman / puffin_motifs
#
# Analogo a prom_df_features.R / remap_tf.R (misma carpeta) pero con los
# scores de efecto por motivo de PUFFIN (Dudnyk et al. 2024) en la
# posicion del TSS de nuestra propia library como predictor continuo en
# vez de una feature booleana: para cada motivo, dif_signed ~ avg_rank +
# motivo + (1|replica). A diferencia de TFs/features (factor TRUE/FALSE,
# con el cuidado de signo documentado en remap_tf.R), estos son
# predictores CONTINUOS - el signo del coeficiente se lee directo, sin
# contraste de referencia: positivo = a mas efecto de ese motivo, mas
# sobreestima el reportero; negativo = mas subestima (dif_signed = rr -
# re, convencion elegida por el autor 2026-07-31).
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv, data/processed/fantom_endo_activity_summary.tsv
# (EXCEPCION, ver R/00_prom_features/analysis_tables_exceptions.R) y
# path_puffin_pred (Analysis/Tables/Dudnyk_puffin_prediction_summ.tsv,
# EXCEPCION - el pipeline que la genera desde cero esta portado pero NO
# ejecutado en R/25_puffin/puffin_prediction_model.R). Run from the
# TesisDoc repo root.

library(tidyverse)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/analysis_tables_exceptions.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance()

# --- Scores de efecto por motivo de PUFFIN (EXCEPCION) ---------------------

df_pred <- read_tsv(path_puffin_pred, show_col_types = FALSE)
motif_cols_raw <- c(
  "CREB+ motif effect", "CREB- motif effect", "ETS+ motif effect", "ETS- motif effect",
  "NFY+ motif effect", "NFY- motif effect", "NRF1+ motif effect", "NRF1- motif effect",
  "SP+ motif effect", "SP- motif effect", "TATA+ motif effect", "TATA- motif effect",
  "U1 snRNP+ motif effect", "U1 snRNP- motif effect", "YY1+ motif effect", "YY1- motif effect",
  "ZNF143+ motif effect", "ZNF143- motif effect",
  "Sum of motif effect", "Sum of initiator effect", "Sum of trinucleotide effect", "Sum of total effect"
)
# Nombres validos para R (sin espacios/+/-), guardando el nombre original
# legible aparte para los graficos.
motif_lookup <- tibble(raw = motif_cols_raw) %>%
  mutate(clean = raw %>% str_replace_all("\\+", "_pos") %>% str_replace_all("-", "_neg") %>% str_replace_all("[ /]", "_"))

df_pred_clean <- df_pred %>%
  select(seq_id, all_of(motif_cols_raw)) %>%
  rename_with(~ motif_lookup$clean[match(.x, motif_lookup$raw)], .cols = all_of(motif_cols_raw))

data <- data %>% left_join(df_pred_clean, by = "seq_id") %>% mutate(rep = as.factor(rep))

# --- LMM: efecto de cada score de motivo PUFFIN sobre la discordancia -----

fit_motif <- function(var) {
  df <- data %>% filter(!is.na(.data[[var]]))
  formula <- as.formula(paste("dif_signed ~ avg_rank +", paste0("`", var, "`"), "+ (1 | rep)"))
  fit <- lmerTest::lmer(formula, data = df, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = var, n = nrow(df))
}

resultados_lmm <- map(motif_lookup$clean, fit_motif) %>%
  list_rbind() %>%
  left_join(motif_lookup, by = c("variable" = "clean")) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

write_tsv(resultados_lmm, file.path(out_dir, "puffin_motifs_rank_dif_lmm.tsv"))

# --- Grafico: estimado (beta) +- IC 95%, coloreado por significancia (BH) -

p <- resultados_lmm %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "FDR < 0.001",
      p_adj < 0.01 ~ "FDR < 0.01",
      p_adj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "ns"
    ),
    sig = factor(sig, levels = c("FDR < 0.001", "FDR < 0.01", "FDR < 0.05", "ns")),
    raw = fct_reorder(raw, estimate)
  ) %>%
  ggplot(aes(x = estimate, y = raw, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, alpha = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70"
  )) +
  labs(
    title = "Efecto de cada score de motivo PUFFIN (en el TSS) sobre la discordancia\ncon signo entre actividad del reportero y endógena",
    subtitle = "Discordancia (reportero - endo) ~ score de motivo (continuo) + nivel de actividad + (1|réplica)",
    x = "Estimado (β) con IC 95% — positivo: a más score, más sobreestima el reportero", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "puffin_motifs_rank_dif_lmm.jpg"), p, width = 9.5, height = 7, units = "in")

message("Figuras guardadas en ", out_dir)
