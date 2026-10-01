# discordancia_reportero_endogeno / 2_dif_signed_bland_altman / element_core_promoter
#
# EXPLORATORIO (2026-07-31, a pedido del autor) - NO forma parte de
# prom_df/build_prom_features.R todavia. Analogo a prom_df_features.R
# (misma carpeta) pero con los elementos core de promotor llamados por
# ElemeNT sobre EPDnew (Sloutskin et al., CORE_EPDnew_Aug2023.xlsx,
# hoja "human") en vez de las ~32 features curadas: para cada elemento
# (TATA, Inr, DPE, MTE, BRE, etc.), dif_signed ~ avg_rank + elemento +
# (1|replica), correccion BH sobre el efecto fijo.
#
# Portado de la seccion "### ElemeNT" de old/seq_features.qmd (nunca
# incorporada al pipeline numerado): alli cada celda es "no" o una lista
# "posicion,score;..." de sitios encontrados: se colapsa a presencia/
# ausencia (T si hay al menos un sitio, F si "no"), igual que el script
# viejo. No se usa la posicion/score de cada sitio - solo presencia.
#
# `estimate` YA esta en la convencion directa TRUE-menos-FALSE (ver nota
# de signo en prom_df_features.R/remap_tf.R, misma carpeta, sesion
# 2026-07-31): dif_signed = rr - re, asi que positivo = promotores CON
# ese elemento tienen dif_signed mas alto = el reportero SOBREESTIMA;
# negativo = SUBESTIMA.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/
# prom_df.tsv, data/processed/fantom_endo_activity_summary.tsv
# (EXCEPCION, ver R/00_prom_features/analysis_tables_exceptions.R) y
# path_element_epd (heavy_data_paths.R - EPDnew, no copiado a TesisDoc).
# Run from the TesisDoc repo root.

library(tidyverse)
library(readxl)
library(lmerTest)
library(broom.mixed)
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  add_rank_discordance()

# --- Elementos core de promotor (ElemeNT/EPDnew, EXCEPCION) ---------------

element <- read_excel(path_element_epd, sheet = "human")
element_cols <- setdiff(names(element), "name")

element_bool <- element %>%
  mutate(across(all_of(element_cols), ~ .x != "no")) %>%
  rename_with(~ paste0(.x, "_element"), all_of(element_cols))

data <- data %>%
  left_join(element_bool, by = "name") %>%
  mutate(
    rep = as.factor(rep),
    across(ends_with("_element"), ~ factor(.x, levels = c("TRUE", "FALSE")))
  )

vars_element <- paste0(element_cols, "_element")

# --- LMM: efecto de cada elemento sobre la discordancia con signo --------

resultados_lmm <- map(vars_element, function(var) {
  df <- data %>% filter(!is.na(.data[[var]]))
  n_rows <- nrow(df)
  n_true_rows <- sum(as.character(df[[var]]) == "TRUE")
  formula <- as.formula(paste(
    "dif_signed ~ avg_rank +",
    paste0("`", var, "`"), "+ (1 | rep)"
  ))
  fit <- lmerTest::lmer(formula, data = df, REML = FALSE)
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    rename_with(~"p_value", matches("p[._]value")) %>%
    filter(grepl(var, term, fixed = TRUE)) %>%
    mutate(variable = str_remove(var, "_element$"), n = n_rows, n_true = n_true_rows)
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

write_tsv(resultados_lmm, file.path(out_dir, "element_core_promoter_rank_dif_lmm.tsv"))
message("n de promotores con cada elemento (TRUE):")
print(resultados_lmm %>% select(variable, n_true, n))

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
    title = "Efecto de cada elemento core de promotor (ElemeNT/EPDnew) sobre la\ndiscordancia (con signo) entre actividad del reportero y endógena",
    subtitle = "Discordancia (reportero - endo) ~ elemento (presencia/ausencia) + nivel de actividad + (1|réplica) — EXPLORATORIO",
    x = "Estimado (β) con IC 95% — positivo: el reportero sobreestima", y = NULL, color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "element_core_promoter_rank_dif_lmm.jpg"), p, width = 9, height = 6.5, units = "in")

message("Figuras guardadas en ", out_dir)
