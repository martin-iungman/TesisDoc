# summary_features_secuencia = Fig. R7 (ver docs/mapping_figuras.csv)
# Resumen del efecto de features de secuencia (seq) y de contexto
# endogeno/genomico (endo) sobre la actividad transcripcional media:
# para cada feature booleana, test de Wilcoxon (mean ~ feature, pareado
# por replica) con correccion BH: las que dan significativas se grafican
# como barras con el tamano de efecto (estimate) e IC 95%.
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y
# R/01_activity_stats - incluye las features de
# R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(coin)
source("R/functions/fig_paths.R")

slug <- "summary_features_secuencia"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

# --- Tabla de features booleanas (todas en espanol) -----------------------

tidy_data <- data %>%
  fastDummies::dummy_cols("TE_superclass", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = FALSE) %>%
  fastDummies::dummy_cols("sample_specificity_class", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = TRUE) %>%
  group_by(rep) %>%
  mutate(
    `Alto contenido G+C` = (cut_number(g_c, n = 3) %>% as.numeric()) == 3,
    across(c(DNA, SINE, LINE, LTR), ~ replace_na(.x, 0)),
    `Sin actividad en ratón` = turnover %in% c("expression-turnover", "mouse-diminished"),
    `Insertado en humanos` = turnover == "human-inserted",
    `Elementos transponibles` = !is.na(TE_superclass),
    `Repeticiones de baja complejidad` = LCR_overlap > 0,
    `Alta conservación (16 a -50pb)` = (cut_number(phylop100_close, n = 3) %>% as.numeric()) == 3,
    `Alta conservación (-50 a -150)` = (cut_number(phylop100_intermediate, n = 3) %>% as.numeric()) == 3,
    `Alta conservación (-150 a -235)` = (cut_number(phylop100_far, n = 3) %>% as.numeric()) == 3,
    `Alta especificidad tisular` = (cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) == 3,
    `Baja especificidad tisular` = (cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) == 1,
    shape_class_n = cut_number(interquantile_width, 3) %>% as_factor() %>% as.numeric(),
    `Promotores angostos` = shape_class_n == 1,
    `Promotores anchos` = shape_class_n == 3,
    `Alta actividad en HEK293` = hek_tpm > median(data$hek_tpm[data$hek_tpm > 0], na.rm = TRUE),
    `Sin módulo cis-regulatorio` = N_TF_CRM == 0,
    `Alta accesibilidad de cromatina (DNase-seq)` = (cut_number(mean_dnase, n = 3) %>% as.numeric()) == 3,
    across(starts_with("enh"), ~ .x > 0)
  ) %>%
  rename(
    `TATA-box` = TATA_EPD,
    CCAAT = CCAAT_EPD,
    `GC-box` = GCbox_EPD,
    `Islas CpG` = CGI,
    Retrotransposón = DNA,
    TCT = TCT_TSS,
    `CG en TSS` = CG_TSS,
    `TA en TSS` = TA_TSS,
    `TG en TSS` = TG_TSS,
    `CA en TSS` = CA_TSS,
    `TSS no canónico` = other_TSS,
    `TSS fuerte` = INR_strong_TSS,
    `No detectado (FANTOM5)` = non_detected,
    `Enhancers a 10kb` = enh10kb,
    `Enhancers a 50kb` = enh50kb,
    `Enhancers a 100kb` = enh100kb
  ) %>%
  select(
    mean, rep,
    LINE, SINE, LTR,
    `Alto contenido G+C`, `Sin actividad en ratón`, `Insertado en humanos`,
    `Elementos transponibles`, `Repeticiones de baja complejidad`,
    `Alta conservación (16 a -50pb)`, `Alta conservación (-50 a -150)`, `Alta conservación (-150 a -235)`,
    `TATA-box`, CCAAT, `GC-box`, `Islas CpG`, Retrotransposón, TCT,
    `CG en TSS`, `TA en TSS`, `TG en TSS`, `CA en TSS`,
    `TSS no canónico`, `TSS fuerte`,
    `No detectado (FANTOM5)`,
    `Alta especificidad tisular`, `Baja especificidad tisular`,
    `Promotores angostos`, `Promotores anchos`,
    `Alta actividad en HEK293`, `Alta accesibilidad de cromatina (DNase-seq)`,
    `Enhancers a 50kb`, `Sin módulo cis-regulatorio`
  ) %>%
  mutate(
    across(c(where(is.numeric), -contains("mean")), as.logical),
    rep = as.factor(rep),
    across(where(is.logical), ~ .x %>% factor(levels = c("TRUE", "FALSE")))
  ) %>%
  ungroup()

# seq = de la secuencia del promotor; endo = de contexto endogeno/genomico
# (clasificacion fija, migrada de transcriptional_library/Analysis/Tables/tidy_names.tsv)
feature_groups <- tribble(
  ~feature, ~group,
  "LINE", "seq", "SINE", "seq", "LTR", "seq",
  "Alto contenido G+C", "seq",
  "Elementos transponibles", "seq",
  "Repeticiones de baja complejidad", "seq",
  "TATA-box", "seq", "CCAAT", "seq", "GC-box", "seq", "Islas CpG", "seq",
  "Retrotransposón", "seq", "TCT", "seq",
  "CG en TSS", "seq", "TA en TSS", "seq", "TG en TSS", "seq", "CA en TSS", "seq",
  "TSS no canónico", "seq", "TSS fuerte", "seq",
  "Sin actividad en ratón", "endo", "Insertado en humanos", "endo",
  "Alta conservación (16 a -50pb)", "endo", "Alta conservación (-50 a -150)", "endo", "Alta conservación (-150 a -235)", "endo",
  "No detectado (FANTOM5)", "endo",
  "Alta especificidad tisular", "endo", "Baja especificidad tisular", "endo",
  "Promotores angostos", "endo", "Promotores anchos", "endo",
  "Alta actividad en HEK293", "endo",
  "Alta accesibilidad de cromatina (DNase-seq)", "endo",
  "Enhancers a 50kb", "endo",
  "Sin módulo cis-regulatorio", "endo"
)
vbles_split <- map(c("seq", "endo"), ~ feature_groups$feature[feature_groups$group == .x])

# --- Wilcoxon: efecto de cada feature sobre la actividad media ----------

wilcox <- map(tidy_data %>% select(-rep, -contains("mean")), ~ coin::wilcox_test(formula = mean ~ .x | rep, data = tidy_data) %>% pvalue())
wilcox <- tibble(feature = names(wilcox), pval = list_c(wilcox), pval_corr = p.adjust(pval, "BH", length(wilcox)))

wilcox_by_rep <- function(rep_id) {
  map(tidy_data %>% filter(rep == rep_id) %>% select(-mean, -rep), ~ coin::wilcox_test(formula = mean ~ .x, data = tidy_data %>% filter(rep == rep_id), conf.int = TRUE))
}
wilcox_rep1 <- wilcox_by_rep("Rep 1")
wilcox_rep2 <- wilcox_by_rep("Rep 2")

wilcox_df <- map2(
  list(wilcox_rep1, wilcox_rep2), c("Rep 1", "Rep 2"),
  ~ tibble(
    feature = names(.x),
    estimate = map(.x, ~ confint(.x)$estimate) %>% list_c(),
    P2.5 = map(.x, ~ confint(.x)$conf.int[1]) %>% list_c(),
    P97.5 = map(.x, ~ confint(.x)$conf.int[2]) %>% list_c(),
    rep = .y
  ) %>%
    pivot_longer(c(starts_with("estimate"), starts_with("P2.5"), starts_with("P97.5")), names_to = "val", values_to = "estimate") %>%
    arrange(desc(estimate)) %>%
    mutate(feature = fct_inorder(feature))
) %>%
  list_rbind() %>%
  left_join(wilcox, by = "feature")

# --- Graficos resumen: seq y endo ------------------------------------------

plot_summary <- function(df, features) {
  df %>%
    filter(feature %in% features, pval_corr < 0.05) %>%
    mutate(act = ifelse(sign(estimate) == 1, "Actividad alta", "Actividad baja")) %>%
    group_by(feature) %>%
    mutate(n = length(unique(act))) %>%
    filter(n == 1) %>%
    arrange(estimate, desc(rep)) %>%
    pivot_wider(names_from = val, values_from = estimate) %>%
    filter(pval_corr < 0.05) %>%
    ggplot(aes(x = estimate, y = fct_inorder(feature), group = fct_inorder(rep))) +
    geom_col(orientation = "y", position = "dodge", aes(fill = act, alpha = rep)) +
    scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7)) +
    theme_pubr() +
    theme(text = element_text(size = 20), legend.position = "top") +
    labs(fill = "Efecto", x = "Efecto sobre la actividad", y = "Features", alpha = "") +
    geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5) +
    xlim(c(-0.3, 0.55)) +
    scale_fill_manual(values = c("Actividad alta" = "#D6741F", "Actividad baja" = "#7FB800")) +
    guides(alpha = "none")
}

panel_seq <- plot_summary(wilcox_df, vbles_split[[1]])
ggsave(file.path(out_dir, "summary_seq.jpg"), panel_seq, width = 9, height = 6.75, units = "in")

panel_endo <- plot_summary(wilcox_df, vbles_split[[2]])
ggsave(file.path(out_dir, "summary_endo.jpg"), panel_endo, width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
