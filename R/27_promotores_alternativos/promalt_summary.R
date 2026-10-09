# promalt_summary (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Fig. 4C del paper (centro + derecha - falta el esquema a mano de la
# izquierda, correlated/independent/switch). Efecto de ser Principal o
# Secundario, dentro de cada clase de par (switch/correlated/
# independent), sobre la actividad media (barplot Wilcoxon, similar a
# summary_features_secuencia) y sobre el ruido (barplot AUC-ROC, similar a R3.2/R3.3/
# R3.4). Las diferencias son maximas para pares independientes y casi
# nulas para pares switch, tanto en actividad como en ruido.
#
# CORRECCION: el original tiene el mismo typo/bug de direccion en pROC
# que R3.3/R3.4 - factor(level=c("TRUE","FALSE")) (typo "level" en vez
# de "levels", que no hace nada y deja el orden alfabetico por defecto
# de un logical, invirtiendo direction=">"). Corregido con
# factor(.x, levels=c("TRUE","FALSE")), mismo patron que R3.3/R3.4.
#
# SIMPLIFICACION: el original arma df_split_main uniendo wide_df con
# `data %>% left_join(promalt_df) %>% select(-prom_alt)` - ese left_join
# no tiene ningun efecto util (la columna se descarta enseguida) salvo
# duplicar filas para los ~62 promotores que son a la vez "unclassified"
# y "independent/correlated/switch" (ver build_prom_alt_classification.R).
# Se omite ese left_join sin efecto y se une directo contra
# activity_stats_highconf.tsv.
#
# Requiere: data/processed/activity_stats_highconf.tsv y data/processed/
# prom_alt_pairs.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (seccion "Alternative promoters", lineas 1440-1558).

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(coin)
library(pROC)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "promalt_summary"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
pairs_df <- read_tsv("data/processed/prom_alt_pairs.tsv", show_col_types = FALSE) %>%
  mutate(prom_alt = recode(prom_alt, switch = "Alternancia", correlated = "Correlacionado", independent = "Independiente"))

# Un promotor Main puede tener varios secundarios de distinta clase (p.ej.
# "Alternancia principal" para uno y "Correlacionado principal" para otro)
# - cada combinacion se cuenta por separado, como en el original.
promalt_pairs <- pairs_df %>%
  distinct(main_name, prom_alt) %>%
  rename(name = main_name) %>%
  mutate(prom_alt = paste(prom_alt, "principal")) %>%
  bind_rows(pairs_df %>% distinct(name, prom_alt))

# --- Panel centro: efecto sobre la actividad (Wilcoxon) --------------------

tidy_palt_act <- promalt_pairs %>%
  inner_join(stats_highconf %>% select(name, mean, rep), by = "name", relationship = "many-to-many") %>%
  select(prom_alt, rep, mean) %>%
  dummy_cols("prom_alt", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = TRUE) %>%
  mutate(
    across(c(where(is.numeric), -contains("mean")), as.logical),
    across(where(is.logical), ~ .x %>% factor(levels = c("TRUE", "FALSE"))),
    rep = as.factor(rep)
  )

wilcox_df <- wilcox_effect_summary(tidy_palt_act)

panel_act <- wilcox_df %>%
  pivot_wider(names_from = val, values_from = estimate) %>%
  mutate(
    act = ifelse(sign(estimate) == 1, "Actividad alta", "Actividad baja"),
    cat = str_remove(feature, " principal"), cat2 = ifelse(str_detect(feature, "principal"), "main", "sec")
  ) %>%
  arrange(desc(cat), desc(rep), desc(cat2)) %>%
  ggplot(aes(x = estimate, y = fct_inorder(feature), group = fct_inorder(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = act, alpha = rep)) +
  scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7)) +
  theme(text = element_text(size = 20), legend.position = "top") +
  labs(fill = "Efecto", x = "Efecto sobre la actividad", y = "Categoría", alpha = "") +
  geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5) +
  xlim(c(-0.3, 0.3)) +
  scale_fill_manual(values = c("Actividad alta" = "#D6741F", "Actividad baja" = "#7FB800")) +
  guides(alpha = "none") +
  theme_pubr()
ggsave(file.path(out_dir, "summary_actividad_promalt.jpg"), panel_act, width = 13.5, height = 6.75, units = "in")

# --- Panel derecha: efecto sobre el ruido (AUC-ROC) -------------------------
# var_rank_sw sobre el set completo de alta confianza (add_noise_rank()),
# mismo criterio que R3.x - no el subconjunto filtrado que usa R5.2/Fig4B
# para el scatter.

data_noise <- add_noise_rank(stats_highconf)

tidy_palt_noise <- promalt_pairs %>%
  inner_join(data_noise %>% select(name, var_rank_sw, rep), by = "name", relationship = "many-to-many") %>%
  select(prom_alt, rep, var_rank_sw) %>%
  dummy_cols("prom_alt", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = TRUE) %>%
  mutate(across(c(where(is.numeric), -contains("var_rank_sw")), as.logical))

vbles <- setdiff(names(tidy_palt_noise), c("rep", "var_rank_sw"))
repname <- unique(tidy_palt_noise$rep)

auc_df <- map(repname, function(r) {
  df <- tidy_palt_noise %>%
    filter(rep == r) %>%
    mutate(across(all_of(vbles), ~ factor(.x, levels = c("TRUE", "FALSE"))))
  map(vbles, function(f) {
    roc_obj <- pROC::roc(df[[f]], df[["var_rank_sw"]], ci = TRUE, direction = ">", quiet = TRUE)
    tibble(feature = f, rep = r, AUC = as.numeric(roc_obj$ci[2]), ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]))
  }) %>% list_rbind()
}) %>% list_rbind()

panel_noise <- auc_df %>%
  mutate(
    noise = ifelse(AUC > 0.5, "Ruido alto", "Ruido bajo"), cat = str_remove(feature, " principal"),
    a = ifelse(str_detect(feature, "Independiente"), 6, ifelse(str_detect(feature, "Correlacionado"), 4, 2)),
    b = ifelse(str_detect(feature, "principal"), 1, 0), c = a + b
  ) %>%
  arrange(c, desc(rep)) %>%
  ggplot(aes(x = AUC - 0.5, y = fct_inorder(feature), group = fct_inorder(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = noise, alpha = rep)) +
  scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7)) +
  theme(text = element_text(size = 20), legend.position = "top") +
  labs(fill = "Efecto", x = "AUC (efecto sobre el ruido)", y = "Categoría", alpha = "") +
  scale_x_continuous(labels = function(x) x + 0.5) +
  geom_errorbarh(aes(xmax = abs(ci2.5) - 0.5, xmin = abs(ci97.5) - 0.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5) +
  scale_fill_manual(values = c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800")) +
  theme_pubr() +
  guides(alpha = "none")
ggsave(file.path(out_dir, "summary_ruido_promalt.jpg"), panel_noise, width = 13.5, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
