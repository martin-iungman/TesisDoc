# promalt_enriquecimiento (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Que features del promotor estan enriquecidas/depletadas en cada
# categoria de promotor alternativo (Principal/Secundario/Sin
# promotores alternativos, mismo agrupamiento que R5.1/R5.2). Matriz
# rectangular: filas = categoria, columnas = features. Cada celda es el
# log2 fold-enrichment de la prevalencia del feature DENTRO de esa
# categoria respecto a su prevalencia general en el conjunto de
# promotores clasificados (normalizacion por columna - cada feature
# contra su propia tasa basal, ya que features raros como TATA-box y
# comunes como Islas CpG no son comparables en escala absoluta).
# Significancia por test exacto de Fisher (categoria vs. resto, feature
# presente vs. ausente), corregido por Benjamini-Hochberg.
#
# No es una figura del paper - agregada a pedido del autor, numero
# provisorio (mismo criterio que el resto de R3/R4/R5 antes de
# confirmar contra un pptx). Reemplaza a un intento anterior
# (promalt_coocurrencia.R, matriz cuadrada de co-ocurrencia
# feature-vs-feature) que no respondia la pregunta pedida.
#
# Requiere: data/processed/prom_df.tsv y data/processed/
# prom_alt_classification.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R). Run from the TesisDoc repo root.

library(tidyverse)
library(fastDummies)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")

slug <- "promalt_enriquecimiento"
out_dir <- fig_dir(slug)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
promalt <- read_tsv("data/processed/prom_alt_classification.tsv", show_col_types = FALSE)

# Mismas 3 categorias que R5.1/R5.2 (Fig 4A/4B) - se excluyen
# non_detected y unclassified, que no pertenecen a ninguna de las tres.
promalt_cat <- promalt %>%
  filter(prom_alt %in% c("unique", "Main promoter", "independent", "correlated", "switch")) %>%
  mutate(prom_alt_cat = case_when(
    prom_alt == "unique" ~ "Sin promotores alternativos",
    prom_alt == "Main promoter" ~ "Promotor principal",
    TRUE ~ "Promotor secundario"
  )) %>%
  distinct(name, prom_alt_cat) %>%
  left_join(prom_df %>% distinct(seq_id, name), by = "name", relationship = "many-to-many") %>%
  filter(!is.na(seq_id)) %>%
  distinct(seq_id, prom_alt_cat)

tidy_data <- build_cooccurrence_features(prom_df) %>%
  inner_join(promalt_cat, by = "seq_id")

feature_cols <- setdiff(names(tidy_data), c("seq_id", "prom_alt_cat"))
# "No detectado (FANTOM5)" es constante (FALSE) en este subconjunto -
# promalt_cat ya excluyo los promotores "non_detected" - no aporta nada,
# se descarta (ver commit anterior, promalt_coocurrencia.R).
feature_cols <- setdiff(feature_cols, "No detectado (FANTOM5)")

# --- Enriquecimiento (log2 fold-change vs. prevalencia general) + Fisher --

enrichment_cell <- function(df, feature, categoria) {
  x <- as.logical(df[[feature]])
  in_group <- df$prom_alt_cat == categoria
  overall_prop <- mean(x, na.rm = TRUE)
  group_prop <- mean(x[in_group], na.rm = TRUE)
  tab <- table(in_group, x)
  pval <- if (all(dim(tab) == c(2, 2))) fisher.test(tab)$p.value else NA_real_
  tibble(
    feature = feature, categoria = categoria,
    log2fc = log2(group_prop / overall_prop),
    group_prop = group_prop, overall_prop = overall_prop, pval = pval
  )
}

enrichment_df <- map(feature_cols, function(f) {
  map(unique(tidy_data$prom_alt_cat), ~ enrichment_cell(tidy_data, f, .x)) %>% list_rbind()
}) %>%
  list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, "BH"))

write_tsv(enrichment_df, file.path(out_dir, "promalt_enrichment.tsv"))

# --- Heatmap: filas = categoria, columnas = features (por log2fc de
# Promotor secundario, de mas a menos enriquecido) -------------------------

feature_order <- enrichment_df %>%
  filter(categoria == "Promotor secundario") %>%
  arrange(desc(log2fc)) %>%
  pull(feature)

categoria_order <- c("Sin promotores alternativos", "Promotor principal", "Promotor secundario")

p <- enrichment_df %>%
  mutate(
    feature = factor(feature, levels = feature_order),
    categoria = factor(categoria, levels = categoria_order),
    signif = ifelse(pval_corr < 0.05, "*", "")
  ) %>%
  ggplot(aes(x = feature, y = categoria, fill = log2fc)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = signif), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0, name = "log2\nenriquecimiento") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 11),
    axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
  ) +
  labs(title = "Enriquecimiento de features por categoría de promotor alternativo", subtitle = "* = FDR < 0.05 (test exacto de Fisher, categoría vs. resto)")
ggsave(file.path(out_dir, "promalt_enriquecimiento_features.jpg"), p, width = 14, height = 4.5, units = "in")

message("Figura guardada en ", out_dir)
