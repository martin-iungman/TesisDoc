# promalt_coocurrencia (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Co-ocurrencia entre las tres categorias principales de promotor
# alternativo (Principal/Secundario/Sin promotores alternativos) y las
# features booleanas del promotor (secuencia + contexto endogeno,
# mismas que M14/coocurrencia_motivos) - heatmaps de coeficiente phi y
# similitud de Jaccard.
#
# No es una figura del paper - agregada a pedido del autor como
# extension de R5.1/R5.2 (categorias) y M14 (features), reusando
# build_cooccurrence_features()/plot_cooccurrence_phi()/
# plot_cooccurrence_jaccard() de R/functions/plot_helpers.R.
#
# Requiere: data/processed/prom_df.tsv y data/processed/
# prom_alt_classification.tsv (ver R/00_prom_features/
# build_prom_alt_classification.R). Run from the TesisDoc repo root.

library(tidyverse)
library(fastDummies)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")

slug <- "promalt_coocurrencia"
out_dir <- fig_dir(slug)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
promalt <- read_tsv("data/processed/prom_alt_classification.tsv", show_col_types = FALSE)

# Mismas 3 categorias que R5.1/R5.2 (Fig 4A/4B) - se excluyen non_detected
# y unclassified, que no pertenecen a ninguna de las tres.
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
  distinct(seq_id, prom_alt_cat) %>%
  dummy_cols("prom_alt_cat", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = TRUE) %>%
  mutate(across(-seq_id, ~ replace_na(.x, FALSE)))

tidy_data <- build_cooccurrence_features(prom_df) %>%
  inner_join(promalt_cat, by = "seq_id")

feat_mat <- tidy_data %>%
  select(-seq_id) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# "No detectado (FANTOM5)" queda constante (siempre FALSE) en este
# subconjunto - promalt_cat ya excluyo los promotores "non_detected" al
# construir las 3 categorias, asi que ninguno de los que quedan puede
# ser "no detectado". Una columna sin varianza da NA en toda su fila/
# columna de la matriz de correlacion, lo que arrastra a cooccurrence_
# cluster_levels() a marcar TODAS las demas columnas como "con NA"
# (contaminadas por esa unica columna) y falla el clustering. Se
# descarta cualquier columna sin varianza antes de armar la matriz.
zero_var <- apply(feat_mat, 2, function(x) length(unique(x)) <= 1)
if (any(zero_var)) {
  message("Excluyendo columnas sin varianza en este subconjunto: ", paste(colnames(feat_mat)[zero_var], collapse = ", "))
  feat_mat <- feat_mat[, !zero_var]
}

titulo_phi <- "Co-ocurrencia: categoría de promotor alternativo vs. features (coeficiente phi)"
titulo_jaccard <- "Co-ocurrencia: categoría de promotor alternativo vs. features (similitud Jaccard)"

ggsave(file.path(out_dir, "promalt_cooccurrence_phi.jpg"), plot_cooccurrence_phi(feat_mat, titulo_phi), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "promalt_cooccurrence_jaccard.jpg"), plot_cooccurrence_jaccard(feat_mat, titulo_jaccard), width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
