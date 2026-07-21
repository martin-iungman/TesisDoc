# coocurrencia_motivos (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Co-ocurrencia entre features de promotores: heatmaps de coeficiente phi
# (correlacion de Pearson sobre variables binarias) y similitud de
# Jaccard, cada uno con todas las features y solo con las
# estadisticamente significativas (misma direccion de efecto en ambas
# replicas y ambos limites del IC, ver R7/activity_summary.tsv) - 4
# figuras en total.
#
# Requiere: data/processed/prom_df.tsv (ver R/00_prom_features) y
# data/processed/activity_summary.tsv (ver
# R/15_summary_features/summary_features.R). Run from the TesisDoc repo root.

library(tidyverse)
library(fastDummies)
source("R/functions/fig_paths.R")

slug <- "coocurrencia_motivos"
out_dir <- fig_dir(slug)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
act_summary <- read_tsv("data/processed/activity_summary.tsv", show_col_types = FALSE)
act_signif_features <- act_summary %>%
  mutate(sign = sign(estimate)) %>%
  filter(pval_corr < 0.05) %>%
  count(feature, sign) %>%
  filter(n == 6) %>%
  pull(feature)

# --- Tabla de features booleanas (mismos nombres en espanol que R7) ------

tidy_data <- prom_df %>%
  fastDummies::dummy_cols("TE_superclass", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = FALSE) %>%
  fastDummies::dummy_cols("sample_specificity_class", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = TRUE) %>%
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
    `Alta actividad en HEK293` = hek_tpm > median(prom_df$hek_tpm[prom_df$hek_tpm > 0], na.rm = TRUE),
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
    seq_id,
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
  )

# --- Helpers: matrices de phi / Jaccard + heatmap -------------------------

phi_matrix <- function(feat_mat) cor(feat_mat, use = "pairwise.complete.obs")

jaccard_matrix <- function(feat_mat) {
  m <- matrix(NA, nrow = ncol(feat_mat), ncol = ncol(feat_mat), dimnames = list(colnames(feat_mat), colnames(feat_mat)))
  for (i in seq_len(ncol(feat_mat))) {
    for (j in seq_len(ncol(feat_mat))) {
      a <- feat_mat[, i]
      b <- feat_mat[, j]
      intersection <- sum(a == 1 & b == 1, na.rm = TRUE)
      union <- sum(a == 1 | b == 1, na.rm = TRUE)
      m[i, j] <- if (union == 0) NA else intersection / union
    }
  }
  m
}

cluster_levels <- function(mat) {
  has_na <- colSums(is.na(mat)) > 0
  complete <- mat[!has_na, !has_na]
  clust_order <- hclust(as.dist(1 - complete))$order
  c(colnames(complete)[clust_order], colnames(mat)[has_na])
}

tidy_matrix <- function(mat, value_name) {
  levels <- cluster_levels(mat)
  mat %>%
    as_tibble(rownames = "feature1") %>%
    pivot_longer(-feature1, names_to = "feature2", values_to = value_name) %>%
    mutate(feature1 = factor(feature1, levels = levels), feature2 = factor(feature2, levels = levels))
}

plot_phi <- function(feat_mat) {
  tidy_matrix(phi_matrix(feat_mat), "phi") %>%
    ggplot(aes(x = feature2, y = feature1, fill = phi)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0, limits = c(-1, 1), name = "Coeficiente\nphi") +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
      axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
    ) +
    labs(title = "Co-ocurrencia de features del promotor (coeficiente phi)")
}

plot_jaccard <- function(feat_mat) {
  tidy_matrix(jaccard_matrix(feat_mat), "jaccard") %>%
    ggplot(aes(x = feature2, y = feature1, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#d6604d", limits = c(0, 1), na.value = "grey80", name = "Similitud\nJaccard") +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
      axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
    ) +
    labs(title = "Co-ocurrencia de features del promotor (similitud Jaccard)")
}

save_pair <- function(data_for_matrix, suffix) {
  feat_mat <- data_for_matrix %>%
    select(-seq_id) %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  ggsave(file.path(out_dir, paste0("cooccurrence_phi", suffix, ".jpg")), plot_phi(feat_mat), width = 9, height = 6.75, units = "in")
  ggsave(file.path(out_dir, paste0("cooccurrence_jaccard", suffix, ".jpg")), plot_jaccard(feat_mat), width = 9, height = 6.75, units = "in")
}

# --- 4 figuras: todas las features / solo las significativas -------------

save_pair(tidy_data, "")
save_pair(tidy_data %>% select(seq_id, any_of(act_signif_features)), "_signif")

message("Figuras guardadas en ", out_dir)
