# promalt_enrichment_en - English/PDF alternate version of R5.7
# (promalt_enriquecimiento), for external use (e.g. sharing with
# collaborators, paper drafts). Same analysis and data as
# promalt_enriquecimiento.R (see that script for the full
# methodological explanation) - only the labels and output format
# differ.
#
# Feature names match the exact English labels used in
# transcriptional_library/Analysis/scripts/final_github.R (the script
# that generated the actual paper figures), not a fresh translation -
# e.g. "CpG islands", "Non canonical TSS", "Transposable elements",
# "Low Complexity Repeats", "High conservation (16 to -50pb)".
# Category names ("Main promoter"/"Secondary promoter"/"No alternative
# promoters") match final_github.R's own prom_alt2 labels (Fig. 4A/4B).
#
# Not tracked as its own row in docs/mapping_figuras.csv - same figure
# as R5.7 (promalt_enriquecimiento), just relabeled/re-exported.
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

promalt_cat <- promalt %>%
  filter(prom_alt %in% c("unique", "Main promoter", "independent", "correlated", "switch")) %>%
  mutate(prom_alt_cat = case_when(
    prom_alt == "unique" ~ "No alternative promoters",
    prom_alt == "Main promoter" ~ "Main promoter",
    TRUE ~ "Secondary promoter"
  )) %>%
  distinct(name, prom_alt_cat) %>%
  left_join(prom_df %>% distinct(seq_id, name), by = "name", relationship = "many-to-many") %>%
  filter(!is.na(seq_id)) %>%
  distinct(seq_id, prom_alt_cat)

# Mismos nombres en espanol de build_cooccurrence_features(), traducidos
# al ingles EXACTO de final_github.R (no una traduccion libre).
feature_labels_en <- c(
  "Alto contenido G+C" = "High G+C content",
  "Sin actividad en ratón" = "No activity in mouse",
  "Insertado en humanos" = "Human inserted",
  "Elementos transponibles" = "Transposable elements",
  "Repeticiones de baja complejidad" = "Low Complexity Repeats",
  "Alta conservación (16 a -50pb)" = "High conservation (16 to -50pb)",
  "Alta conservación (-50 a -150)" = "High conservation (-50 to -150)",
  "Alta conservación (-150 a -235)" = "High conservation (-150 to -235)",
  "TATA-box" = "TATA-box",
  "CCAAT" = "CCAAT",
  "GC-box" = "GC-box",
  "Islas CpG" = "CpG islands",
  "Retrotransposón" = "Retrotransposon",
  "TCT" = "TCT",
  "CG en TSS" = "CG at TSS",
  "TA en TSS" = "TA at TSS",
  "TG en TSS" = "TG at TSS",
  "CA en TSS" = "CA at TSS",
  "TSS no canónico" = "Non canonical TSS",
  "TSS fuerte" = "Strong TSS",
  "Alta especificidad tisular" = "High tissue specificity",
  "Baja especificidad tisular" = "Low tissue specificity",
  "Promotores angostos" = "Narrow promoters",
  "Promotores anchos" = "Broad promoters",
  "Alta actividad en HEK293" = "High activity in HEK293",
  "Alta accesibilidad de cromatina (DNase-seq)" = "High chromatin accessibility (DNase-seq)",
  "Enhancers a 50kb" = "Enhancers at 50kb window",
  "Sin módulo cis-regulatorio" = "None Cis Regulatory Module",
  "LINE" = "LINE", "SINE" = "SINE", "LTR" = "LTR"
)

# Unidirectional promoter (PRO-seq Orientation Index) - not yet in
# data/processed/prom_df.tsv (TesisDoc), read from transcriptional_library/
# Analysis/Tables/prom_df.tsv instead (see prom_df_features.R, same
# session). not_detected/NA orientation stays NA here (na.rm=TRUE in
# enrichment_cell()/table()'s default NA-drop handle the exclusion).
orientation_df <- read_tsv("../transcriptional_library/Analysis/Tables/prom_df.tsv", show_col_types = FALSE) %>%
  select(seq_id, orientation) %>%
  distinct() %>%
  mutate(`Unidirectional promoter` = case_when(
    orientation == "unidirectional" ~ TRUE,
    orientation == "bidirectional" ~ FALSE,
    TRUE ~ NA
  )) %>%
  select(-orientation)

tidy_data <- build_cooccurrence_features(prom_df) %>%
  select(-`No detectado (FANTOM5)`) %>% # constante en este subconjunto, ver promalt_enriquecimiento.R
  rename(!!!set_names(names(feature_labels_en), feature_labels_en)) %>%
  left_join(orientation_df, by = "seq_id") %>%
  inner_join(promalt_cat, by = "seq_id")

feature_cols <- setdiff(names(tidy_data), c("seq_id", "prom_alt_cat"))

enrichment_cell <- function(df, feature, categoria) {
  x <- as.logical(df[[feature]])
  in_group <- df$prom_alt_cat == categoria
  rest_prop <- mean(x[!in_group], na.rm = TRUE)
  group_prop <- mean(x[in_group], na.rm = TRUE)
  tab <- table(in_group, x)
  pval <- if (all(dim(tab) == c(2, 2))) fisher.test(tab)$p.value else NA_real_
  tibble(feature = feature, categoria = categoria, log2fc = log2(group_prop / rest_prop), pval = pval)
}

enrichment_df <- map(feature_cols, function(f) {
  map(unique(tidy_data$prom_alt_cat), ~ enrichment_cell(tidy_data, f, .x)) %>% list_rbind()
}) %>%
  list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, "BH"))

feature_order <- enrichment_df %>%
  filter(categoria == "Secondary promoter") %>%
  arrange(desc(log2fc)) %>%
  pull(feature)

categoria_order <- c("No alternative promoters", "Main promoter", "Secondary promoter")

p <- enrichment_df %>%
  mutate(
    feature = factor(feature, levels = feature_order),
    categoria = factor(categoria, levels = categoria_order),
    signif = ifelse(pval_corr < 0.05, "*", "")
  ) %>%
  ggplot(aes(x = feature, y = categoria, fill = log2fc)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = signif), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0, name = "log2\nenrichment") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 11),
    axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
  ) +
  labs(title = "Feature enrichment by alternative-promoter category", subtitle = "* = FDR < 0.05 (Fisher's exact test, category vs. rest)")
ggsave(file.path(out_dir, "promalt_enrichment_features_en.pdf"), p, width = 14, height = 4.5, units = "in")

message("PDF guardado en ", out_dir)
