# Co-occurrence of promoter orientation (Unidirectional, from PRO-seq
# Orientation Index - see ../transcriptional_library/Analysis/scripts/
# promoter_orientation.R) with the rest of the curated promoter features
# (build_cooccurrence_features(), same feature set as R/09_coocurrencia/
# coocurrencia.R). Phi coefficient (Pearson correlation on 0/1 variables)
# between "Unidirectional" and each other feature.
#
# not_detected promoters (no PRO-seq signal) are NA for "Unidirectional"
# and excluded pairwise from each correlation (cor(..., use =
# "pairwise.complete.obs")).
#
# Requires: data/processed/prom_df.tsv (see R/00_prom_features) and
# ../transcriptional_library/External_data/PRO-seq/promoter_OI.tsv (see
# promoter_orientation.R, transcriptional_library repo). Run from the
# TesisDoc repo root.
#
# Exploratory - not part of the numbered pipeline (no docs/mapping_
# figuras.csv entry), same as top_repressors_panel.R.

library(tidyverse)
library(fastDummies)
source("R/functions/plot_helpers.R")

out_dir <- "figures/09_coocurrencia"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
oi_df <- read_tsv("../transcriptional_library/External_data/PRO-seq/promoter_OI.tsv", show_col_types = FALSE)

tidy_data <- build_cooccurrence_features(prom_df) %>%
  left_join(
    oi_df %>%
      mutate(Unidirectional = case_when(
        class == "unidirectional" ~ TRUE,
        class == "bidirectional" ~ FALSE,
        TRUE ~ NA
      )) %>%
      select(seq_id, Unidirectional),
    by = "seq_id"
  )

message("Promoters with orientation call: ", sum(!is.na(tidy_data$Unidirectional)), " / ", nrow(tidy_data))

feat_mat <- tidy_data %>%
  select(-seq_id) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

phi_vals <- cor(feat_mat[, "Unidirectional"], feat_mat[, colnames(feat_mat) != "Unidirectional"], use = "pairwise.complete.obs")[1, ]

phi_df <- tibble(feature = names(phi_vals), phi = phi_vals) %>%
  arrange(phi) %>%
  mutate(feature = fct_inorder(feature))

p <- ggplot(phi_df, aes(x = phi, y = feature, fill = phi)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0, limits = c(-1, 1), guide = "none") +
  labs(x = "Phi coefficient with \"Unidirectional\"", y = NULL) +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "cooccurrence_orientation_phi.jpg"), p, width = 8, height = 9, units = "in")
write_tsv(phi_df %>% arrange(desc(phi)), file.path(out_dir, "cooccurrence_orientation_phi.tsv"))

message("Figura y tabla guardadas en ", out_dir)
