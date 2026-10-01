# Endogenous HEK293 CAGE activity (hek_tpm, prom_df.tsv) vs. PRO-seq
# promoter orientation (Unidirectional/Bidirectional class, and the
# continuous Orientation Index - see ../transcriptional_library/
# Analysis/scripts/promoter_orientation.R). Both measurements come from
# HEK293 (matched cell line).
#
# not_detected promoters (no PRO-seq signal) are excluded from the
# violin comparison but shown as their own group would require a 3rd
# level - kept out for the same reason as coocurrencia_orientation.R.
#
# Requires: data/processed/prom_df.tsv (see R/00_prom_features) and
# ../transcriptional_library/External_data/PRO-seq/promoter_OI.tsv (see
# promoter_orientation.R, transcriptional_library repo). Run from the
# TesisDoc repo root.
#
# Exploratory - not part of the numbered pipeline (no docs/mapping_
# figuras.csv entry), same as top_repressors_panel.R.

library(tidyverse)
library(ggpubr)

out_dir <- "figures/09_coocurrencia"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>%
  select(seq_id, hek_tpm)
oi_df <- read_tsv("../transcriptional_library/External_data/PRO-seq/promoter_OI.tsv", show_col_types = FALSE)

data <- oi_df %>%
  select(seq_id, OI_mean, class) %>%
  inner_join(prom_df, by = "seq_id") %>%
  filter(!is.na(hek_tpm))

message("Promoters with hek_tpm + OI class: ", nrow(data))
print(table(data$class))

# --- Violin: hek_tpm by orientation class (unidirectional vs bidirectional) --

viol_data <- data %>%
  filter(class %in% c("unidirectional", "bidirectional")) %>%
  mutate(orientation = factor(
    ifelse(class == "unidirectional", "Unidirectional", "Bidirectional"),
    levels = c("Unidirectional", "Bidirectional")
  ))

orient_colors <- c(Unidirectional = "#4292C6", Bidirectional = "#E45C3A")

p_violin <- viol_data %>%
  ggviolin(x = "orientation", y = "hek_tpm", fill = "orientation",
           draw_quantiles = 0.5, alpha = 0.8) +
  scale_y_log10() +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = orient_colors) +
  labs(x = NULL, y = "Endogenous HEK293 activity (CAGE TPM)") +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "orientation_hek_tpm_violin.jpg"), p_violin, width = 6, height = 6, units = "in")

# --- Scatter: hek_tpm vs. Orientation Index (continuous) -------------------

p_scatter <- data %>%
  filter(class != "not_detected") %>%
  ggplot(aes(x = OI_mean, y = hek_tpm)) +
  geom_point(size = 0.6, alpha = 0.3, color = "#358AAA") +
  geom_smooth(method = "loess", color = "#216869") +
  scale_y_log10() +
  labs(x = "Orientation Index (mean reps)", y = "Endogenous HEK293 activity (CAGE TPM)") +
  theme_bw(base_size = 13)

cor_test <- cor.test(data$OI_mean[data$class != "not_detected"],
                     log10(data$hek_tpm[data$class != "not_detected"]),
                     method = "spearman")
message(sprintf("Spearman rho = %.3f, p = %.2e (n=%d)", cor_test$estimate, cor_test$p.value,
                sum(data$class != "not_detected")))

ggsave(file.path(out_dir, "orientation_hek_tpm_scatter.jpg"), p_scatter, width = 7, height = 6, units = "in")

message("Figuras guardadas en ", out_dir)
