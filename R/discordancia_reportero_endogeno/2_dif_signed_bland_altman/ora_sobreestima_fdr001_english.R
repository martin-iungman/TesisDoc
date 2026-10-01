# English, editable-text figure for the paper: ORA terms for the
# "reporter overestimates" (sobreestima) set @ FDR<0.001, from
# remap_tf_maxtpm_hek_active_ora_thresholds.R (max_tpm + all ReMap ChIP,
# HEK293-active promoters only). No title/subtitle/threshold text on the
# plot itself, per convention for paper-ready figures. Companion to
# ora_subestima_fdr001_english.R.
#
# Output goes to ../transcriptional_library/Plots/Review/.

library(tidyverse)
library(readxl)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
plots_review_dir <- "../transcriptional_library/Plots/Review"

datos <- read_excel(file.path(out_dir, "ORA_remap_rank_dif_maxtpm_hek_active_thresholds.xlsx"), sheet = "sobreestima_fdr0.001") %>%
  arrange(p.adjust) %>%
  slice_head(n = 12) %>%
  mutate(Description = fct_reorder(Description, -p.adjust))

p <- ggplot(datos, aes(x = -log10(p.adjust), y = Description, fill = Count)) +
  geom_col() +
  scale_fill_gradient(low = "#7FB800", high = "#0D2C54") +
  labs(x = "-log10(adjusted p-value)", y = NULL, fill = "Number of TFs") +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())

dir.create(plots_review_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plots_review_dir, "ORA_overestimate_terms_maxtpm_hek_active_english.jpg"), p, width = 8, height = 6, units = "in")
ggsave(file.path(plots_review_dir, "ORA_overestimate_terms_maxtpm_hek_active_english.pdf"), p, width = 8, height = 6, units = "in", device = cairo_pdf)

message("Figures saved to ", plots_review_dir)
