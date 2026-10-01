# English, editable-text figure for the paper: ORA terms for BOTH the
# "reporter overestimates" and "reporter underestimates" sets @
# FDR<0.001 (max_tpm + all ReMap ChIP, HEK293-active promoters only),
# as two panels of the same plot with an aligned (shared) x-axis.
# Supersedes the separate ora_sobreestima_fdr001_english.R /
# ora_subestima_fdr001_english.R single-panel figures.
#
# Output goes to ../transcriptional_library/Plots/Review/.

library(tidyverse)
library(readxl)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
plots_review_dir <- "../transcriptional_library/Plots/Review"
xlsx_path <- file.path(out_dir, "ORA_remap_rank_dif_maxtpm_hek_active_thresholds.xlsx")

read_top <- function(sheet, label, n_top = 12) {
  read_excel(xlsx_path, sheet = sheet) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top) %>%
    mutate(direction = label)
}

datos <- bind_rows(
  read_top("sobreestima_fdr0.001", "Reporter overestimates"),
  read_top("subestima_fdr0.001", "Reporter underestimates")
) %>%
  mutate(
    direction = factor(direction, levels = c("Reporter overestimates", "Reporter underestimates")),
    neglog10_p = -log10(p.adjust)
  ) %>%
  mutate(Description = str_wrap(Description, width = 30)) %>%
  group_by(direction) %>%
  mutate(Description = fct_reorder(paste0(Description, "  "), neglog10_p)) %>%
  ungroup()

p <- ggplot(datos, aes(x = neglog10_p, y = Description, fill = Count)) +
  geom_col() +
  facet_wrap(~direction, ncol = 1, scales = "free_y") +
  scale_fill_gradient(low = "#7FB800", high = "#0D2C54") +
  labs(x = "-log10(adjusted p-value)", y = NULL, fill = "Number of TFs") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = rel(0.9), lineheight = 0.85),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )

dir.create(plots_review_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plots_review_dir, "ORA_terms_maxtpm_hek_active_combined_english.jpg"), p, width = 11, height = 9.5, units = "in")
ggsave(file.path(plots_review_dir, "ORA_terms_maxtpm_hek_active_combined_english.pdf"), p, width = 11, height = 9.5, units = "in", device = cairo_pdf)

message("Figures saved to ", plots_review_dir)
