# Faceted comparison plot for remap_tf_maxtpm_hek_active_ora_thresholds.R:
# top GO terms (by p.adjust) for sobreestima/subestima x FDR<0.05/FDR<0.001.

library(tidyverse)
library(readxl)
library(ggpubr)

out_dir <- "figures/discordancia_reportero_endogeno/2_dif_signed_bland_altman"
xlsx_path <- file.path(out_dir, "ORA_remap_rank_dif_maxtpm_hek_active_thresholds.xlsx")

sheets <- c("sobreestima_fdr0.05", "sobreestima_fdr0.001", "subestima_fdr0.05", "subestima_fdr0.001")

read_top <- function(sheet, n_top = 12) {
  read_excel(xlsx_path, sheet = sheet) %>%
    arrange(p.adjust) %>%
    slice_head(n = n_top) %>%
    mutate(panel = sheet)
}

datos <- map_dfr(sheets, read_top) %>%
  separate(panel, into = c("direction", "fdr"), sep = "_fdr") %>%
  mutate(
    direction = recode(direction, "sobreestima" = "Reportero sobreestima", "subestima" = "Reportero subestima"),
    direction = factor(direction, levels = c("Reportero sobreestima", "Reportero subestima")),
    fdr = paste0("FDR<", fdr),
    label = paste(Description, fdr, direction, sep = "___")
  ) %>%
  group_by(direction, fdr) %>%
  mutate(Description = fct_reorder(paste0(Description, "  "), -p.adjust)) %>%
  ungroup()

p <- ggplot(datos, aes(x = -log10(p.adjust), y = Description, fill = Count)) +
  geom_col() +
  facet_wrap(direction ~ fdr, scales = "free_y", ncol = 2) +
  theme_pubr(base_size = 11) +
  theme(legend.position = "right", strip.text = element_text(face = "bold")) +
  labs(
    x = "-log10(p ajustado)", y = NULL, fill = "N° TFs",
    title = "ORA - max_tpm + todo ChIP, solo HEK-activos",
    subtitle = "Comparacion FDR<0.05 vs FDR<0.001 (top 12 terminos por panel)"
  ) +
  scale_fill_gradient(low = "#7FB800", high = "#0D2C54")

ggsave(file.path(out_dir, "ORA_terminos_comparacion_thresholds.jpg"), p, width = 16, height = 11, units = "in", dpi = 200)
message("Guardado en ", out_dir)
