# composicion_secuencia_library (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Panel A: motivos EPD (TATA, INR, CCAAT, GC-box)
# Panel B: islas CpG
# Panel C: contenido G+C
# Panel D: patrones de motivo en el TSS
#
# Requiere: data/processed/prom_df.tsv (generado por
# R/00_prom_features/build_prom_features.R). Run from the TesisDoc repo root.

library(tidyverse)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "composicion_secuencia_library"

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE)
promoters_df <- prom_df %>% filter(type == "promoter")
n_promoters <- nrow(promoters_df)

# --- Panel A: motivos EPD -------------------------------------------------

panel_a <- promoters_df %>%
  select(TATA_EPD, CCAAT_EPD, INR_EPD, GCbox_EPD) %>%
  pivot_longer(everything()) %>%
  mutate(name = str_remove(name, "_EPD")) %>%
  count(name, value) %>%
  filter(value) %>%
  count_bar_labeled(name, n, n_promoters) +
  labs(y = "Frecuencia", x = "Motivo")

# --- Panel B: islas CpG ----------------------------------------------------

panel_b <- promoters_df %>%
  count(CGI) %>%
  filter(CGI) %>%
  mutate(CGI = "Isla CpG") %>%
  count_bar_labeled(CGI, n, n_promoters) +
  labs(y = "Frecuencia", x = "") +
  ylim(c(0, n_promoters))

# --- Panel C: contenido G+C --------------------------------------------

panel_c <- promoters_df %>%
  ggplot(aes(g_c * 100)) +
  geom_density(col = thesis_clr, linewidth = 1.2, fill = thesis_clr, alpha = 0.6) +
  xlab("Contenido [G+C] (%)") +
  ylab("Densidad") +
  theme_pubclean(base_size = 16)

# --- Panel D: patrones de motivo en el TSS ------------------------------

panel_d <- promoters_df %>%
  select(ends_with("TSS")) %>%
  pivot_longer(ends_with("TSS")) %>%
  mutate(
    name = ifelse(name == "other_TSS", "No\ncanónico", ifelse(name == "INR_strong_TSS", "INR\nfuerte", name)) %>%
      str_remove("_TSS")
  ) %>%
  count(name, value) %>%
  filter(value) %>%
  arrange(desc(n)) %>%
  mutate(name = fct_inorder(name)) %>%
  count_bar_labeled(name, n, n_promoters) +
  labs(y = "Frecuencia", x = "Motivo en TSS")

# --- Save individual panels + combined figure -----------------------------

out_dir <- fig_dir(slug)
ggsave(file.path(out_dir, "panel_a_motivos_epd.jpg"), panel_a, width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "panel_b_cgi.jpg"), panel_b, width = 3, height = 6.75, units = "in")
ggsave(file.path(out_dir, "panel_c_gc.jpg"), panel_c, width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "panel_d_tss_motifs.jpg"), panel_d, width = 9, height = 6.75, units = "in")

combined <- (panel_a + panel_b) / (panel_c + panel_d) + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, paste0(slug, ".jpg")), combined, width = 14, height = 12, units = "in")
