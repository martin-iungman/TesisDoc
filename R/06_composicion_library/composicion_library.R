# composicion_secuencia_library (ver docs/mapping_figuras.csv para el numero
# de figura vigente)
# Panel A: motivos EPD (TATA, INR, CCAAT, GC-box)
# Panel B: islas CpG (definidas por la composicion del fragmento, CGI_frag)
# Panel C: contenido G+C
# Panel C2: razon CpG observado/esperado del fragmento (panel_cpg_oe.jpg)
# Panel D: patrones de motivo en el TSS
# Las letras de arriba (y de los nombres de variable/archivo) son historicas;
# en la compuesta el orden es A EPD, B CpG o/e, C islas CpG, D G+C, E TSS.
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
# Usa CGI_frag (el fragmento de 252 pb cumple G+C >= 50% y CpG o/e > 0.6),
# no la anotacion genomica de UCSC (CGI); ver R/00_prom_features/build_prom_features.R.

panel_b <- promoters_df %>%
  count(CGI = CGI_frag) %>%
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

# --- Panel C2: razon CpG o/e ----------------------------------------------
# cpg_oe del fragmento de 252 pb (ver R/00_prom_features/build_prom_features.R);
# la linea marca el corte de 0.6 de CGI_frag.

panel_cpg_oe <- promoters_df %>%
  filter(!is.na(cpg_oe)) %>%
  ggplot(aes(cpg_oe)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.025, fill = thesis_clr_bg, colour = "white", linewidth = 0.2) +
  geom_density(col = thesis_clr, linewidth = 1.2) +
  geom_vline(xintercept = 0.6, linetype = "dashed", colour = "grey30", linewidth = 0.8) +
  annotate("text", x = 0.62, y = Inf, vjust = 1.5, hjust = 0, label = "corte isla CpG", colour = "grey30", size = 5) +
  xlab("Razón CpG observado/esperado") +
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
ggsave(file.path(out_dir, "panel_cpg_oe.jpg"), panel_cpg_oe, width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "panel_d_tss_motifs.jpg"), panel_d, width = 9, height = 6.75, units = "in")

combined <- (panel_a + panel_cpg_oe) / (panel_b + panel_c + plot_layout(widths = c(1, 3))) / panel_d + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, paste0(slug, ".jpg")), combined, width = 14, height = 18, units = "in")
