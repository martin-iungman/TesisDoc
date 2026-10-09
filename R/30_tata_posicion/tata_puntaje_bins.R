# tata_puntaje (ver docs/mapping_figuras.csv para el numero de figura
# vigente; los otros paneles de la figura los genera
# R/30_tata_posicion/tata_posicion.R): actividad segun el puntaje de la
# TATA-box canonica, por bins FIJOS de puntaje (misma idea que
# tata_puntaje_dosis.jpg, pero con los grupos armados por valor de puntaje y
# no por cantidad de promotores; mismo formato que cgi_oe_deciles de
# R/14_tata_cgi/tata_cgi.R).
#
# Para cada promotor, el mejor puntaje (pesos W de EPD, misma hebra) con la
# columna de referencia entre -31 y -26 del TSS, sin corte; y lo mismo en una
# ventana control del mismo ancho en -131 a -126. Bins fijos, no acumulados
# (cada promotor cae en uno solo, iguales en ambas replicas y ventanas):
# < -16, -16 a -14, -14 a -12, -12 a -10, -10 a -8.17, -8.17 a -7, -7 a -6,
# -6 a -5, >= -5. Eje y: percentil de actividad medio del bin (dentro de cada
# replica), IC 95%. Eje x: mediana del puntaje del bin. En el grafico se
# omiten los bins con menos de N_MIN promotores (los de puntaje alto de la
# ventana control); la tabla los incluye a todos con su n.
# Referencias: vertical en el corte de EPD (-8.17), horizontal en 50.
#
# Salida en fig_dir("tata_puntaje"): tata_puntaje_bins.jpg y
# tata_puntaje_bins.tsv.
#
# Requiere: data/raw/library.bed, BSgenome.Hsapiens.UCSC.hg38,
# data/processed/activity_stats_highconf.tsv y data/processed/prom_df.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
source("R/functions/fig_paths.R")
source("R/functions/tata_scan.R")

select <- dplyr::select
filter <- dplyr::filter

out_dir <- fig_dir("tata_puntaje")
CORTE <- TATA_CORTE_EPD
CONTROL <- TATA_CANON - 100
LIMITES <- c(-Inf, -16, -14, -12, -10, CORTE, -7, -6, -5, Inf)
ETIQUETAS <- c("< −16", "−16 a −14", "−14 a −12", "−12 a −10", "−10 a −8,17", "−8,17 a −7", "−7 a −6", "−6 a −5", "≥ −5")
N_MIN <- 20
VENTANAS <- c("Canónica (−31 a −26)", "Control (−131 a −126)")

# --- Mejor puntaje por promotor, en la ventana canonica y en la control -----

lib <- tata_load_library()
code <- lib$code[!lib$has_n, ] # 1 fragmento con N fuera
scores_misma <- tata_score_matrix(TATA_W, code)
best <- tibble(
  seq_id = rownames(code),
  canonica = tata_best_in_window(scores_misma, TATA_CANON[1], TATA_CANON[2]),
  control = tata_best_in_window(scores_misma, CONTROL[1], CONTROL[2])
) %>%
  pivot_longer(-seq_id, names_to = "ventana", values_to = "score") %>%
  mutate(ventana = factor(if_else(ventana == "canonica", VENTANAS[1], VENTANAS[2]), levels = VENTANAS))

# Actividad media por replica, como en tata_cgi.R; percentil dentro de cada replica
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  filter(seq_id %in% rownames(code)) %>%
  group_by(rep) %>%
  mutate(act_pct = 100 * percent_rank(mean)) %>%
  ungroup() %>%
  inner_join(best, by = "seq_id", relationship = "many-to-many") %>%
  # intervalos cerrados a izquierda: [a, b); se suma la tolerancia para que
  # un puntaje igual al limite (p. ej. -8.17) caiga en el bin superior
  mutate(bin = cut(score + TATA_TOL, breaks = LIMITES, labels = ETIQUETAS, right = FALSE))

# --- Resumen por bin ------------------------------------------------------

por_bin <- data %>%
  group_by(rep, ventana, bin) %>%
  summarise(
    n = n(), puntaje_mediana = median(score), puntaje_min = min(score), puntaje_max = max(score),
    act_pct_media = mean(act_pct), se = sd(act_pct) / sqrt(n()), .groups = "drop"
  ) %>%
  mutate(conf.low = act_pct_media - 1.96 * se, conf.high = act_pct_media + 1.96 * se)
write_tsv(por_bin, file.path(out_dir, "tata_puntaje_bins.tsv"))
options(width = 200)
print(por_bin %>% mutate(across(where(is.double), ~ round(.x, 2))), n = Inf, width = Inf)

# --- Panel ------------------------------------------------------------------

panel <- ggplot(filter(por_bin, n >= N_MIN), aes(puntaje_mediana, act_pct_media, color = ventana)) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = CORTE, linetype = "dashed", color = "grey40") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0) +
  geom_line() +
  geom_point(size = 3) +
  annotate("text", x = CORTE, y = Inf, label = "−8,17 (EPD)", hjust = 1.05, vjust = 1.4, size = 6, color = "grey30") +
  scale_color_manual(values = setNames(c("#216869", "grey55"), VENTANAS)) +
  scale_x_continuous(labels = function(x) sub("^-", "−", x), expand = expansion(add = c(1, 0.8))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  theme_pubclean() +
  theme(text = element_text(size = 30), panel.spacing = unit(3, "lines"), legend.text = element_text(size = 22)) +
  facet_wrap(~rep) +
  labs(x = "Mejor puntaje TATA en la ventana (mediana del bin)", y = "Percentil de\nactividad media", color = NULL) +
  ggtitle("TATA-box canónica según su puntaje")
ggsave(file.path(out_dir, "tata_puntaje_bins.jpg"), panel, width = 16, height = 7.5, units = "in")
