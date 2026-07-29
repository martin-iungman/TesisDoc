# distribuciones_expresion (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Validacion por citometria de 2 promotores de ejemplo (KIAA0753_1,
# TMEM87A_1) contra la actividad estimada por el ensayo high-throughput
# (library) y por el conteo de gates de sorting.
# Panel A ("aka 1C"): densidad de fluorescencia EGFP (citometria) para los
#   2 promotores de ejemplo.
# Panel B: histograma de reconstruccion de distribucion (cuentas
#   normalizadas por gate de sorting), ambos promotores, solo Rep 2 -
#   mismo dato que Fig. R1.4 (histograma_gates_kiaa0753) pero con
#   TMEM87A_1 tambien, para comparar forma de distribucion entre
#   promotores.
# Panel C ("aka 1D"): correlacion entre la media del ensayo high-throughput
#   y la media del ensayo especifico (citometria) por promotor.
#
# Resync 2026-07-29: separado en 3 figuras (docs/Fig R1.pptx paso a tener
# slides propios para cada una) - el histograma de KIAA0753_1 solo (Fig.
# R1.4) y la densidad individual de los 8 promotores (Fig. R1.5) ahora
# viven en R/12_distribuciones/histograma_gates_kiaa0753.R y
# densidad_promotores_individuales.R respectivamente; este script se quedo
# con los paneles A/B/C y paso de "Fig. R1.4" a "Fig. R1.6".
#
# Requiere: data/processed/activity_stats_full.tsv (generado por
# R/01_activity_stats/build_activity_stats.R), data/processed/data_long.tsv
# y ~519MB de datos de citometria leidos directo de transcriptional_library
# (ver R/00_prom_features/heavy_data_paths.R). Run from the TesisDoc repo
# root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "distribuciones_expresion"
out_dir <- fig_dir(slug)

df <- load_citometry_stable_validation(path_citometry_stable_validation)

# --- Panel A: densidad de EGFP, 2 promotores de ejemplo -------------------

panel_a <- df %>%
  filter(name %in% c("KIAA0753_1", "TMEM87A_1")) %>%
  ggplot(aes(Comp_FL2_A, fill = name, group = name)) +
  geom_density(alpha = 0.6) +
  scale_x_log10(limits = c(1, 1000)) +
  theme_pubr() +
  labs(x = "Señal de fluorescencia de EGFP", y = "Densidad", fill = "Promotor", alpha = NULL) +
  scale_fill_manual(values = c(KIAA0753_1 = "#AD343E", TMEM87A_1 = "#FFB400"))

ggsave(file.path(out_dir, "panel_a_densidad_ejemplo.jpg"), panel_a, width = 9, height = 6.75, units = "in")

# --- Panel B: histograma de reconstruccion, ambos promotores, solo Rep 2 --

gates_hist <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE) %>%
  filter(name %in% c("KIAA0753_1", "TMEM87A_1"), rep == "Rep 2") %>%
  group_by(rep, name) %>%
  mutate(counts_rel = counts_norm / sum(counts_norm)) %>%
  ungroup()

panel_b <- gates_hist %>%
  ggplot(aes(factor(sample), counts_rel, fill = name)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(KIAA0753_1 = "#AD343E", TMEM87A_1 = "#FFB400")) +
  labs(x = "Gate de fluorescencia EGFP", y = "Cuentas relativas", fill = "Promotor") +
  theme_pubclean(base_size = 20)

ggsave(file.path(out_dir, "panel_b_histograma_gates_ejemplo.jpg"), panel_b, width = 9, height = 6.75, units = "in")

# --- Panel C: correlacion library vs citometria (media) -------------------

act_mean <- read_tsv("data/processed/activity_stats_full.tsv", show_col_types = FALSE) %>%
  group_by(seq_id, name) %>%
  summarise(
    mean_hist = mean(mean), sd_mean_hist = sd(mean) %>% replace_na(0),
    var_hist = mean(var), sd_var_hist = sd(var) %>% replace_na(0),
    cv_hist = mean(sqrt(var) / mean), sd_cv_hist = sd(sqrt(var) / mean) %>% replace_na(0),
    fano_hist = mean(var / mean) %>% replace_na(0), sd_fano_hist = sd(var / mean) %>% replace_na(0),
    .groups = "drop"
  )

density_mean <- df %>%
  group_by(sample_name, name) %>%
  summarise(mean_EGFP = mean(log10(Comp_FL2_A), na.rm = TRUE), .groups = "drop") %>%
  group_by(name) %>%
  summarise(mean_density = mean(mean_EGFP), sd_density = sd(mean_EGFP), .groups = "drop")

panel_c <- density_mean %>%
  left_join(act_mean, by = "name") %>%
  ggplot(aes(mean_hist, mean_density, col = name)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = mean_density - sd_density, ymax = mean_density + sd_density), width = .125) +
  geom_errorbarh(aes(xmin = mean_hist - sd_mean_hist, xmax = mean_hist + sd_mean_hist), height = .05) +
  labs(x = "Media (Ensayo high-throughput)", y = "Media (Ensayo especifico)", col = "Promotor") +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = c(
    KIAA0753_1 = "#AD343E", TMEM87A_1 = "#FFB400", BTG1_1 = "#7FB800", METAP2_1 = "#D6741F",
    PPP1R14B_3 = "#1B8C8E", ZKSCAN2_1 = "#9EDBCB", LSM1_1 = "#0D2C54", ETS1_1 = "#25afe9"
  )) +
  ylim(c(1.2, 1.8)) +
  ggtitle("Validación de actividad")

ggsave(file.path(out_dir, "panel_c_correlacion_media.jpg"), panel_c, width = 9, height = 6.75, units = "in")

message("Paneles A, B y C guardados en ", out_dir, ".")
