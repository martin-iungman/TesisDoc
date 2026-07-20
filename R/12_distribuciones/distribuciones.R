# distribuciones_expresion (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Validacion por citometria de 8 promotores individuales + 2 controles
# (Control="US", Strong) contra la actividad estimada por el ensayo
# high-throughput (library).
# Panel B ("aka 1C"): densidad de fluorescencia EGFP para 2 promotores de
#   ejemplo (KIAA0753_1, TMEM87A_1).
# Panel C ("aka 1D"): correlacion entre la media del ensayo high-throughput
#   y la media del ensayo especifico (citometria) por promotor.
# Extra (sin mapear a panel): densidad de EGFP para cada uno de los 8
#   promotores, contra el control ("US"/Control) de fondo.
#
# PENDIENTE: "Panel A: histogramas de reconstruccion de distribuciones
# (aka 1B)" mencionado en mapping_figuras.csv no esta en
# transcriptional_library/Tesis/tesis.R - el codigo debe estar en otro
# script, buscarlo mas adelante.
#
# Requiere: data/processed/activity_stats_full.tsv (generado por
# R/01_activity_stats/build_activity_stats.R) y ~519MB de datos de
# citometria leidos directo de transcriptional_library (ver
# R/00_prom_features/heavy_data_paths.R). Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "distribuciones_expresion"
out_dir <- fig_dir(slug)

# --- Cargar datos de citometria (FlowJo export, un .cells.csv por muestra) --

files <- list.files(path_citometry_stable_validation, pattern = ".cells", full.names = TRUE, recursive = TRUE)
prom_name <- files %>%
  str_remove("^.+Tables/(lvs-)?") %>%
  str_remove("_Data .+$")
df <- map2(files, prom_name, ~ read_csv(.x, show_col_types = FALSE) %>% mutate(sample_name = .y)) %>% list_rbind()
names(df) <- str_replace_all(names(df), "-", "_")
df$prom_name <- str_remove(df$sample_name, " -.+$")

# El orden de este vector tiene que coincidir con el orden en que
# unique(df$prom_name) devuelve los short_name - fragil pero heredado tal
# cual de transcriptional_library/Tesis/tesis.R.
name_df <- tibble(
  short_name = unique(df$prom_name),
  name = c("BTG1_1", "ETS1_1", "KIAA0753_1", "LSM1_1", "METAP2_1", "PPP1R14B_3", "TMEM87A_1", "ZKSCAN2_1", "Control", "Strong")
)
df <- inner_join(name_df, df, by = c("short_name" = "prom_name"))

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

# --- Panel B: densidad de EGFP, 2 promotores de ejemplo -------------------

panel_b <- df %>%
  filter(name %in% c("KIAA0753_1", "TMEM87A_1")) %>%
  ggplot(aes(Comp_FL2_A, fill = name, group = name)) +
  geom_density(alpha = 0.6) +
  scale_x_log10(limits = c(1, 1000)) +
  theme_pubr() +
  labs(x = "Señal de fluorescencia de EGFP", y = "Densidad", fill = "Promotor", alpha = NULL) +
  scale_fill_manual(values = c(KIAA0753_1 = "#AD343E", TMEM87A_1 = "#FFB400"))

ggsave(file.path(out_dir, "panel_b_densidad_ejemplo.jpg"), panel_b, width = 9, height = 6.75, units = "in")

# --- Extra: densidad de EGFP para los 8 promotores individuales -----------

panel_extra <- df %>%
  filter(name != "Control") %>%
  ggplot() +
  geom_density(mapping = aes(Comp_FL2_A, col = name)) +
  scale_x_log10() +
  facet_wrap(~name) +
  geom_density(data = df %>% filter(name == "Control") %>% select(-name), mapping = aes(Comp_FL2_A), col = "grey") +
  xlab("Señal de EGFP\n(unidades de fluorescencia relativa)") +
  ylab("Densidad") +
  theme_pubclean() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "extra_densidad_individual.jpg"), panel_extra, width = 9, height = 6.75, units = "in")

message("Paneles B, C y el extra guardados en ", out_dir, ". Panel A (histogramas de reconstrucción) sigue pendiente.")
