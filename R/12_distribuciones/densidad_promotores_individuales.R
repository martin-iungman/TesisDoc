# Extra de distribuciones_expresion (no se usa en la tesis): densidad de
# fluorescencia EGFP (citometria) para cada uno de los 8 promotores
# validados individualmente, contra el control ("US"/Control) de fondo.
# Se guarda en la carpeta de distribuciones_expresion.
#
# Requiere: ~519MB de datos de citometria leidos directo de
# transcriptional_library (ver R/00_prom_features/heavy_data_paths.R).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")
source("R/functions/plot_helpers.R")
source("R/00_prom_features/heavy_data_paths.R")

slug <- "distribuciones_expresion"
out_dir <- fig_dir(slug)

df <- load_citometry_stable_validation(path_citometry_stable_validation)

# Mismos colores que distribuciones_expresion; sus 2 promotores de ejemplo
# (KIAA0753_1, TMEM87A_1) van primero y el control positivo al final.
cmv <- "CMV (control positivo)"
prom_colors <- c(
  KIAA0753_1 = "#AD343E", TMEM87A_1 = "#FFB400", BTG1_1 = "#7FB800", METAP2_1 = "#D6741F",
  PPP1R14B_3 = "#1B8C8E", ZKSCAN2_1 = "#9EDBCB", LSM1_1 = "#0D2C54", ETS1_1 = "#25afe9",
  setNames("black", cmv)
)
prom_order <- c("KIAA0753_1", "TMEM87A_1", sort(setdiff(names(prom_colors), c("KIAA0753_1", "TMEM87A_1", cmv))), cmv)

p <- df %>%
  filter(name != "Control") %>%
  mutate(name = factor(recode(name, Strong = cmv), levels = prom_order)) %>%
  ggplot() +
  geom_density(mapping = aes(Comp_FL2_A, col = name)) +
  scale_x_log10() +
  scale_color_manual(values = prom_colors) +
  facet_wrap(~name) +
  geom_density(data = df %>% filter(name == "Control") %>% select(-name), mapping = aes(Comp_FL2_A), col = "grey") +
  xlab("Señal de EGFP\n(unidades de fluorescencia relativa)") +
  ylab("Densidad") +
  theme_pubclean() +
  theme(legend.position = "none", axis.title = element_text(size = 16), plot.margin = margin(5.5, 15, 5.5, 5.5))
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("densidad_individual_promotores.", ext)), p, width = 9, height = 6.75, units = "in")
}

message("Figura guardada en ", out_dir)
