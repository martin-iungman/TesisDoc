# metrica_ruido (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Construccion de la metrica de ruido desacoplada de la actividad
# (var_rank_sw: rango de varianza dentro de cada bin de 100 promotores
# ordenados por actividad media, mean_sw) y su reproducibilidad entre
# replicas. Abre la seccion de ruido; la dependencia media-varianza que
# motiva el desacople esta en correccion_sesgo_muestreo (panel B).
#
# A: construccion (Rep 1) - scatter media vs. varianza con los bins de
#    actividad en bandas alternadas y un bin intermedio resaltado; a la
#    derecha, los 100 promotores de ese bin ordenados por varianza con el
#    rango asignado (1-100).
# B: resultado (Rep 1) - var_rank_sw vs. bin de actividad (formato de
#    noise_scatter() sin colorear por feature) con mediana por bin y
#    loess, plano alrededor de 50.
# C: reproducibilidad - rango de varianza medio en Rep 2 (IC95) por decil
#    del rango en Rep 1, por tercil de actividad y para todos.
#
# var_rank_sw se calcula con add_noise_rank() sobre el set completo de
# cada replica por separado (igual que el resto de las figuras de ruido);
# los paneles/estadisticos entre replicas usan los promotores presentes
# en ambas (inner join por seq_id).
# Actividad para terciles/quintiles = promedio del percentil de mean en
# ambas replicas (percent_rank dentro del set unido).
#
# Requiere: data/processed/activity_stats_highconf.tsv (ver
# R/01_activity_stats). Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "metrica_ruido"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
data <- add_noise_rank(stats_highconf)

# bin de actividad resaltado en A y B: el bin mediano de Rep 1
rep1 <- data %>% filter(rep == "Rep 1")
bin_sel <- ceiling(max(rep1$mean_sw) / 2)
clr_sel <- "#AD343E"

# --- Panel A: construccion ----------------------------------------------

bandas <- rep1 %>%
  group_by(mean_sw) %>%
  summarise(xmin = min(mean), xmax = max(mean), .groups = "drop") %>%
  # los bins comparten bordes cuando hay empates en mean; se extiende cada
  # banda hasta el punto medio con la siguiente para que no queden huecos
  mutate(
    xmax = (xmax + lead(xmin, default = last(xmax))) / 2,
    xmin = lag(xmax, default = first(xmin)),
    banda = ifelse(mean_sw %% 2 == 0, "par", "impar")
  )

panel_a_izq <- ggplot() +
  geom_rect(
    data = bandas, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = banda),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(par = "white", impar = "#cfe3ef")) +
  geom_rect(
    data = filter(bandas, mean_sw == bin_sel),
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = clr_sel, alpha = 0.25
  ) +
  geom_point(data = filter(rep1, mean_sw != bin_sel), aes(mean, var), size = 0.3, col = "grey55") +
  geom_point(data = filter(rep1, mean_sw == bin_sel), aes(mean, var), size = 0.8, col = clr_sel) +
  labs(
    x = "Media", y = "Varianza",
    title = paste0("Bins de 100 promotores por actividad\n(bin ", bin_sel, " resaltado)")
  ) +
  theme_pubr(base_size = 16)

bin_df <- rep1 %>% filter(mean_sw == bin_sel)

panel_a_der <- bin_df %>%
  ggplot(aes(var, var_rank_sw, col = var_rank_sw)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "#f2c4c8", high = clr_sel, guide = "none") +
  scale_y_continuous(breaks = c(1, 25, 50, 75, 100)) +
  annotate(
    "segment", x = min(bin_df$var), xend = min(bin_df$var), y = 15, yend = 90,
    arrow = arrow(length = unit(0.25, "cm")), col = "grey30"
  ) +
  annotate(
    "text", x = min(bin_df$var), y = 52, label = "rango de varianza\ndentro del bin",
    hjust = -0.08, col = "grey30", size = 4.5
  ) +
  labs(x = "Varianza", y = "Rango de varianza (1-100)", title = paste0("Bin ", bin_sel, " ordenado por varianza")) +
  theme_pubr(base_size = 16)

panel_a <- panel_a_izq | panel_a_der

ggsave(file.path(out_dir, "panel_a_construccion.jpg"), panel_a, width = 13, height = 6, units = "in")

# --- Panel B: desacople ----------------------------------------------------

mediana_bin <- rep1 %>%
  group_by(mean_sw) %>%
  summarise(var_rank_med = median(var_rank_sw), .groups = "drop")

panel_b <- rep1 %>%
  ggplot(aes(mean_sw, var_rank_sw)) +
  geom_point(size = 0.2, alpha = 0.5, col = "grey") +
  geom_point(data = bin_df, size = 0.8, col = clr_sel) +
  geom_hline(yintercept = 50, linetype = "dashed", col = "grey40") +
  geom_line(data = mediana_bin, aes(y = var_rank_med), col = thesis_clr, linewidth = 0.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2, col = "#136869") +
  labs(x = "Bin de actividad media", y = "Rango de varianza (dentro del bin)", title = "Ruido desacoplado de la actividad") +
  theme_pubr(base_size = 16)

ggsave(file.path(out_dir, "panel_b_desacople.jpg"), panel_b, width = 9, height = 6.75, units = "in")

# --- Datos unidos por seq_id ------------------------------------------------

wide <- data %>%
  select(seq_id, rep, mean, var, var_rank_sw) %>%
  mutate(rep = recode(rep, "Rep 1" = "r1", "Rep 2" = "r2")) %>%
  pivot_wider(names_from = rep, values_from = c(mean, var, var_rank_sw)) %>%
  drop_na() %>%
  mutate(
    actividad = (percent_rank(mean_r1) + percent_rank(mean_r2)) / 2,
    tercil = ntile(actividad, 3),
    quintil = ntile(actividad, 5)
  )
message(nrow(wide), " promotores en ambas replicas")

# --- Panel C: reproducibilidad ----------------------------------------------

tercil_lab <- c("1" = "Baja", "2" = "Media", "3" = "Alta")

# rango medio (IC95, t) en la replica `y` por decil del rango en la replica `x`
decil_summary <- function(df, x, y) {
  df <- df %>% mutate(decil = ceiling(.data[[x]] / 10), valor = .data[[y]])
  bind_rows(
    df %>% mutate(grupo = tercil_lab[as.character(tercil)]),
    df %>% mutate(grupo = "Todos")
  ) %>%
    group_by(grupo, decil) %>%
    summarise(
      n = n(), media = mean(valor), sd = sd(valor), .groups = "drop"
    ) %>%
    mutate(
      ic95_inf = media - qt(0.975, n - 1) * sd / sqrt(n),
      ic95_sup = media + qt(0.975, n - 1) * sd / sqrt(n)
    ) %>%
    select(-sd)
}

panel_c_df <- bind_rows(
  decil_summary(wide, "var_rank_sw_r1", "var_rank_sw_r2") %>% mutate(direccion = "Rep 1 -> Rep 2"),
  decil_summary(wide, "var_rank_sw_r2", "var_rank_sw_r1") %>% mutate(direccion = "Rep 2 -> Rep 1")
) %>%
  mutate(grupo = factor(grupo, levels = c("Baja", "Media", "Alta", "Todos"))) %>%
  relocate(direccion)

write_tsv(panel_c_df, file.path(out_dir, "panel_c_reproducibilidad.tsv"))

panel_c <- panel_c_df %>%
  filter(direccion == "Rep 1 -> Rep 2") %>%
  ggplot(aes(decil, media, col = grupo)) +
  geom_hline(yintercept = 50, linetype = "dashed", col = "grey40") +
  geom_errorbar(aes(ymin = ic95_inf, ymax = ic95_sup), width = 0.25, linewidth = 0.6) +
  geom_line(aes(linewidth = grupo == "Todos")) +
  geom_point(size = 2.2) +
  scale_linewidth_manual(values = c(`FALSE` = 0.8, `TRUE` = 1.4), guide = "none") +
  scale_color_manual(values = c(Baja = "#9ECAE1", Media = "#4292C6", Alta = "#08306B", Todos = "grey20")) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    x = "Decil del rango de varianza en la réplica 1",
    y = "Rango de varianza medio\nen la réplica 2 (IC95)",
    col = "Actividad", title = "Reproducibilidad entre réplicas"
  ) +
  theme_pubr(base_size = 16)

ggsave(file.path(out_dir, "panel_c_reproducibilidad.jpg"), panel_c, width = 9, height = 6.75, units = "in")

# --- Figura compuesta ---------------------------------------------------------

compuesta <- wrap_elements(full = panel_a) / (panel_b | panel_c) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(file.path(out_dir, "metrica_ruido.jpg"), compuesta, width = 16, height = 13, units = "in")

# --- Estadisticos -------------------------------------------------------------

cor_entre <- map(c("mean", "var", "var_rank_sw"), function(v) {
  x <- wide[[paste0(v, "_r1")]]
  y <- wide[[paste0(v, "_r2")]]
  tibble(
    categoria = "entre_replicas", estadistico = c("pearson", "spearman"), variable = v,
    subgrupo = "todos", valor = c(cor(x, y), cor(x, y, method = "spearman")), n = length(x)
  )
}) %>% list_rbind()

cor_dentro <- data %>%
  group_by(rep) %>%
  summarise(
    `mean~var` = cor(mean, var, method = "spearman"),
    `mean~var_rank_sw` = cor(mean, var_rank_sw, method = "spearman"),
    n = n(), .groups = "drop"
  ) %>%
  pivot_longer(c(`mean~var`, `mean~var_rank_sw`), names_to = "variable", values_to = "valor") %>%
  transmute(categoria = "dentro_replica", estadistico = "spearman", variable, subgrupo = rep, valor, n)

# P(rango en `y` >= 75 o <= 25 | rango en `x` >= 75)
concordancia <- function(x, y, dir) {
  alto <- wide[[x]] >= 75
  tibble(
    categoria = "concordancia",
    estadistico = c("P(>=75 | >=75)", "P(<=25 | >=75)"),
    variable = "var_rank_sw", subgrupo = dir,
    valor = c(mean(wide[[y]][alto] >= 75), mean(wide[[y]][alto] <= 25)),
    n = sum(alto)
  )
}

cor_por_grupo <- function(col, nombre) {
  wide %>%
    group_by(g = .data[[col]]) %>%
    summarise(valor = cor(var_rank_sw_r1, var_rank_sw_r2, method = "spearman"), n = n(), .groups = "drop") %>%
    transmute(
      categoria = paste0("entre_replicas_por_", nombre), estadistico = "spearman",
      variable = "var_rank_sw", subgrupo = paste0(nombre, " ", g), valor, n
    )
}

estadisticos <- bind_rows(
  cor_entre,
  cor_dentro,
  concordancia("var_rank_sw_r1", "var_rank_sw_r2", "Rep 1 -> Rep 2"),
  concordancia("var_rank_sw_r2", "var_rank_sw_r1", "Rep 2 -> Rep 1"),
  cor_por_grupo("quintil", "quintil"),
  cor_por_grupo("tercil", "tercil")
)

write_tsv(estadisticos, file.path(out_dir, "metrica_ruido_estadisticos.tsv"))
print(estadisticos, n = Inf)

message("Figuras guardadas en ", out_dir)
