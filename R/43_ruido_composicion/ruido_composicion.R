# ruido_composicion (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Ruido transcripcional (var_rank_sw, rango de varianza dentro de cada bin de
# actividad media) por decil de la composicion del fragmento de 252 pb:
# razon CpG o/e (cpg_oe, arriba) y contenido G+C (g_c, abajo), por replica.
# Mismo formato que el panel de deciles de CpG o/e de cgi_actividad
# (cgi_oe_deciles.jpg): punto = media del decil, barra = IC 95%, rho de
# Spearman por replica; linea discontinua vertical = umbral de CGI_frag
# (CpG o/e > 0.6, G+C >= 50%); linea punteada horizontal = 50 (~media nula
# del rango, 50.5; mismo valor que en cgi_oe_deciles).
# Tabla de soporte del texto: ruido_composicion_correlaciones.tsv (Spearman
# del ruido con cpg_oe y g_c, en total y dentro de cada grupo de CGI_frag,
# por replica - mismo formato que cgi_correlaciones.tsv de cgi_actividad).
# Analisis completo (bins finos, modelos, dependencia de la actividad) en el
# exploratorio R/42_exploratorio_ruido_composicion.
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y R/01_activity_stats).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "ruido_composicion"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  mutate(CGI_frag = as.logical(as.character(CGI_frag)))

# var_rank_sw sobre el set completo de cada replica (igual que ruido_cgi)
data <- add_noise_rank(data) %>% filter(!is.na(cpg_oe), !is.na(g_c))

# media e IC 95% del ruido por decil de `var`, por replica
deciles_ruido <- function(var) {
  data %>%
    group_by(rep) %>%
    mutate(decil = ntile(.data[[var]], 10)) %>%
    group_by(rep, decil) %>%
    summarise(x = median(.data[[var]]), m = mean(var_rank_sw), se = sd(var_rank_sw) / sqrt(n()), .groups = "drop") %>%
    mutate(ymin = m - 1.96 * se, ymax = m + 1.96 * se)
}
dec_cpg <- deciles_ruido("cpg_oe")
dec_gc <- deciles_ruido("g_c")
# misma escala en Y para los dos paneles (rango de los IC de ambos, con 50)
y_lims <- range(c(dec_cpg$ymin, dec_cpg$ymax, dec_gc$ymin, dec_gc$ymax, 50))

panel_deciles <- function(dec, var, umbral, titulo, xlab) {
  rho <- data %>%
    group_by(rep) %>%
    summarise(rho = cor(.data[[var]], var_rank_sw, method = "spearman"), .groups = "drop") %>%
    mutate(label = sprintf("rho == %.2f", rho))
  ggplot(dec, aes(x, m)) +
    geom_hline(yintercept = 50, linetype = "dotted", col = "grey60") +
    geom_vline(xintercept = umbral, linetype = "dashed", col = "grey40") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, col = "#216869") +
    geom_line(col = "#216869") +
    geom_point(col = "#14AFB2", size = 3) +
    geom_text(data = rho, aes(x = Inf, y = Inf, label = label), parse = TRUE,
      hjust = 1.1, vjust = 1.5, size = 9, inherit.aes = FALSE) +
    scale_y_continuous(limits = y_lims) +
    theme_pubclean() +
    theme(text = element_text(size = 30), panel.spacing = unit(3, "lines")) +
    facet_wrap(~rep) +
    labs(x = xlab, y = "Rango de ruido\nmedio") +
    ggtitle(titulo)
}

panel_cpg <- panel_deciles(dec_cpg, "cpg_oe", 0.6, "Razón CpG observado/esperado", "Razón CpG o/e del fragmento (mediana del decil)")
panel_gc <- panel_deciles(dec_gc, "g_c", 0.5, "Contenido G+C", "Contenido G+C del fragmento (mediana del decil)")

for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("ruido_cpg_oe_deciles.", ext)), panel_cpg, width = 12, height = 6.75, units = "in")
  ggsave(file.path(out_dir, paste0("ruido_gc_deciles.", ext)), panel_gc, width = 12, height = 6.75, units = "in")
  # fila vacia intermedia para separar los dos paneles (mismo formato que cgi_composicion_deciles)
  ggsave(file.path(out_dir, paste0("ruido_composicion_deciles.", ext)),
    panel_cpg / plot_spacer() / panel_gc + plot_layout(heights = c(1, 0.08, 1)),
    width = 12, height = 14.2, units = "in")
}

# Correlaciones de Spearman del ruido con la razon CpG o/e y con el
# contenido de G+C, por replica: en todos los promotores y dentro de los que
# cumplen / no cumplen CGI_frag (var_rank_sw del set completo, sin
# re-binning dentro del grupo).
cor_spearman <- function(df, var) {
  ct <- suppressWarnings(cor.test(df[[var]], df$var_rank_sw, method = "spearman", exact = FALSE))
  tibble(rho = unname(ct$estimate), pval = ct$p.value)
}
ruido_cor <- data %>%
  mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")) %>%
  bind_rows(mutate(., grupo = "Todos")) %>%
  group_by(rep, grupo) %>%
  reframe(
    n = n(),
    bind_rows(
      mutate(cor_spearman(pick(everything()), "cpg_oe"), variable = "cpg_oe"),
      mutate(cor_spearman(pick(everything()), "g_c"), variable = "g_c")
    )
  ) %>%
  mutate(grupo = factor(grupo, levels = c("Todos", "Con isla CpG", "Sin isla CpG"))) %>%
  arrange(rep, variable, grupo) %>%
  select(rep, variable, grupo, n, rho, pval)
write_tsv(ruido_cor, file.path(out_dir, "ruido_composicion_correlaciones.tsv"))
print(ruido_cor, n = Inf)

message("Figuras guardadas en ", out_dir)
