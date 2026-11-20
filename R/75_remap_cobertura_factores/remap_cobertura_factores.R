# remap_cobertura_factores (ver docs/mapping_figuras.csv para el numero de
# figura vigente). Figura de Metodos ("Analisis de datos masivos de
# ChIP-seq"): por que los efectos por factor (unido vs. no unido) estan
# confundidos con la cobertura de cada factor. No usa ningun dato del ensayo
# (ni actividad ni ruido): activity_stats_highconf.tsv solo define que
# promotores forman el universo (alta confianza por replica).
#
# Matriz de union: data/processed/remap_tf_hits.tsv (ReMap2022 no redundante,
# unido = >= 1 pb de solapamiento con el fragmento de 252 pb; ver
# R/00_prom_features/build_remap_tf_hits.R), la misma que usan
# R/17_remap_activity_gsea y R/56_ruido_contexto_endogeno. Criterio de
# inclusion identico al de esos scripts: > 100 promotores unidos en ambas
# replicas, sin marcas de histona (ReMap no tiene ninguna; el filtro no saca
# nada) -> 982 factores.
#
# Panel A: histograma del numero de promotores unidos por factor (log10) en
#   el universo de la replica 1; factores incluidos en color, el resto en
#   gris; linea punteada en 100. Eje superior: % del universo.
# Panel B: por factor incluido, proporcion de promotores con isla CpG
#   (CGI_frag) y de promotores no detectados en FANTOM5
#   (sample_specificity_class == "non_detected") entre los NO unidos, vs. %
#   de promotores unidos; linea punteada = proporcion en todo el universo;
#   loess. Replica 1.
#
# Conteos por factor y replica cacheados en
# data/processed/remap_cobertura_factores_v1.tsv (leer remap_tf_hits.tsv
# tarda varios minutos; borrar el cache para recalcular). Numeros para el
# texto en remap_cobertura_numeros.tsv dentro de la carpeta de la figura.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "remap_cobertura_factores"
out_dir <- fig_dir(slug)
base_size <- 16
umbral <- 100

save_panel <- function(p, name, width = 9, height = 6.75) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), p, width = width, height = height, units = "in")
  }
}
num <- function(x, digits = 2) formatC(x, format = "f", digits = digits, big.mark = "") %>% str_replace("^-", "−")

# --- Universo y conteos por factor ------------------------------------------

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
universo <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE) %>%
  distinct(seq_id, rep) %>%
  semi_join(prom_df, by = "seq_id") %>%
  left_join(prom_df %>%
    transmute(seq_id, cgi = CGI_frag, no_det = sample_specificity_class == "non_detected"), by = "seq_id")

totales <- universo %>%
  group_by(rep) %>%
  summarise(
    N = n(),
    N_cgi = sum(!is.na(cgi)), cgi_tot = sum(cgi, na.rm = TRUE),
    N_nodet = sum(!is.na(no_det)), nodet_tot = sum(no_det, na.rm = TRUE)
  )

cache <- "data/processed/remap_cobertura_factores_v1.tsv"
if (file.exists(cache)) {
  conteos <- read_tsv(cache, show_col_types = FALSE)
} else {
  hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE, col_select = c(seq_id, TF)) %>%
    distinct(seq_id, TF)
  conteos <- universo %>%
    inner_join(hits, by = "seq_id", relationship = "many-to-many") %>%
    group_by(rep, TF) %>%
    summarise(
      n_unidos = n(),
      n_unidos_cgi = sum(!is.na(cgi)), cgi_unidos = sum(cgi, na.rm = TRUE),
      n_unidos_nodet = sum(!is.na(no_det)), nodet_unidos = sum(no_det, na.rm = TRUE),
      .groups = "drop"
    )
  write_tsv(conteos, cache)
}

# Mismo criterio que R/17 y R/56: > 100 unidos en ambas replicas, sin
# marcas de histona.
incluidos <- conteos %>%
  filter(n_unidos > umbral, !str_detect(TF, "^H\\d")) %>%
  count(TF) %>%
  filter(n == 2) %>%
  pull(TF)

por_factor <- conteos %>%
  left_join(totales, by = "rep") %>%
  mutate(
    incluido = TF %in% incluidos,
    cobertura = n_unidos / N,
    n_no_unidos = N - n_unidos,
    prop_cgi_no_unidos = (cgi_tot - cgi_unidos) / (N_cgi - n_unidos_cgi),
    prop_nodet_no_unidos = (nodet_tot - nodet_unidos) / (N_nodet - n_unidos_nodet)
  )

# --- Panel A: distribucion de la cobertura ----------------------------------

rep_a <- "Rep 1"
N1 <- totales$N[totales$rep == rep_a]
dat_a <- por_factor %>% filter(rep == rep_a)
clr_excl <- "grey70"

panel_a <- ggplot(dat_a, aes(n_unidos, fill = incluido)) +
  geom_histogram(binwidth = 0.1, boundary = 2, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = umbral, linetype = "dashed") +
  # n arriba y % del universo debajo en cada marca (sec_axis falla con la
  # combinacion scales 1.3.0 / ggplot2 3.5.2 de este entorno)
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, N1),
    labels = function(x) {
      pct <- x / N1 * 100
      paste0(num(x, 0), "\n", map2_chr(pct, if_else(pct < 0.1, 2, if_else(pct < 10, 1, 0)), num), "%")
    }
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(`TRUE` = thesis_clr, `FALSE` = clr_excl),
    labels = c(`TRUE` = "Incluidos", `FALSE` = "Excluidos"), breaks = c("TRUE", "FALSE"), name = NULL) +
  labs(x = "Promotores unidos por el factor\n(n y % de los promotores de la library)", y = "Factores") +
  theme_pubclean(base_size = base_size) +
  theme(legend.position = "inside", legend.position.inside = c(0.2, 0.85), legend.direction = "vertical")

save_panel(panel_a, "panel_a_cobertura")

# --- Panel B: composicion de los no unidos -----------------------------------

ref_b <- totales %>%
  filter(rep == rep_a) %>%
  transmute(`Islas CpG` = cgi_tot / N_cgi, `No detectados en FANTOM5` = nodet_tot / N_nodet) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "ref")

dat_b <- dat_a %>%
  filter(incluido) %>%
  select(TF, cobertura, n_no_unidos, `Islas CpG` = prop_cgi_no_unidos, `No detectados en FANTOM5` = prop_nodet_no_unidos) %>%
  pivot_longer(c(`Islas CpG`, `No detectados en FANTOM5`), names_to = "variable", values_to = "prop")

plot_b <- function(d, ref) {
  ggplot(d, aes(cobertura * 100, prop * 100)) +
    geom_hline(data = ref, aes(yintercept = ref * 100), linetype = "dashed") +
    geom_point(colour = thesis_clr, alpha = 0.4, size = 1.5) +
    geom_smooth(method = "loess", formula = y ~ x, colour = "black", se = TRUE, linewidth = 0.8) +
    scale_x_continuous(limits = c(0, 100), labels = function(x) num(x, 0)) +
    scale_y_continuous(labels = function(x) num(x, 0)) +
    labs(x = "Promotores unidos por el factor (%)", y = "Promotores entre los no unidos (%)") +
    theme_pubclean(base_size = base_size) +
    theme(strip.background = element_blank(), strip.text = element_text(size = base_size))
}

panel_b <- plot_b(dat_b, ref_b) + facet_wrap(~variable, scales = "free_y")
panel_b_cgi <- plot_b(filter(dat_b, variable == "Islas CpG"), filter(ref_b, variable == "Islas CpG")) +
  labs(y = "Promotores con isla CpG entre los no unidos (%)")
panel_b_nodet <- plot_b(filter(dat_b, variable != "Islas CpG"), filter(ref_b, variable != "Islas CpG")) +
  labs(y = "No detectados en FANTOM5 entre los no unidos (%)")

save_panel(panel_b, "panel_b_composicion_no_unidos", 12, 6.75)
save_panel(panel_b_cgi, "panel_b_cgi_no_unidos")
save_panel(panel_b_nodet, "panel_b_no_detectados_no_unidos")

figura <- (panel_a | panel_b) + plot_layout(widths = c(1, 1.6)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = base_size + 4, face = "bold"))
save_panel(figura, slug, 20, 6.75)

# --- Numeros para el texto ---------------------------------------------------

numeros <- por_factor %>%
  filter(incluido) %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    p95 <- quantile(d$cobertura, 0.95)
    alta <- d %>% filter(cobertura >= p95)
    ct_cgi <- cor.test(d$cobertura, d$prop_cgi_no_unidos, method = "spearman", exact = FALSE)
    ct_nodet <- cor.test(d$cobertura, d$prop_nodet_no_unidos, method = "spearman", exact = FALSE)
    tibble(
      n_factores = nrow(d),
      N_universo = d$N[1],
      cobertura_mediana = median(d$cobertura),
      cobertura_p95 = p95,
      n_no_unidos_p95 = d$N[1] - ceiling(p95 * d$N[1]),
      n_no_unidos_mediana_top5 = median(alta$n_no_unidos),
      n_no_unidos_min = min(d$n_no_unidos),
      tf_max_cobertura = d$TF[which.max(d$cobertura)],
      rho_cgi = unname(ct_cgi$estimate), p_cgi = ct_cgi$p.value,
      rho_nodet = unname(ct_nodet$estimate), p_nodet = ct_nodet$p.value,
      ref_cgi = d$cgi_tot[1] / d$N_cgi[1],
      ref_nodet = d$nodet_tot[1] / d$N_nodet[1],
      cgi_no_unidos_top5 = median(alta$prop_cgi_no_unidos),
      nodet_no_unidos_top5 = median(alta$prop_nodet_no_unidos),
      n_factores_total = n_distinct(conteos$TF[conteos$rep == k$rep])
    )
  }) %>%
  ungroup()
write_tsv(numeros, file.path(out_dir, "remap_cobertura_numeros.tsv"))
print(as.data.frame(numeros))
