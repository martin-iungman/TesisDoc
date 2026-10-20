# EXPLORATORIO - relacion del ruido desacoplado de la actividad
# (var_rank_sw) con la composicion del fragmento de 252 pb: razon CpG o/e
# (cpg_oe) y contenido G+C (g_c), y si esa relacion difiere de la que tienen
# esas mismas variables con la actividad (motivos_cgi_actividad: aumento
# gradual de la actividad con cpg_oe, que sigue dentro de las islas, y g_c
# sin efecto dentro de las islas).
#
# Datos: activity_stats_highconf + prom_df (type == "promoter"), cada
# replica por separado. var_rank_sw con add_noise_rank() sobre el set
# completo de cada replica (igual que ruido_cgi_tata); add_noise_rank()
# descarta el resto de n %% 100 promotores de menor actividad, por eso n =
# 9800 / 11600. Isla CpG = CGI_frag. Percentil de actividad = 100 *
# percent_rank(mean) dentro de cada replica (igual que motivos_cgi_actividad).
# La referencia nula del rango de ruido es 50.5 (media de 1..100).
#
# Puntos (prefijo de los archivos de salida):
#   p1  deciles de cpg_oe y g_c: rango de ruido medio (IC95 t) + Spearman.
#       Tambien con el percentil de actividad como respuesta, para comparar.
#   p2  escalon o gradiente: bins finos de cpg_oe (0.05) y g_c (0.025);
#       pendientes a cada lado del umbral de isla y barrido del umbral de
#       un modelo escalon (cual ajusta mejor).
#   p3  estratificado por CGI_frag: Spearman y deciles dentro de cada
#       grupo, con var_rank_sw original y re-calculado dentro del subgrupo.
#   p4  modelos lineales (a)-(f), coeficientes estandarizados, R2 y F de
#       modelos anidados, para ruido y para actividad.
#   p5  por tercil de actividad: Spearman, AUC de CGI_frag (y TATA_EPD para
#       comparar) con roc_curve()/noise_auc() + IC DeLong (pROC), curvas ROC
#       en el formato de ROC_TATA_expr, y AUC por decil de actividad.
#   p6  reproducibilidad entre replicas por tercil de actividad, AUC sobre
#       el ruido consenso de ambas replicas y AUC corregida por atenuacion.
#
# Salidas en figures/exploratorio_ruido_composicion/ (sin fig_dir: no es
# figura de la tesis). Requiere: data/processed/activity_stats_highconf.tsv
# y data/processed/prom_df.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
library(pROC)
source("R/functions/plot_helpers.R")

out_dir <- "figures/exploratorio_ruido_composicion"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

REF_NULA <- 50.5
UMBRAL_CPG <- 0.6
UMBRAL_GC <- 0.5
clr_rep <- c("Rep 1" = "#0D2C54", "Rep 2" = "#AD343E")
var_lab <- c(cpg_oe = "CpG o/e", g_c = "Contenido G+C")
resp_lab <- c(ruido = "Ruido (rango de varianza)", actividad = "Actividad (percentil)")
tercil_lab <- c("1" = "Baja", "2" = "Media", "3" = "Alta")

# --- Datos ---------------------------------------------------------------------

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
base <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  # algunas columnas logicas de prom_df se leen con tipo mixto
  mutate(across(c(CGI_frag, TATA_EPD), ~ as.logical(as.character(.x))))

preparar <- function(df) {
  df %>%
    add_noise_rank() %>%
    group_by(rep) %>%
    mutate(act_pct = 100 * percent_rank(mean), tercil = ntile(act_pct, 3), decil_act = ntile(act_pct, 10)) %>%
    ungroup()
}

data <- preparar(base)
n_na <- sum(is.na(data$cpg_oe) | is.na(data$g_c) | is.na(data$CGI_frag))
message("Promotores por replica: ", paste(table(data$rep), collapse = " / "), " (", n_na, " con NA en composicion)")
data <- data %>% filter(!is.na(cpg_oe), !is.na(g_c), !is.na(CGI_frag))

# formato largo con las dos respuestas (ruido y actividad), misma escala 0-100
largo_resp <- function(df) {
  df %>%
    transmute(seq_id, rep, CGI_frag, cpg_oe, g_c, tercil, ruido = var_rank_sw, actividad = act_pct) %>%
    pivot_longer(c(ruido, actividad), names_to = "resp", values_to = "y")
}

# --- Helpers estadisticos -------------------------------------------------------

# Spearman con IC95 (transformacion z de Fisher con el SE de Bonett & Wright)
rho_ci <- function(x, y) {
  ok <- complete.cases(x, y)
  n <- sum(ok)
  r <- cor(x[ok], y[ok], method = "spearman")
  se <- sqrt((1 + r^2 / 2) / (n - 3))
  p <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)$p.value)
  tibble(n = n, rho = r, ic95_inf = tanh(atanh(r) - qnorm(0.975) * se), ic95_sup = tanh(atanh(r) + qnorm(0.975) * se), pval = p)
}

media_ci <- function(y) {
  n <- length(y)
  m <- mean(y)
  h <- qt(0.975, n - 1) * sd(y) / sqrt(n)
  tibble(n = n, media = m, ic95_inf = m - h, ic95_sup = m + h)
}

# AUC con el mismo metodo que ruido_cgi_tata (roc_curve/noise_auc) + IC y
# p-valor de DeLong (pROC, como ruido_summary). direction = "<": AUC =
# P(ruido | feature TRUE > ruido | feature FALSE), misma convencion que
# roc_curve() (AUC > 0.5 = la feature se asocia a ruido alto).
auc_ci <- function(df, feature) {
  resp <- df[[feature]]
  ok <- !is.na(resp)
  d <- df[ok, ]
  roc_obj <- pROC::roc(resp[ok], d$var_rank_sw, levels = c(FALSE, TRUE), direction = "<", ci = TRUE, quiet = TRUE)
  se <- sqrt(pROC::var(roc_obj, method = "delong"))
  auc_p <- as.numeric(roc_obj$auc)
  tibble(
    n = nrow(d), n_true = sum(resp[ok]), prevalencia = mean(resp[ok]),
    AUC = noise_auc(roc_curve(d, feature)), AUC_pROC = auc_p,
    ic95_inf = as.numeric(roc_obj$ci[1]), ic95_sup = as.numeric(roc_obj$ci[3]),
    pval = 2 * pnorm(-abs((auc_p - 0.5) / se))
  )
}

# =============================================================================
# p1 - Deciles de cpg_oe y g_c
# =============================================================================

deciles_resumen <- function(df_largo, grupos = NULL) {
  df_largo %>%
    pivot_longer(c(cpg_oe, g_c), names_to = "variable", values_to = "x") %>%
    group_by(across(all_of(c(grupos, "resp", "rep", "variable")))) %>%
    mutate(decil = ntile(x, 10)) %>%
    group_by(across(all_of(c(grupos, "resp", "rep", "variable", "decil")))) %>%
    summarise(x_mediana = median(x), media_ci(y), .groups = "drop")
}

rho_resumen <- function(df_largo, grupos = NULL) {
  df_largo %>%
    pivot_longer(c(cpg_oe, g_c), names_to = "variable", values_to = "x") %>%
    group_by(across(all_of(c(grupos, "resp", "rep", "variable")))) %>%
    summarise(rho_ci(x, y), .groups = "drop")
}

d_largo <- largo_resp(data)
p1_deciles <- deciles_resumen(d_largo)
p1_rho <- rho_resumen(d_largo)
write_tsv(p1_deciles, file.path(out_dir, "p1_deciles.tsv"))
write_tsv(p1_rho, file.path(out_dir, "p1_spearman.tsv"))

plot_deciles <- function(dec, rho, titulo, facet_formula = ~variable, ylab = "Rango de ruido medio (IC95)") {
  rho_lab <- rho %>%
    group_by(across(-c(n, rho, ic95_inf, ic95_sup, pval, rep))) %>%
    summarise(lab = paste0("ρ = ", paste(sprintf("%.3f", rho), collapse = " / ")), .groups = "drop")
  ggplot(dec, aes(x_mediana, media, col = rep)) +
    geom_hline(yintercept = REF_NULA, linetype = "dashed", col = "grey40") +
    geom_errorbar(aes(ymin = ic95_inf, ymax = ic95_sup), width = 0, linewidth = 0.6) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_text(data = rho_lab, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -0.05, vjust = 1.4, size = 4.2) +
    facet_grid(facet_formula, scales = "free", labeller = labeller(variable = var_lab, resp = resp_lab)) +
    scale_color_manual(values = clr_rep) +
    labs(x = "Mediana del decil", y = ylab, col = "", title = titulo, subtitle = "ρ de Spearman: Rep 1 / Rep 2") +
    theme_pubr(base_size = 16) +
    theme(strip.text = element_text(size = 14))
}

ggsave(
  file.path(out_dir, "p1_deciles_cpg_gc.jpg"),
  plot_deciles(filter(p1_deciles, resp == "ruido"), filter(p1_rho, resp == "ruido"), "Ruido por decil de composición", . ~ variable),
  width = 11, height = 5.5, units = "in"
)
ggsave(
  file.path(out_dir, "p1_deciles_ruido_vs_actividad.jpg"),
  plot_deciles(p1_deciles, p1_rho, "Composición: ruido vs. actividad", resp ~ variable, ylab = "Media del decil (IC95)"),
  width = 11, height = 9, units = "in"
)

# =============================================================================
# p2 - Escalon o gradiente
# =============================================================================

bins_def <- list(
  cpg_oe = list(breaks = seq(0.2, 1.0, 0.05), umbral = UMBRAL_CPG),
  g_c = list(breaks = seq(0.40, 0.80, 0.025), umbral = UMBRAL_GC)
)

p2_bins <- imap(bins_def, function(def, v) {
  d_largo %>%
    mutate(x = .data[[v]], bin = cut(x, def$breaks, include.lowest = TRUE)) %>%
    filter(!is.na(bin)) %>%
    group_by(resp, rep, bin) %>%
    summarise(x_centro = mean(range(def$breaks[as.integer(bin[1]) + 0:1])), media_ci(y), .groups = "drop") %>%
    mutate(variable = v)
}) %>% list_rbind()
write_tsv(p2_bins, file.path(out_dir, "p2_bins_finos.tsv"))

plot_bins <- function(v, resp_sel = "ruido") {
  d <- p2_bins %>% filter(variable == v, resp == resp_sel)
  u <- bins_def[[v]]$umbral
  ancho <- diff(bins_def[[v]]$breaks)[1]
  top <- ggplot(d, aes(x_centro, media, col = rep)) +
    geom_hline(yintercept = REF_NULA, linetype = "dashed", col = "grey40") +
    geom_vline(xintercept = u, linetype = "dotted", linewidth = 0.8) +
    geom_errorbar(aes(ymin = ic95_inf, ymax = ic95_sup), width = 0, position = position_dodge(ancho * 0.4)) +
    geom_line(position = position_dodge(ancho * 0.4)) +
    geom_point(size = 2, position = position_dodge(ancho * 0.4)) +
    scale_color_manual(values = clr_rep) +
    labs(x = NULL, y = paste0(resp_lab[resp_sel], "\nmedio (IC95)"), col = "", title = var_lab[v]) +
    theme_pubr(base_size = 15) +
    theme(axis.text.x = element_blank())
  bottom <- ggplot(d, aes(x_centro, n, fill = rep)) +
    geom_vline(xintercept = u, linetype = "dotted", linewidth = 0.8) +
    geom_col(position = position_dodge(ancho * 0.9), width = ancho * 0.9) +
    scale_fill_manual(values = clr_rep, guide = "none") +
    labs(x = var_lab[v], y = "n") +
    theme_pubr(base_size = 15)
  top / bottom + plot_layout(heights = c(3, 1))
}

ggsave(
  file.path(out_dir, "p2_bins_finos_ruido.jpg"),
  (plot_bins("cpg_oe") | plot_bins("g_c")) + plot_annotation(title = "Ruido en bins finos de composición (línea punteada = umbral de isla CpG)"),
  width = 14, height = 7, units = "in"
)
ggsave(
  file.path(out_dir, "p2_bins_finos_actividad.jpg"),
  (plot_bins("cpg_oe", "actividad") | plot_bins("g_c", "actividad")) + plot_annotation(title = "Actividad en bins finos de composición (comparación)"),
  width = 14, height = 7, units = "in"
)

# Pendientes a cada lado del umbral (por 0.1 de la variable) y salto en el
# umbral, en el rango de los bins: y ~ lado + lado:x (modelo segmentado
# discontinuo). Si el cambio es un escalon, las pendientes son ~0 y el salto
# grande; si es gradual, pendientes parecidas y salto ~0.
p2_segmentos <- imap(bins_def, function(def, v) {
  d_largo %>%
    mutate(x = .data[[v]]) %>%
    filter(x >= min(def$breaks), x <= max(def$breaks)) %>%
    mutate(lado = factor(ifelse(x > def$umbral, "arriba", "abajo"), levels = c("abajo", "arriba")), xc = (x - def$umbral) / 0.1) %>%
    group_by(resp, rep) %>%
    group_modify(function(dd, k) {
      m <- lm(y ~ lado + lado:xc, data = dd)
      ci <- confint(m)
      tibble(
        termino = c("salto_en_umbral", "pendiente_abajo", "pendiente_arriba"),
        estimado = coef(m)[c("ladoarriba", "ladoabajo:xc", "ladoarriba:xc")],
        ic95_inf = ci[c("ladoarriba", "ladoabajo:xc", "ladoarriba:xc"), 1],
        ic95_sup = ci[c("ladoarriba", "ladoabajo:xc", "ladoarriba:xc"), 2],
        n_abajo = sum(dd$lado == "abajo"), n_arriba = sum(dd$lado == "arriba")
      )
    }) %>%
    ungroup() %>%
    mutate(variable = v, umbral = def$umbral)
}) %>% list_rbind()
write_tsv(p2_segmentos, file.path(out_dir, "p2_segmentos_umbral.tsv"))

# Barrido del umbral de un modelo escalon y ~ I(x > t): R2 por t. El t de R2
# maximo dice donde esta el cambio, independientemente de la definicion de
# isla; se compara con el R2 de un modelo lineal y ~ x en el mismo rango.
p2_barrido <- imap(bins_def, function(def, v) {
  grilla <- seq(quantile(data[[v]], 0.05), quantile(data[[v]], 0.95), by = 0.01)
  d_largo %>%
    mutate(x = .data[[v]]) %>%
    group_by(resp, rep) %>%
    group_modify(function(dd, k) {
      r2_lin <- summary(lm(y ~ x, data = dd))$r.squared
      map(grilla, ~ tibble(umbral = .x, r2_escalon = summary(lm(y ~ I(x > .x), data = dd))$r.squared)) %>%
        list_rbind() %>%
        mutate(r2_lineal = r2_lin)
    }) %>%
    ungroup() %>%
    mutate(variable = v)
}) %>% list_rbind()
p2_barrido_max <- p2_barrido %>%
  group_by(variable, resp, rep) %>%
  slice_max(r2_escalon, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(umbral_optimo = umbral)
write_tsv(p2_barrido_max, file.path(out_dir, "p2_umbral_optimo.tsv"))

p2_barrido_plot <- p2_barrido %>%
  ggplot(aes(umbral, r2_escalon, col = rep)) +
  geom_line(linewidth = 1) +
  geom_hline(aes(yintercept = r2_lineal, col = rep), linetype = "dashed") +
  geom_vline(data = tibble(variable = names(bins_def), u = map_dbl(bins_def, "umbral")), aes(xintercept = u), linetype = "dotted") +
  facet_grid(resp ~ variable, scales = "free", labeller = labeller(variable = var_lab, resp = resp_lab)) +
  scale_color_manual(values = clr_rep) +
  labs(
    x = "Umbral del modelo escalón", y = "R²", col = "",
    title = "¿Dónde está el cambio?", subtitle = "Línea sólida: R² de y ~ I(x > umbral); discontinua: R² de y ~ x; punteada: umbral de isla"
  ) +
  theme_pubr(base_size = 15)
ggsave(file.path(out_dir, "p2_barrido_umbral.jpg"), p2_barrido_plot, width = 11, height = 8, units = "in")

# =============================================================================
# p3 - Estratificado por isla CpG
# =============================================================================

grupo_lab <- function(cgi) ifelse(cgi, "Con isla CpG", "Sin isla CpG")

# var_rank_sw original (bins del set completo) y re-calculado dentro de cada
# subgrupo (re-binning, como el panel no-CGI de chipatlas_histonas)
d_rebin <- map(c(TRUE, FALSE), ~ base %>% filter(CGI_frag == .x) %>% preparar()) %>%
  list_rbind() %>%
  filter(!is.na(cpg_oe), !is.na(g_c))

d_estrat <- bind_rows(
  largo_resp(data) %>% mutate(ranking = "original"),
  largo_resp(d_rebin) %>% filter(resp == "ruido") %>% mutate(ranking = "re-calculado en el subgrupo")
) %>%
  mutate(grupo = grupo_lab(CGI_frag))

p3_rho <- rho_resumen(d_estrat, grupos = c("ranking", "grupo"))
p3_deciles <- deciles_resumen(d_estrat, grupos = c("ranking", "grupo"))
write_tsv(p3_rho, file.path(out_dir, "p3_spearman_por_isla.tsv"))
write_tsv(p3_deciles, file.path(out_dir, "p3_deciles_por_isla.tsv"))

ggsave(
  file.path(out_dir, "p3_deciles_por_isla_ruido.jpg"),
  plot_deciles(
    p3_deciles %>% filter(resp == "ruido", ranking == "original"),
    p3_rho %>% filter(resp == "ruido", ranking == "original"),
    "Ruido por decil de composición, dentro de cada grupo", grupo ~ variable
  ),
  width = 11, height = 9, units = "in"
)
ggsave(
  file.path(out_dir, "p3_deciles_por_isla_ruido_rebin.jpg"),
  plot_deciles(
    p3_deciles %>% filter(resp == "ruido", ranking != "original"),
    p3_rho %>% filter(resp == "ruido", ranking != "original"),
    "Ruido re-calculado dentro de cada grupo (re-binning)", grupo ~ variable
  ),
  width = 11, height = 9, units = "in"
)
ggsave(
  file.path(out_dir, "p3_deciles_por_isla_actividad.jpg"),
  plot_deciles(
    p3_deciles %>% filter(resp == "actividad"), p3_rho %>% filter(resp == "actividad"),
    "Actividad por decil de composición, dentro de cada grupo (comparación)", grupo ~ variable,
    ylab = "Percentil de actividad medio (IC95)"
  ),
  width = 11, height = 9, units = "in"
)

# =============================================================================
# p4 - Modelos lineales
# =============================================================================
# Respuesta y predictores continuos estandarizados (z) dentro de los datos de
# cada modelo; CGI_frag como 0/1 (su coeficiente queda en DE de la respuesta).

modelos_def <- list(
  a = list(f = y ~ cpg_oe, sub = "todos"),
  b = list(f = y ~ g_c, sub = "todos"),
  c = list(f = y ~ cpg_oe + g_c, sub = "todos"),
  d = list(f = y ~ cgi, sub = "todos"),
  e = list(f = y ~ cgi + cpg_oe + g_c, sub = "todos"),
  f = list(f = y ~ cpg_oe + g_c, sub = "isla"),
  f_cpg = list(f = y ~ cpg_oe, sub = "isla"),
  f_gc = list(f = y ~ g_c, sub = "isla")
)
comparaciones <- tribble(
  ~reducido, ~completo, ~pregunta,
  "a", "c", "g_c suma sobre cpg_oe",
  "b", "c", "cpg_oe suma sobre g_c",
  "d", "e", "continuas suman sobre la isla binaria",
  "c", "e", "isla binaria suma sobre las continuas",
  "f_cpg", "f", "dentro de islas: g_c suma sobre cpg_oe",
  "f_gc", "f", "dentro de islas: cpg_oe suma sobre g_c"
)

z <- function(x) (x - mean(x)) / sd(x)

ajustar <- function(dd) {
  map(modelos_def, function(def) {
    d <- if (def$sub == "isla") filter(dd, CGI_frag) else dd
    d <- d %>% transmute(y = z(y), cpg_oe = z(cpg_oe), g_c = z(g_c), cgi = as.numeric(CGI_frag))
    lm(def$f, data = d)
  })
}

p4 <- d_largo %>%
  group_by(resp, rep) %>%
  group_map(function(dd, k) {
    ms <- ajustar(dd)
    coefs <- imap(ms, function(m, id) {
      ci <- confint(m)
      tibble(modelo = id, termino = names(coef(m))[-1], beta = coef(m)[-1], ic95_inf = ci[-1, 1], ic95_sup = ci[-1, 2], pval = summary(m)$coefficients[-1, 4])
    }) %>% list_rbind()
    r2 <- imap(ms, ~ tibble(modelo = .y, n = nobs(.x), r2 = summary(.x)$r.squared, r2_aj = summary(.x)$adj.r.squared)) %>% list_rbind()
    f <- pmap(comparaciones, function(reducido, completo, pregunta) {
      a <- anova(ms[[reducido]], ms[[completo]])
      tibble(reducido, completo, pregunta, delta_r2 = summary(ms[[completo]])$r.squared - summary(ms[[reducido]])$r.squared, F = a$F[2], df = a$Df[2], pval = a$`Pr(>F)`[2])
    }) %>% list_rbind()
    list(coefs = bind_cols(k, coefs), r2 = bind_cols(k, r2), f = bind_cols(k, f))
  })
p4_coefs <- map(p4, "coefs") %>% list_rbind()
p4_r2 <- map(p4, "r2") %>% list_rbind()
p4_f <- map(p4, "f") %>% list_rbind()
write_tsv(p4_coefs, file.path(out_dir, "p4_coeficientes.tsv"))
write_tsv(p4_r2, file.path(out_dir, "p4_r2.tsv"))
write_tsv(p4_f, file.path(out_dir, "p4_F_anidados.tsv"))

p4_cor_pred <- data %>%
  group_by(rep, grupo = grupo_lab(CGI_frag)) %>%
  summarise(r_cpg_gc = cor(cpg_oe, g_c), .groups = "drop") %>%
  bind_rows(data %>% group_by(rep) %>% summarise(r_cpg_gc = cor(cpg_oe, g_c)) %>% mutate(grupo = "Todos"))
write_tsv(p4_cor_pred, file.path(out_dir, "p4_correlacion_predictores.tsv"))

term_lab <- c(cpg_oe = "CpG o/e", g_c = "G+C", cgi = "Isla CpG")
mod_lab <- c(
  a = "(a) CpG o/e", b = "(b) G+C", c = "(c) CpG o/e + G+C", d = "(d) Isla",
  e = "(e) Isla + CpG o/e + G+C", f = "(f) islas: CpG o/e + G+C"
)
p4_plot <- p4_coefs %>%
  filter(modelo %in% names(mod_lab)) %>%
  left_join(p4_r2, by = c("resp", "rep", "modelo")) %>%
  mutate(modelo = factor(modelo, levels = rev(names(mod_lab)))) %>%
  ggplot(aes(beta, modelo, col = term_lab[termino], shape = rep)) +
  geom_vline(xintercept = 0, col = "grey50") +
  geom_pointrange(aes(xmin = ic95_inf, xmax = ic95_sup), position = position_dodge(0.6)) +
  facet_wrap(~resp, labeller = labeller(resp = resp_lab), scales = "free_x") +
  scale_y_discrete(labels = mod_lab) +
  scale_color_manual(values = c("CpG o/e" = "#216869", "G+C" = "#D6741F", "Isla CpG" = thesis_clr)) +
  labs(x = "Coeficiente estandarizado (IC95)", y = NULL, col = "Término", shape = "", title = "Modelos lineales de composición") +
  theme_pubr(base_size = 15)
ggsave(file.path(out_dir, "p4_coeficientes.jpg"), p4_plot, width = 13, height = 7, units = "in")

# =============================================================================
# p5 - Dependencia de la actividad (terciles de percentil de actividad)
# =============================================================================

p5_rho <- d_largo %>%
  filter(resp == "ruido") %>%
  mutate(tercil = tercil_lab[as.character(tercil)]) %>%
  rho_resumen(grupos = "tercil") %>%
  bind_rows(
    d_largo %>% filter(resp == "ruido") %>% mutate(tercil = tercil_lab[as.character(tercil)], grupo = grupo_lab(CGI_frag)) %>%
      rho_resumen(grupos = c("tercil", "grupo"))
  ) %>%
  mutate(grupo = replace_na(grupo, "Todos"))
write_tsv(p5_rho, file.path(out_dir, "p5_spearman_por_tercil.tsv"))

p5_auc <- data %>%
  group_by(rep, tercil) %>%
  group_modify(~ bind_rows(auc_ci(.x, "CGI_frag") %>% mutate(feature = "CGI_frag"), auc_ci(.x, "TATA_EPD") %>% mutate(feature = "TATA_EPD"))) %>%
  ungroup() %>%
  mutate(tercil = factor(tercil_lab[as.character(tercil)], levels = tercil_lab)) %>%
  bind_rows(
    data %>% group_by(rep) %>%
      group_modify(~ bind_rows(auc_ci(.x, "CGI_frag") %>% mutate(feature = "CGI_frag"), auc_ci(.x, "TATA_EPD") %>% mutate(feature = "TATA_EPD"))) %>%
      ungroup() %>% mutate(tercil = "Todos")
  )
write_tsv(p5_auc, file.path(out_dir, "p5_auc_por_tercil.tsv"))

p5_auc_decil <- data %>%
  group_by(rep, decil_act) %>%
  group_modify(~ bind_rows(auc_ci(.x, "CGI_frag") %>% mutate(feature = "CGI_frag"), auc_ci(.x, "TATA_EPD") %>% mutate(feature = "TATA_EPD"))) %>%
  ungroup()
write_tsv(p5_auc_decil, file.path(out_dir, "p5_auc_por_decil_actividad.tsv"))

# Curvas ROC por tercil de actividad, en el formato de noise_roc_by_expr()
# (ROC_TATA_expr): color = replica, linetype = nivel de actividad.
noise_roc_by_tercil <- function(data, feature, titulo) {
  curvas <- data %>%
    mutate(expr_level = factor(tercil_lab[as.character(tercil)], levels = tercil_lab)) %>%
    group_split(rep, expr_level) %>%
    purrr::map(~ roc_curve(.x, feature) %>% mutate(rep = unique(.x$rep), expr_level = unique(.x$expr_level), AUC = noise_auc(.))) %>%
    purrr::list_rbind()
  etiquetas <- curvas %>%
    distinct(AUC, rep, expr_level) %>%
    arrange(rep, expr_level) %>%
    mutate(AUC = paste0(rep, ", ", tolower(expr_level), ": ", round(AUC, 3)), TPR = seq(0.95, by = -0.07, length.out = n()), FPR = 0.02)
  ggplot(curvas, aes(FPR, TPR, col = rep, linetype = expr_level)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", col = "grey50") +
    geom_text(data = etiquetas, aes(label = AUC), hjust = "left", size = 3.4, show.legend = FALSE) +
    ggtitle(titulo) +
    scale_color_manual(values = clr_rep) +
    scale_linetype_manual(values = c(Baja = "dotted", Media = "dashed", Alta = "solid")) +
    theme_pubr(base_size = 16) +
    theme(legend.box = "vertical") +
    labs(col = "", linetype = "Actividad")
}

ggsave(file.path(out_dir, "p5_ROC_CGI_tercil.jpg"), noise_roc_by_tercil(data, "CGI_frag", "Curva ROC - Isla CpG"), width = 6, height = 6, units = "in")
ggsave(file.path(out_dir, "p5_ROC_TATA_tercil.jpg"), noise_roc_by_tercil(data, "TATA_EPD", "Curva ROC - TATA-box"), width = 6, height = 6, units = "in")

feat_lab <- c(CGI_frag = "Isla CpG", TATA_EPD = "TATA-box")
p5_auc_plot <- p5_auc %>%
  mutate(tercil = factor(tercil, levels = c(tercil_lab, "Todos"))) %>%
  ggplot(aes(tercil, AUC, col = rep)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "grey40") +
  geom_pointrange(aes(ymin = ic95_inf, ymax = ic95_sup), position = position_dodge(0.4)) +
  facet_wrap(~feature, labeller = labeller(feature = feat_lab)) +
  scale_color_manual(values = clr_rep) +
  labs(x = "Tercil de actividad", y = "AUC (ruido alto) - IC95 DeLong", col = "", title = "Efecto sobre el ruido por nivel de actividad") +
  theme_pubr(base_size = 16)
ggsave(file.path(out_dir, "p5_auc_por_tercil.jpg"), p5_auc_plot, width = 9, height = 5.5, units = "in")

p5_auc_decil_plot <- p5_auc_decil %>%
  ggplot(aes(decil_act, AUC, col = rep)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "grey40") +
  geom_pointrange(aes(ymin = ic95_inf, ymax = ic95_sup), position = position_dodge(0.4)) +
  geom_line(position = position_dodge(0.4)) +
  facet_wrap(~feature, labeller = labeller(feature = feat_lab)) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_manual(values = clr_rep) +
  labs(x = "Decil de actividad", y = "AUC (ruido alto) - IC95 DeLong", col = "", title = "Efecto sobre el ruido por decil de actividad") +
  theme_pubr(base_size = 16)
ggsave(file.path(out_dir, "p5_auc_por_decil_actividad.jpg"), p5_auc_decil_plot, width = 11, height = 5.5, units = "in")

# =============================================================================
# p6 - Reproducibilidad del rango de ruido por nivel de actividad
# =============================================================================
# Promotores presentes en ambas replicas; tercil por actividad promedio
# (promedio del percentil en ambas replicas, como metrica_ruido).

wide <- data %>%
  select(seq_id, rep, CGI_frag, TATA_EPD, mean, var_rank_sw) %>%
  mutate(rep = recode(rep, "Rep 1" = "r1", "Rep 2" = "r2")) %>%
  pivot_wider(names_from = rep, values_from = c(mean, var_rank_sw)) %>%
  drop_na() %>%
  mutate(
    actividad = (percent_rank(mean_r1) + percent_rank(mean_r2)) / 2,
    tercil = ntile(actividad, 3),
    # ruido consenso: rango percentil del promedio de ambas replicas (mas
    # confiable que cualquiera de las dos solas)
    var_rank_cons = 100 * percent_rank((var_rank_sw_r1 + var_rank_sw_r2) / 2)
  )

p6_repro <- bind_rows(
  wide %>% group_by(tercil) %>% summarise(rho_ci(var_rank_sw_r1, var_rank_sw_r2), .groups = "drop") %>% mutate(grupo = "Todos"),
  wide %>% group_by(tercil, grupo = grupo_lab(CGI_frag)) %>% summarise(rho_ci(var_rank_sw_r1, var_rank_sw_r2), .groups = "drop")
) %>%
  mutate(tercil = tercil_lab[as.character(tercil)])
write_tsv(p6_repro, file.path(out_dir, "p6_reproducibilidad_por_tercil.tsv"))

# AUC de CGI_frag por tercil (actividad promedio) con el ruido de cada
# replica y con el consenso, y AUC corregida por atenuacion: d = sqrt(2) *
# qnorm(AUC) (d de Cohen bajo normalidad), d_corr = d / sqrt(confiabilidad),
# con confiabilidad de una replica = rho entre replicas del tercil (para el
# consenso, Spearman-Brown: 2 rho / (1 + rho)). AUC_corr = pnorm(d_corr / sqrt(2)).
p6_auc <- wide %>%
  group_by(tercil) %>%
  group_modify(function(dd, k) {
    rel <- cor(dd$var_rank_sw_r1, dd$var_rank_sw_r2, method = "spearman")
    map(c(r1 = "var_rank_sw_r1", r2 = "var_rank_sw_r2", consenso = "var_rank_cons"), function(col) {
      auc_ci(dd %>% mutate(var_rank_sw = .data[[col]]), "CGI_frag")
    }) %>%
      list_rbind(names_to = "ruido") %>%
      mutate(
        confiabilidad = ifelse(ruido == "consenso", 2 * rel / (1 + rel), rel),
        d = sqrt(2) * qnorm(AUC_pROC),
        AUC_corr_atenuacion = pnorm(d / sqrt(pmax(confiabilidad, 0.01)) / sqrt(2))
      )
  }) %>%
  ungroup() %>%
  mutate(tercil = tercil_lab[as.character(tercil)])
write_tsv(p6_auc, file.path(out_dir, "p6_auc_cgi_consenso_atenuacion.tsv"))

# --- Resumen en consola ----------------------------------------------------------

cat("\n== p1 Spearman ==\n"); print(p1_rho %>% select(resp, rep, variable, n, rho, ic95_inf, ic95_sup))
cat("\n== p1 deciles (ruido) ==\n"); print(p1_deciles %>% filter(resp == "ruido") %>% select(rep, variable, decil, media) %>% pivot_wider(names_from = decil, values_from = media) %>% mutate(across(where(is.numeric), ~ round(.x, 1))))
cat("\n== p2 segmentos ==\n"); print(p2_segmentos %>% select(variable, resp, rep, termino, estimado, ic95_inf, ic95_sup), n = Inf)
cat("\n== p2 umbral optimo ==\n"); print(p2_barrido_max)
cat("\n== p3 Spearman por isla ==\n"); print(p3_rho %>% select(ranking, grupo, resp, rep, variable, n, rho, ic95_inf, ic95_sup), n = Inf)
cat("\n== p4 R2 ==\n"); print(p4_r2 %>% mutate(across(c(r2, r2_aj), ~ round(.x, 4))) %>% select(-r2_aj) %>% pivot_wider(names_from = modelo, values_from = c(r2)), n = Inf)
cat("\n== p4 coeficientes ==\n"); print(p4_coefs %>% filter(modelo %in% c("c", "e", "f")) %>% mutate(across(where(is.double), ~ signif(.x, 3))), n = Inf)
cat("\n== p4 F ==\n"); print(p4_f %>% mutate(across(where(is.double), ~ signif(.x, 3))), n = Inf)
cat("\n== p4 cor(cpg_oe, g_c) ==\n"); print(p4_cor_pred)
cat("\n== p5 Spearman por tercil ==\n"); print(p5_rho %>% filter(grupo == "Todos") %>% select(tercil, rep, variable, rho, ic95_inf, ic95_sup), n = Inf)
cat("\n== p5 AUC por tercil ==\n"); print(p5_auc %>% select(feature, rep, tercil, n, prevalencia, AUC, AUC_pROC, ic95_inf, ic95_sup, pval), n = Inf)
cat("\n== p6 reproducibilidad ==\n"); print(p6_repro %>% select(tercil, grupo, n, rho, ic95_inf, ic95_sup), n = Inf)
cat("\n== p6 AUC consenso / atenuacion ==\n"); print(p6_auc %>% select(tercil, ruido, n, AUC_pROC, ic95_inf, ic95_sup, confiabilidad, AUC_corr_atenuacion), n = Inf)

message("Salidas guardadas en ", out_dir)
