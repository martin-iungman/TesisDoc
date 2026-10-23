# EXPLORATORIO - efecto ajustado de las caracteristicas de contexto endogeno
# (grupo "endo" de feature_groups, las de summary_features_endogeno) sobre la
# actividad: modelo lineal multivariable, por replica, del percentil de
# actividad con todas las caracteristicas a la vez como efectos principales,
# ajustado por composicion del fragmento (razon CpG o/e y contenido de G+C,
# continuas) en lugar de estratificar por isla CpG (ver R/39_contexto_cgi).
#
# Respuesta: percentil de actividad dentro de cada replica (100 * (rango
# promedio - 0.5) / n, sobre todos los promotores high-conf; mismo que
# R/38_aditividad). Coeficientes en puntos de percentil.
#
# Caracteristicas (15 en el summary) -> 13 en el modelo:
#   - Especificidad tisular: alta y baja son los dos terciles extremos del
#     mismo indice (phi -0.50). Se deja "Baja especificidad tisular"
#     (referencia = terciles medio y alto): menos colineal con la razon CpG o/e
#     (r 0.36 vs. -0.48 la alta).
#   - Forma: angostos y anchos son los terciles extremos del mismo ancho
#     interquantil (phi -0.50). Se deja "Promotores anchos" (referencia =
#     angostos e intermedios): menos colineal con CpG o/e (0.24 vs. -0.33) y con
#     especificidad tisular.
#   - El resto se mantiene (las tres ventanas de conservacion son regiones
#     distintas, phi 0.47-0.49 entre vecinas; HEK293 y DNase-seq miden cosas
#     distintas, phi 0.36). La colinealidad se reporta con el VIF.
#
# Datos faltantes (no son al azar; borrarlos sesga hacia promotores con mas
# expresion endogena):
#   - Especificidad tisular y forma no estan definidas para los promotores no
#     detectados en FANTOM5 (indice de Gini NA exactamente en los "No
#     detectado"): ahi se codifican FALSE; "No detectado (FANTOM5)" ya los
#     distingue.
#   - Modelo principal: los NA restantes (turnover sin anotacion de Young et
#     al., forma NA en promotores detectados, orientacion sin senal PRO-cap
#     suficiente) se codifican FALSE y se agrega un indicador de dato faltante
#     por variable ("Sin dato: ...").
#   - Sensibilidad: el mismo modelo en casos completos (sin indicadores).
#
# Efecto crudo = lm(percentil ~ caracteristica) en la misma muestra y con la
# misma codificacion que el ajustado. BH dentro de cada replica y modelo
# (familia = las 13 caracteristicas; covariables e indicadores fuera).
# "Sobrevive" = p BH < 0.05 en ambas replicas e IC del mismo signo.
#
# Ajuste solo por secuencia: cada caracteristica de contexto (las 15, incluidas
# las dos redundantes, porque entran de a una) en un modelo propio con todas
# las caracteristicas de secuencia del summary (grupo "seq" de feature_groups:
# TATA, CCAAT, GC-box, isla CpG, YY1, CRE, E-box, poli(dA:dT), dinucleotido e
# INR del TSS, TCT, TEs, alto G+C) + razon CpG o/e y G+C continuas, sin las
# demas caracteristicas de contexto. Mismos datos y codificacion que el modelo
# principal (indicador de faltante de la propia caracteristica si tiene NA).
# Se reporta junto al crudo y al ajustado completo; delta R2 = R2 con la
# caracteristica - R2 solo secuencia (+ su indicador).
#
# Ademas: distribucion de la razon CpG o/e de los fragmentos de la library
# (todos los promotores EPD de prom_df, no solo high-conf) y prueba de
# bimodalidad (mezcla gaussiana de 1 vs. 2 componentes por BIC, D de Ashman y
# modos de la densidad kernel).
#
# Salidas (figures/exploratorio_contexto_ajustado/, sin fig_dir: no es figura
# de la tesis):
#   ajustado_efectos.tsv        crudo y ajustado por caracteristica y replica
#                               (principal y casos completos)
#   ajustado_covariables.tsv    coeficientes de CpG o/e, G+C e indicadores
#   ajustado_r2.tsv             R2 (completo, solo CpG o/e + G+C, solo
#                               caracteristicas + indicadores) y n por
#                               replica y modelo
#   ajustado_vif.tsv            VIF generalizado de cada termino
#   ajustado_crudo_vs_ajustado.jpg
#   secuencia_efectos.tsv       crudo / ajustado por secuencia / ajustado
#                               completo, por caracteristica y replica
#   secuencia_r2.tsv            R2 del modelo solo secuencia y VIF maximo
#   secuencia_crudo_vs_ajustado.jpg
# Ademas guarda el panel B de la figura de tesis summary_features_endogeno
# (fig_dir("summary_features_endogeno"), summary_endo_ajustado_secuencia.jpg/.pdf):
# crudo vs. ajustado por secuencia, sin el ajustado completo, y arma la
# figura completa summary_features_endogeno.jpg/.pdf (A + B). El panel A
# (summary_endo) lo genera R/15_summary_features/summary_features.R, que
# hay que correr antes (deja data/processed/panel_summary_endo.rds).
#   cpg_oe_distribucion.tsv, cpg_oe_distribucion.jpg
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(patchwork)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

out_dir <- "figures/exploratorio_contexto_ajustado"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

EXCLUIDAS <- c("Alta especificidad tisular", "Promotores angostos")

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  mutate(percentil = 100 * (rank(mean, ties.method = "average") - 0.5) / n(), .by = rep)

endo_features <- feature_groups$feature[feature_groups$group == "endo"]
features <- setdiff(endo_features, EXCLUIDAS)
seq_features <- feature_groups$feature[feature_groups$group == "seq"]

# build_tidy_features() conserva el orden de filas de `data`, asi que se pueden
# pegar columnas de `data` por posicion.
tidy <- build_tidy_features(data, keep = c("percentil", "seq_id")) %>%
  mutate(across(all_of(c(endo_features, seq_features)), ~ as.numeric(.x == "TRUE"))) %>%
  mutate(cpg_oe = data$cpg_oe, g_c = data$g_c)
stopifnot(identical(tidy$seq_id, data$seq_id))

nd <- "No detectado (FANTOM5)"
faltantes <- c(
  "Insertado en humanos" = "Sin dato: turnover humano-ratón",
  "Promotores anchos" = "Sin dato: forma del promotor",
  "Promotor unidireccional" = "Sin dato: orientación (PRO-cap)"
)

tidy <- tidy %>%
  # especificidad tisular y forma no definidas en los no detectados
  mutate(across(
    c("Baja especificidad tisular", "Alta especificidad tisular", "Promotores anchos", "Promotores angostos"),
    ~ if_else(.data[[nd]] %in% 1, 0, .x)
  ))

# Modelo principal: NA -> 0 + indicador de faltante
datos_ind <- tidy %>%
  filter(!is.na(.data[[nd]]), if_all(starts_with("Alta conservación"), ~ !is.na(.x))) %>%
  mutate(
    !!faltantes[[1]] := as.numeric(is.na(`Insertado en humanos`)),
    !!faltantes[[2]] := as.numeric(is.na(`Promotores anchos`)),
    !!faltantes[[3]] := as.numeric(is.na(`Promotor unidireccional`)),
    across(all_of(c(names(faltantes), "Promotores angostos")), ~ replace_na(.x, 0))
  )
# Sensibilidad: casos completos
datos_cc <- tidy %>% filter(if_all(all_of(features), ~ !is.na(.x)))

# Covariables continuas estandarizadas dentro de cada replica (coeficiente =
# puntos de percentil por 1 DE).
estandarizar <- function(df) {
  df %>% mutate(across(c(cpg_oe, g_c), ~ (.x - mean(.x)) / sd(.x)), .by = rep)
}
datos_ind <- estandarizar(datos_ind)
datos_cc <- estandarizar(datos_cc)

bt <- function(x) paste0("`", x, "`")
formula_de <- function(terms) as.formula(paste("percentil ~", paste(bt(terms), collapse = " + ")))

coefs <- function(fit, terms) {
  stopifnot(!anyNA(coef(fit)))
  ci <- confint(fit)
  s <- summary(fit)$coefficients
  # lm solo pone backticks en los nombres no sintacticos
  bt <- function(x) ifelse(make.names(x) == x, x, paste0("`", x, "`"))
  tibble(
    termino = terms,
    estimate = s[bt(terms), "Estimate"],
    IC2.5 = ci[bt(terms), 1],
    IC97.5 = ci[bt(terms), 2],
    pval = s[bt(terms), "Pr(>|t|)"]
  )
}

ajustar <- function(df, modelo, extra) {
  map(c("Rep 1", "Rep 2"), function(rep_id) {
    d <- df %>% filter(rep == rep_id)
    terms_full <- c(features, extra, "cpg_oe", "g_c")
    fit <- lm(formula_de(terms_full), data = d)
    fit_cov <- lm(formula_de(c("cpg_oe", "g_c")), data = d)
    fit_feat <- lm(formula_de(c(features, extra)), data = d)

    crudo <- map(features, ~ coefs(lm(formula_de(.x), data = d), .x)) %>% list_rbind()
    ajust <- coefs(fit, features)
    efectos <- inner_join(crudo, ajust, by = "termino", suffix = c("_crudo", "_ajustado")) %>%
      mutate(
        prevalencia = map_dbl(termino, ~ mean(d[[.x]])),
        pBH_crudo = p.adjust(pval_crudo, "BH"),
        pBH_ajustado = p.adjust(pval_ajustado, "BH")
      )

    cov <- coefs(fit, c(extra, "cpg_oe", "g_c"))

    vif <- car::vif(fit)
    vif <- tibble(termino = str_remove_all(names(vif), "`"), vif = unname(vif))

    r2 <- tibble(
      modelo_r2 = c("completo", "solo CpG o/e + G+C", "solo características (+ indicadores)"),
      r2 = c(summary(fit)$r.squared, summary(fit_cov)$r.squared, summary(fit_feat)$r.squared),
      r2_ajustado = c(summary(fit)$adj.r.squared, summary(fit_cov)$adj.r.squared, summary(fit_feat)$adj.r.squared),
      n = nobs(fit)
    )
    list(efectos = efectos, cov = cov, vif = vif, r2 = r2) %>% map(~ mutate(.x, rep = rep_id, modelo = modelo))
  }) %>%
    transpose() %>%
    map(list_rbind)
}

res <- list(
  ajustar(datos_ind, "principal (indicadores de faltante)", unname(faltantes)),
  ajustar(datos_cc, "casos completos", character(0))
) %>%
  transpose() %>%
  map(list_rbind)

sobrevive <- function(df, sufijo) {
  df %>%
    group_by(modelo, termino) %>%
    mutate(
      "sobrevive_{sufijo}" := all(.data[[paste0("pBH_", sufijo)]] < 0.05) &
        n_distinct(sign(c(.data[[paste0("IC2.5_", sufijo)]], .data[[paste0("IC97.5_", sufijo)]]))) == 1
    ) %>%
    ungroup()
}

efectos <- res$efectos %>%
  sobrevive("crudo") %>%
  sobrevive("ajustado") %>%
  select(modelo, rep, caracteristica = termino, prevalencia,
         ends_with("_crudo"), ends_with("_ajustado"))

write_tsv(efectos, file.path(out_dir, "ajustado_efectos.tsv"))
write_tsv(res$cov %>% rename(covariable = termino) %>% select(modelo, rep, everything()), file.path(out_dir, "ajustado_covariables.tsv"))
write_tsv(res$r2 %>% select(modelo, rep, everything()), file.path(out_dir, "ajustado_r2.tsv"))
write_tsv(res$vif %>% select(modelo, rep, everything()), file.path(out_dir, "ajustado_vif.tsv"))

# --- Plot crudo vs. ajustado (modelo principal) --------------------------------
# Por caracteristica: crudo (hueco) -> ajustado (relleno, con IC 95%).

principal <- efectos %>% filter(modelo == "principal (indicadores de faltante)")
orden <- principal %>%
  filter(rep == "Rep 1") %>%
  arrange(estimate_crudo) %>%
  pull(caracteristica)

plot_cva <- principal %>%
  mutate(y = match(caracteristica, orden)) %>%
  ggplot(aes(y = y)) +
  geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.5) +
  geom_segment(aes(x = estimate_crudo, xend = estimate_ajustado, yend = y), colour = "#bbbbbb", linewidth = 0.6) +
  geom_segment(aes(x = IC2.5_ajustado, xend = IC97.5_ajustado, y = y - 0.18, yend = y - 0.18), colour = "#2a78d6", linewidth = 0.7) +
  geom_segment(aes(x = IC2.5_crudo, xend = IC97.5_crudo, y = y + 0.18, yend = y + 0.18), colour = "#5f5e5a", linewidth = 0.7) +
  geom_point(aes(x = estimate_crudo, y = y + 0.18, colour = "Crudo", fill = ifelse(sobrevive_crudo, "Crudo", "no")), shape = 21, size = 2.6, stroke = 1) +
  geom_point(aes(x = estimate_ajustado, y = y - 0.18, colour = "Ajustado", fill = ifelse(sobrevive_ajustado, "Ajustado", "no")), shape = 21, size = 2.6, stroke = 1) +
  facet_wrap(~rep, nrow = 1) +
  scale_y_continuous(breaks = seq_along(orden), labels = orden, expand = expansion(add = 0.5)) +
  scale_colour_manual(values = c("Crudo" = "#5f5e5a", "Ajustado" = "#2a78d6"), breaks = c("Crudo", "Ajustado")) +
  scale_fill_manual(values = c("Crudo" = "#5f5e5a", "Ajustado" = "#2a78d6", "no" = "white"), guide = "none") +
  theme_pubr(base_size = 14) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "#eeeeee"), strip.background = element_blank()) +
  labs(
    x = "Efecto sobre el percentil de actividad (puntos, IC 95%)", y = NULL, colour = NULL,
    caption = paste0(
      "Ajustado: todas las características juntas + razón CpG o/e y G+C (continuas) + indicadores de dato faltante.\n",
      "Punto relleno: p BH < 0,05 en ambas réplicas e IC del mismo signo; hueco: no."
    )
  )

ggsave(file.path(out_dir, "ajustado_crudo_vs_ajustado.jpg"), plot_cva, width = 13, height = 8, units = "in", dpi = 300)

# --- Distribucion de la razon CpG o/e ------------------------------------------

oe <- prom_df$cpg_oe[!is.na(prom_df$cpg_oe)]

# Mezcla gaussiana de 1 vs. 2 componentes (EM).
em_2normales <- function(x, iter = 1000, tol = 1e-8) {
  q <- quantile(x, c(0.25, 0.75))
  mu <- unname(q); s <- rep(sd(x) / 2, 2); w <- c(0.5, 0.5)
  ll_prev <- -Inf
  for (i in seq_len(iter)) {
    d1 <- w[1] * dnorm(x, mu[1], s[1]); d2 <- w[2] * dnorm(x, mu[2], s[2])
    ll <- sum(log(d1 + d2))
    r <- d1 / (d1 + d2)
    w <- c(mean(r), 1 - mean(r))
    mu <- c(sum(r * x) / sum(r), sum((1 - r) * x) / sum(1 - r))
    s <- sqrt(c(sum(r * (x - mu[1])^2) / sum(r), sum((1 - r) * (x - mu[2])^2) / sum(1 - r)))
    if (abs(ll - ll_prev) < tol) break
    ll_prev <- ll
  }
  list(w = w, mu = mu, s = s, loglik = ll)
}
m2 <- em_2normales(oe)
ll1 <- sum(dnorm(oe, mean(oe), sd(oe) * sqrt((length(oe) - 1) / length(oe)), log = TRUE))
bic1 <- -2 * ll1 + 2 * log(length(oe))
bic2 <- -2 * m2$loglik + 5 * log(length(oe))
ashman_d <- sqrt(2) * abs(diff(m2$mu)) / sqrt(sum(m2$s^2))

dens <- density(oe)
es_max <- which(diff(sign(diff(dens$y))) == -2) + 1
es_min <- which(diff(sign(diff(dens$y))) == 2) + 1
modos <- dens$x[es_max][dens$y[es_max] > 0.05 * max(dens$y)]
antimodo <- dens$x[es_min][dens$x[es_min] > min(modos) & dens$x[es_min] < max(modos)]

bimodal <- tibble(
  n = length(oe),
  media = mean(oe), mediana = median(oe),
  n_modos_kde = length(modos),
  modos_kde = paste(round(modos, 3), collapse = "; "),
  antimodo_kde = paste(round(antimodo, 3), collapse = "; "),
  bic_1_normal = bic1, bic_2_normales = bic2, delta_bic = bic1 - bic2,
  mezcla_peso = paste(round(m2$w, 3), collapse = "; "),
  mezcla_media = paste(round(m2$mu, 3), collapse = "; "),
  mezcla_sd = paste(round(m2$s, 3), collapse = "; "),
  ashman_d = ashman_d,
  frac_oe_mayor_0.6 = mean(oe > 0.6),
  frac_cgi_frag = mean(prom_df$CGI_frag, na.rm = TRUE)
)
write_tsv(bimodal, file.path(out_dir, "cpg_oe_distribucion.tsv"))

curvas <- tibble(x = seq(min(oe), max(oe), length.out = 500)) %>%
  mutate(
    `Componente 1` = m2$w[1] * dnorm(x, m2$mu[1], m2$s[1]),
    `Componente 2` = m2$w[2] * dnorm(x, m2$mu[2], m2$s[2])
  ) %>%
  pivot_longer(-x, names_to = "componente", values_to = "densidad")

plot_oe <- ggplot(tibble(oe = oe), aes(oe)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.025, fill = "#cde2fb", colour = "white", linewidth = 0.2) +
  geom_density(colour = "#184f95", linewidth = 0.8) +
  geom_line(data = curvas, aes(x, densidad, linetype = componente), colour = "#5f5e5a", linewidth = 0.5) +
  geom_vline(xintercept = 0.6, colour = "#eb6834", linewidth = 0.6, linetype = "dashed") +
  annotate("text", x = 0.61, y = Inf, vjust = 1.5, hjust = 0, label = "Corte de isla CpG (o/e > 0,6)", colour = "#eb6834", size = 4) +
  scale_linetype_manual(values = c("Componente 1" = "dotted", "Componente 2" = "dotdash")) +
  theme_pubr(base_size = 14) +
  theme(legend.position = "top") +
  labs(
    x = "Razón CpG observado/esperado del fragmento", y = "Densidad", linetype = "Mezcla de 2 normales",
    caption = sprintf(
      "n = %s promotores de la library. Modos (KDE): %s; antimodo: %s. ΔBIC (1 vs. 2 normales) = %.0f; D de Ashman = %.2f.",
      format(length(oe), big.mark = "."), bimodal$modos_kde, bimodal$antimodo_kde, bimodal$delta_bic, ashman_d
    )
  )
ggsave(file.path(out_dir, "cpg_oe_distribucion.jpg"), plot_oe, width = 10, height = 6, units = "in", dpi = 300)

# --- Ajuste solo por secuencia ---------------------------------------------------

seq_terms <- c(seq_features, "cpg_oe", "g_c")
indicador_de <- c(faltantes, "Promotores angostos" = "Sin dato: forma del promotor")

ajuste_seq <- map(c("Rep 1", "Rep 2"), function(rep_id) {
  d <- datos_ind %>% filter(rep == rep_id)
  base <- lm(formula_de(seq_terms), data = d)
  vif_base <- car::vif(base)
  efectos_rep <- map(endo_features, function(f) {
    ind <- if (f %in% names(indicador_de)) unname(indicador_de[f]) else character(0)
    fit <- lm(formula_de(c(f, ind, seq_terms)), data = d)
    base_f <- lm(formula_de(c(ind, seq_terms)), data = d)
    crudo <- coefs(lm(formula_de(f), data = d), f)
    coefs(fit, f) %>%
      rename_with(~ paste0(.x, "_secuencia"), -termino) %>%
      left_join(crudo %>% rename_with(~ paste0(.x, "_crudo"), -termino), by = "termino") %>%
      mutate(prevalencia = mean(d[[f]]), delta_r2 = summary(fit)$r.squared - summary(base_f)$r.squared)
  }) %>%
    list_rbind() %>%
    mutate(pBH_crudo = p.adjust(pval_crudo, "BH"), pBH_secuencia = p.adjust(pval_secuencia, "BH"), rep = rep_id)
  list(
    efectos = efectos_rep,
    r2 = tibble(rep = rep_id, n = nobs(base), r2_solo_secuencia = summary(base)$r.squared,
                r2_ajustado_solo_secuencia = summary(base)$adj.r.squared,
                vif_max = max(vif_base), vif_max_termino = str_remove_all(names(which.max(vif_base)), "`"))
  )
}) %>%
  transpose() %>%
  map(list_rbind)

efectos_seq <- ajuste_seq$efectos %>%
  mutate(modelo = "secuencia") %>%
  sobrevive("crudo") %>%
  sobrevive("secuencia") %>%
  select(-modelo) %>%
  left_join(
    principal %>% select(rep, termino = caracteristica, ends_with("_ajustado")) %>%
      rename_with(~ str_replace(.x, "_ajustado$", "_completo")),
    by = c("rep", "termino")
  ) %>%
  select(rep, caracteristica = termino, prevalencia, ends_with("_crudo"), ends_with("_secuencia"), delta_r2, ends_with("_completo"))

write_tsv(efectos_seq, file.path(out_dir, "secuencia_efectos.tsv"))
write_tsv(ajuste_seq$r2, file.path(out_dir, "secuencia_r2.tsv"))

# Plot: crudo / ajustado por secuencia / ajustado completo (las dos excluidas
# del modelo completo no tienen ese punto).
niveles_seq <- c("Crudo", "Ajustado por secuencia", "Ajustado completo")
colores_seq <- c("Crudo" = "#5f5e5a", "Ajustado por secuencia" = "#eb6834", "Ajustado completo" = "#2a78d6")
desplaz_seq <- c("Crudo" = 0.25, "Ajustado por secuencia" = 0, "Ajustado completo" = -0.25)
orden_seq <- efectos_seq %>%
  filter(rep == "Rep 1") %>%
  arrange(estimate_crudo) %>%
  pull(caracteristica)

plot_seq <- efectos_seq %>%
  pivot_longer(
    c(matches("^(estimate|IC2.5|IC97.5|sobrevive)_(crudo|secuencia|completo)$")),
    names_to = c(".value", "tipo"), names_pattern = "^(.+)_(crudo|secuencia|completo)$"
  ) %>%
  filter(!is.na(estimate)) %>%
  mutate(
    tipo = recode(tipo, crudo = "Crudo", secuencia = "Ajustado por secuencia", completo = "Ajustado completo"),
    tipo = factor(tipo, levels = niveles_seq),
    y = match(caracteristica, orden_seq) + desplaz_seq[as.character(tipo)],
    relleno = ifelse(sobrevive, as.character(tipo), "no")
  ) %>%
  ggplot(aes(x = estimate, y = y, colour = tipo)) +
  geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.5) +
  geom_segment(aes(x = IC2.5, xend = IC97.5, yend = y), linewidth = 0.7) +
  geom_point(aes(fill = relleno), shape = 21, size = 2.4, stroke = 1) +
  facet_wrap(~rep, nrow = 1) +
  scale_y_continuous(breaks = seq_along(orden_seq), labels = orden_seq, expand = expansion(add = 0.5)) +
  scale_colour_manual(values = colores_seq, breaks = niveles_seq) +
  scale_fill_manual(values = c(colores_seq, "no" = "white"), guide = "none") +
  theme_pubr(base_size = 14) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "#eeeeee"), strip.background = element_blank()) +
  labs(
    x = "Efecto sobre el percentil de actividad (puntos, IC 95%)", y = NULL, colour = NULL,
    caption = paste0(
      "Ajustado por secuencia: la característica sola + todas las de secuencia del summary + razón CpG o/e y G+C.
",
      "Ajustado completo: todas las de contexto juntas + CpG o/e y G+C (sin alta especificidad tisular ni promotores angostos).
",
      "Punto relleno: p BH < 0,05 en ambas réplicas e IC del mismo signo; hueco: no."
    )
  )
ggsave(file.path(out_dir, "secuencia_crudo_vs_ajustado.jpg"), plot_seq, width = 13, height = 9, units = "in", dpi = 300)

# Panel B de summary_features_endogeno: solo crudo y ajustado por secuencia.
# Relleno = p BH < 0.05 en ambas replicas e IC del mismo signo; hueco = no.
niveles_b <- c("Crudo", "Ajustado por secuencia")
colores_b <- colores_seq[niveles_b]
desplaz_b <- c("Crudo" = 0.18, "Ajustado por secuencia" = -0.18)

panel_b <- efectos_seq %>%
  pivot_longer(
    c(matches("^(estimate|IC2.5|IC97.5|sobrevive)_(crudo|secuencia)$")),
    names_to = c(".value", "tipo"), names_pattern = "^(.+)_(crudo|secuencia)$"
  ) %>%
  mutate(
    tipo = factor(recode(tipo, crudo = "Crudo", secuencia = "Ajustado por secuencia"), levels = niveles_b),
    y = match(caracteristica, orden_seq) + desplaz_b[as.character(tipo)],
    relleno = ifelse(sobrevive, as.character(tipo), "no")
  ) %>%
  ggplot(aes(x = estimate, y = y, colour = tipo)) +
  geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.5) +
  geom_segment(aes(x = IC2.5, xend = IC97.5, yend = y), linewidth = 0.8) +
  geom_point(aes(fill = relleno), shape = 21, size = 3, stroke = 1.1) +
  facet_wrap(~rep, nrow = 1) +
  scale_y_continuous(breaks = seq_along(orden_seq), labels = orden_seq, expand = expansion(add = 0.5)) +
  scale_colour_manual(values = colores_b, breaks = niveles_b) +
  scale_fill_manual(values = c(colores_b, "no" = "white"), guide = "none") +
  theme_pubr(base_size = 20) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "#eeeeee"), strip.background = element_blank()) +
  labs(x = "Efecto sobre el percentil de actividad", y = NULL, colour = NULL)

dir_endo <- fig_dir("summary_features_endogeno")
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(dir_endo, paste0("summary_endo_ajustado_secuencia.", ext)), panel_b, width = 13.5, height = 8, units = "in", dpi = 300)
}

# Figura completa de summary_features_endogeno: A = summary_endo (lo genera
# R/15_summary_features/summary_features.R y lo deja en
# data/processed/panel_summary_endo.rds; correr ese script antes), B = panel_b.
panel_a <- readRDS("data/processed/panel_summary_endo.rds")
figura_endo <- wrap_elements(full = panel_a) / wrap_elements(full = panel_b) +
  plot_layout(heights = c(6.75, 8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30, face = "bold"))
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(dir_endo, paste0("summary_features_endogeno.", ext)), figura_endo, width = 13.5, height = 15.5, units = "in", dpi = 300)
}

print(as.data.frame(ajuste_seq$r2))
print(as.data.frame(efectos_seq %>% select(rep, caracteristica, estimate_crudo, estimate_secuencia, IC2.5_secuencia, IC97.5_secuencia, pBH_secuencia, sobrevive_secuencia, delta_r2, estimate_completo) %>% mutate(across(where(is.double), ~ signif(.x, 3)))))
print(as.data.frame(bimodal))
print(as.data.frame(res$r2))
options(width = 250)
print(as.data.frame(efectos %>% mutate(across(where(is.double), ~ signif(.x, 3)))))
message("Salidas guardadas en ", out_dir)
