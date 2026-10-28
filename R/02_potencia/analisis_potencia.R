# analisis_potencia (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Analisis de potencia por simulacion para elegir el numero de fracciones
# de sorting (K). Usa solo informacion disponible antes del sorting: el
# perfil de citometria previo de la library (que define el rango a
# fraccionar) y parametros de diseno (celulas por fraccion, n de
# promotores, replicas, profundidad planificada). NO usa los conteos de
# amplicon-seq del experimento.
#
# Modelo:
# - Expresion por celula de cada promotor ~ lognormal: log10 EGFP ~
#   N(mu, sigma). sigma = ruido del promotor; se suma un error de
#   medicion/clasificacion del sorter (tau) en la asignacion a gates.
# - mu de los promotores: perfil de la library (celulas dentro del rango
#   fraccionado) contraido hacia su media para descontar sigma0 (ajuste
#   por momentos: la mezcla de N(mu, sigma0) reproduce la varianza del
#   perfil).
# - Fracciones: desde el umbral inferior, K gates de igual ancho en log10
#   que cubren siempre el mismo tramo del perfil (el de las 7 fracciones del
#   sorter: 6 gates cerrados de los breaks de gating_facs + uno mas del mismo
#   ancho); el ultimo gate queda abierto hacia arriba. Asi K = 7 reproduce
#   los gates reales (la fraccion abierta = 1% mas brillante del perfil) y
#   con otro K la fraccion abierta tiene el contenido que le corresponde por
#   su ancho (mas del 1% con K < 7, menos con K > 7). Por debajo del umbral
#   no se recolecta.
# - Muestreo: de cada fraccion se recolectan n_cel celulas; las celulas
#   de un promotor en la fraccion k ~ Poisson(n_cel * f_p * P_p(k) /
#   P_lib(k)), con f_p la frecuencia del promotor en la library
#   (lognormal alrededor de 1/N) y P_lib(k) la fraccion de la library en
#   el gate k (del perfil). Lecturas ~ Poisson(celulas * lecturas/celula).
# - Estimacion: igual que el pipeline (R/01_activity_stats): conteos
#   normalizados por esfuerzo de muestreo (P_lib(k)), media y varianza
#   del indice de gate; se pasan a log10 multiplicando por el ancho de
#   gate para comparar entre K.
# - Potencia: pares de promotores (3 replicas c/u) con diferencia de
#   media (dmu, igual sigma) o de ruido (dsigma, igual mu). Estadistico =
#   diferencia estimada promediada entre replicas; valor critico (alfa =
#   0.05, bilateral) calibrado con pares nulos simulados, por estrato de
#   mu. Efecto minimo detectable (EMD) = efecto con potencia 0.8.
#
# Paneles:
# A: un promotor de ejemplo discretizado con K = 3, 7 y 12.
# B: EMD (log10) vs. K, con K = 6-7 resaltado; filas = media / ruido,
#    columnas = 250000 celulas por fraccion (diseno) / mismo total de
#    celulas que K = 7 (aisla el efecto de la resolucion), color = sigma0
#    (0.15/0.25/0.35).
#
# Requiere: data/raw/Citometry_18-09-2023/library 1.fcs (se cachea el
# perfil en data/processed/potencia_perfil_library_v1.tsv; la simulacion
# en data/processed/potencia_sim_v4.rds). Run from the TesisDoc repo root.

library(tidyverse)
library(patchwork)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "analisis_potencia"
out_dir <- fig_dir(slug)
set.seed(20230918)

# --- Parametros de diseno --------------------------------------------------

# rango dinamico fraccionado: el de los breaks del sorter (gating_facs),
# ~1.04 decadas; su posicion se fija sobre el perfil previo (ver abajo)
breaks_sorter <- c(5.56, 8.35, 12.53, 17.98, 27, 40.54, 60.86)
rango_log <- log10(max(breaks_sorter) / min(breaks_sorter))
K_grid <- 3:12
n_cel <- 250000 # celulas recolectadas por fraccion y replica (diseno)
# dos escenarios de esfuerzo: "por_fraccion" = n_cel en cada fraccion (el
# total de celulas, y el tiempo de sorting, crecen con K); "total_fijo" =
# el total de K = 7 (7 x n_cel) repartido entre las K fracciones, para
# separar el efecto de la resolucion del de recolectar mas celulas
escenarios <- c("por_fraccion", "total_fijo")
n_cel_K <- function(K, escenario) if (escenario == "por_fraccion") n_cel else n_cel * 7 / K
n_prom <- 23806 # promotores + enhancers + cancer related de la library
sd_repr <- 0.5 # dispersion (log natural) de la representacion de cada promotor
lect_por_cel <- 20 # profundidad planificada (lecturas por celula recolectada)
tau <- 0.05 # error de medicion/clasificacion del sorter (log10)
n_rep <- 3
min_cel <- 20 # minimo de celulas recolectadas por replica para estimar (analogo al umbral de conteos)
sigma0_grid <- c(0.15, 0.25, 0.35)
sigma0_ref <- 0.25
dmu_grid <- seq(0, 0.3, by = 0.025)
dsigma_grid <- seq(0, 0.2, by = 0.02)
n_sim <- 5000 # pares por condicion
alfa <- 0.05
n_estratos <- 5 # estratos de mu para calibrar el valor critico

# --- Perfil de citometria previo (library) ---------------------------------

perfil_path <- "data/processed/potencia_perfil_library_v1.tsv"
if (!file.exists(perfil_path)) {
  library(flowCore)
  library(flowWorkspace)
  filter <- dplyr::filter
  select <- dplyr::select
  # mismo gating de outliers y singlets que gating_facs
  fs <- read.flowSet(files = "data/raw/Citometry_18-09-2023/library 1.fcs", transformation = FALSE,
    emptyValue = FALSE, alter.names = TRUE, truncate_max_range = FALSE)
  gs <- GatingSet(fs)
  gs_pop_add(gs, rectangleGate(filterId = "-outlier", "FSC1.Height" = c(.4e9, 3.5e9), "SSC1.Height" = c(5e7, 1e9)), parent = "root")
  gs_pop_add(gs, rectangleGate(filterId = "singlet", "FSC1.Height" = c(0, 1e10), "FSC1.Width" = c(0, 2e9)), parent = "-outlier")
  recompute(gs)
  egfp <- exprs(gs_pop_get_data(gs, "singlet")[[1]])[, "FL18.Height"]
  write_tsv(tibble(log10_egfp = log10(egfp[egfp > 0])), perfil_path)
  filter <- dplyr::filter
  select <- dplyr::select
}
perfil <- read_tsv(perfil_path, show_col_types = FALSE)$log10_egfp

# Como en el diseno real: la fraccion superior abierta recolecta el 1% mas
# brillante del perfil previo, y por debajo se cubre el mismo rango
# dinamico que los breaks del sorter. (La escala del citometro analizador
# no coincide exactamente con la del sorter, por eso los limites se ubican
# sobre el perfil y no trasladando los breaks.)
lim_sup <- unname(quantile(perfil, 0.99))
lim_inf <- lim_sup - rango_log
message("Rango fraccionado (log10): ", round(lim_inf, 2), " - ", round(lim_sup, 2),
  "; fraccion de la library recolectada: ", round(mean(perfil >= lim_inf), 3))

# --- Fracciones para cada K --------------------------------------------------

# K gates de igual ancho sobre el tramo de las 7 fracciones del sorter
# (rango_log + un ancho de gate); el ultimo, abierto. Con K = 7 los limites
# coinciden con breaks_sorter y la fraccion abierta empieza en lim_sup.
tramo <- rango_log * 7 / 6
ancho_K <- function(K) tramo / K
limites_K <- function(K) c(lim_inf + ancho_K(K) * 0:(K - 1), Inf)
message("Ancho de gate con K = 7: ", round(ancho_K(7), 3), " log10 (sorter: ",
  paste(round(diff(log10(breaks_sorter)), 3), collapse = ", "), ")")

# fraccion de la poblacion total (library) en cada gate
p_lib_K <- function(K) as.numeric(table(cut(perfil, limites_K(K), right = FALSE))) / length(perfil)

# --- Distribucion de mu de los promotores ------------------------------------

perfil_ord <- perfil[perfil >= lim_inf]
m_perfil <- mean(perfil_ord)
v_perfil <- var(perfil_ord)
sample_mu <- function(n, sigma0 = sigma0_ref) {
  x <- sample(perfil_ord, n, replace = TRUE)
  m_perfil + sqrt(max(1 - sigma0^2 / v_perfil, 0)) * (x - m_perfil)
}

# --- Simulacion de un set de promotores en una replica ------------------------

# mu, sigma, repr: vectores (un elemento por promotor). Devuelve media y sd
# estimadas (log10) y si pasa el minimo de celulas.
simular_medicion <- function(mu, sigma, repr, K, n_cel_frac = n_cel) {
  lim <- limites_K(K)
  p_lib <- p_lib_K(K)
  s_eff <- sqrt(sigma^2 + tau^2)
  # P_p(k): probabilidad de que una celula del promotor caiga en el gate k
  cdf <- sapply(lim, function(b) pnorm(b, mu, s_eff))
  p_prom <- cdf[, -1, drop = FALSE] - cdf[, -ncol(cdf), drop = FALSE]
  lambda <- n_cel_frac * (repr / n_prom) * sweep(p_prom, 2, p_lib, "/")
  celulas <- matrix(rpois(length(lambda), lambda), nrow = nrow(lambda))
  lecturas <- matrix(rpois(length(celulas), celulas * lect_por_cel), nrow = nrow(celulas))
  # normalizacion por esfuerzo de muestreo, como perc_cells en el pipeline
  w <- sweep(lecturas, 2, p_lib / min(p_lib), "*")
  idx <- seq_len(K)
  tot <- rowSums(w)
  media_idx <- as.numeric(w %*% idx) / tot
  var_idx <- as.numeric(w %*% idx^2) / tot - media_idx^2
  tibble(
    media = media_idx * ancho_K(K),
    sd = sqrt(pmax(var_idx, 0)) * ancho_K(K),
    ok = rowSums(celulas) >= min_cel
  )
}

# Pares: promotor 1 (mu1, sigma0) y promotor 2 (mu1 + dmu, sigma0 + dsigma).
# Devuelve la diferencia estimada (media entre replicas) en media y en sd.
simular_pares <- function(K, sigma0, dmu, dsigma, n, n_cel_frac) {
  mu1 <- sample_mu(n, sigma0)
  mu2 <- mu1 + dmu
  repr1 <- exp(rnorm(n, 0, sd_repr))
  repr2 <- exp(rnorm(n, 0, sd_repr))
  reps <- map(seq_len(n_rep), function(r) {
    a <- simular_medicion(mu1, rep(sigma0, n), repr1, K, n_cel_frac)
    b <- simular_medicion(mu2, rep(sigma0 + dsigma, n), repr2, K, n_cel_frac)
    tibble(d_media = b$media - a$media, d_sd = b$sd - a$sd, ok = a$ok & b$ok)
  })
  ok <- reduce(map(reps, "ok"), `&`)
  tibble(
    mu1 = mu1,
    d_media = reduce(map(reps, "d_media"), `+`) / n_rep,
    d_sd = reduce(map(reps, "d_sd"), `+`) / n_rep
  )[ok, ]
}

estratos_mu <- quantile(sample_mu(1e5), seq(0, 1, length.out = n_estratos + 1))
estratos_mu[c(1, length(estratos_mu))] <- c(-Inf, Inf)

# potencia de un tipo de efecto ("media" o "ruido") para un K, sigma0 y escenario
potencia_K <- function(K, sigma0, efecto, escenario) {
  grid <- if (efecto == "media") dmu_grid else dsigma_grid
  col <- if (efecto == "media") "d_media" else "d_sd"
  sims <- map(grid, function(d) {
    simular_pares(K, sigma0, dmu = if (efecto == "media") d else 0,
      dsigma = if (efecto == "ruido") d else 0, n = n_sim, n_cel_frac = n_cel_K(K, escenario)) %>%
      mutate(efecto_tam = d, estrato = cut(mu1, estratos_mu, labels = FALSE))
  }) %>% list_rbind()
  # valor critico por estrato de mu, a partir de los pares nulos (efecto 0)
  critico <- sims %>%
    filter(efecto_tam == 0) %>%
    group_by(estrato) %>%
    summarise(c_val = quantile(abs(.data[[col]]), 1 - alfa), .groups = "drop")
  sims %>%
    left_join(critico, by = "estrato") %>%
    group_by(efecto_tam) %>%
    summarise(potencia = mean(abs(.data[[col]]) > c_val), n = n(), .groups = "drop") %>%
    mutate(K = K, sigma0 = sigma0, efecto = efecto, escenario = escenario)
}

sim_path <- "data/processed/potencia_sim_v4.rds"
if (!file.exists(sim_path)) {
  t0 <- Sys.time()
  potencia <- expand_grid(K = K_grid, sigma0 = sigma0_grid, efecto = c("media", "ruido"), escenario = escenarios) %>%
    pmap(potencia_K) %>%
    list_rbind()
  message("Simulacion: ", round(as.numeric(Sys.time() - t0, units = "secs")), "s")
  saveRDS(potencia, sim_path)
}
potencia <- readRDS(sim_path)

# EMD: efecto con potencia 0.8 segun un ajuste probit de la curva de
# potencia (menos sensible al ruido de Monte Carlo que interpolar la grilla)
emd <- potencia %>%
  group_by(K, sigma0, efecto, escenario) %>%
  summarise(emd = {
    # sin los puntos saturados (potencia ~1), que no aportan y no dejan converger el ajuste
    ok <- potencia < 0.995
    efecto_tam <- efecto_tam[ok]
    n <- n[ok]
    potencia <- potencia[ok]
    fit <- glm(cbind(round(potencia * n), n - round(potencia * n)) ~ efecto_tam, family = binomial("probit"))
    unname((qnorm(0.8) - coef(fit)[1]) / coef(fit)[2])
  }, .groups = "drop")
write_tsv(emd, "tables/analisis_potencia_emd.tsv")
print(emd %>% pivot_wider(names_from = c(efecto, sigma0), values_from = emd) %>% arrange(escenario, K), n = Inf)
# --- Figura ------------------------------------------------------------------

base_txt <- 22
tema <- theme_pubclean(base_size = base_txt) + theme(legend.key.width = unit(2, "lines"))
etiqueta_K <- function(K) paste0("K = ", K, if_else(K == 7, " (diseño)", ""))
k_resaltado <- c(6, 7) # opciones evaluadas como "suficientes"

# A: promotor de ejemplo (mu en el medio del rango, sigma0_ref)
mu_ej <- (lim_inf + lim_sup) / 2
ej <- map(c(3, 7, 12), function(K) {
  lim <- limites_K(K)
  lim_plot <- pmin(lim, lim_sup + 0.35)
  tibble(K = K, gate = seq_len(K), x0 = head(lim_plot, -1), x1 = tail(lim_plot, -1),
    p = diff(pnorm(lim, mu_ej, sigma0_ref))) %>%
    mutate(alto = p / (x1 - x0))
}) %>% list_rbind() %>% mutate(K_lab = factor(etiqueta_K(K), levels = etiqueta_K(c(3, 7, 12))))
curva <- tibble(x = seq(lim_inf - 0.4, lim_sup + 0.35, length.out = 400)) %>%
  mutate(y = dnorm(x, mu_ej, sigma0_ref))
panel_a <- ggplot(ej) +
  geom_rect(aes(xmin = x0, xmax = x1, ymin = 0, ymax = alto), fill = "#AADAD4", col = "#216869") +
  geom_line(data = curva, aes(x, y), linewidth = 1, col = "grey20") +
  # umbral inferior y comienzo de la fraccion abierta de cada K
  geom_vline(data = ej %>% filter(gate == 1 | gate == K) %>% mutate(x = x0),
    aes(xintercept = x), linetype = "dashed", col = "grey40") +
  facet_wrap(~K_lab, nrow = 1) +
  labs(x = "log10 EGFP", y = "Densidad") +
  tema

# B: EMD (log10) vs. K; filas = tipo de efecto (escala Y libre), columnas =
# escenario de esfuerzo (el diseno real, y el mismo total de celulas que K = 7
# repartido entre las K fracciones, que aisla el efecto de la resolucion),
# color = ruido de base sigma0
lab_escenario <- c(por_fraccion = "250000 células por fracción (diseño)",
  total_fijo = "Mismo total de células que K = 7")
emd_plot <- emd %>%
  mutate(efecto = factor(recode(efecto, media = "Actividad media (Δμ)", ruido = "Ruido (Δσ)"),
    levels = c("Actividad media (Δμ)", "Ruido (Δσ)")),
    escenario = factor(lab_escenario[escenario], levels = lab_escenario))
clr_sigma0 <- c("0.15" = "#A8DADC", "0.25" = "#358AAA", "0.35" = "#0B2A3C")
panel_b <- ggplot(emd_plot, aes(K, emd, col = factor(sigma0))) +
  annotate("rect", xmin = min(k_resaltado) - 0.4, xmax = max(k_resaltado) + 0.4,
    ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = clr_sigma0, name = "Ruido de base (σ0)") +
  scale_x_continuous(breaks = K_grid) +
  expand_limits(y = 0) +
  facet_grid(efecto ~ escenario, scales = "free_y") +
  guides(col = guide_legend(override.aes = list(linewidth = 2, size = 4))) +
  labs(x = "Número de fracciones (K)", y = "Efecto mínimo detectable (log10)") +
  tema + theme(panel.spacing = unit(2, "lines"))

for (ext in c("jpg", "pdf")) {
  ggsave(file.path(out_dir, paste0("analisis_potencia.", ext)),
    (panel_a / panel_b) + plot_layout(heights = c(0.7, 1.5)) +
      plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 32, face = "bold")),
    width = 16, height = 16, units = "in", device = if (ext == "pdf") cairo_pdf else ragg::agg_jpeg)
}
message("Figura guardada en ", out_dir)
