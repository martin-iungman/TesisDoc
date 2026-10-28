# tata_ruido_especificidad (ver docs/mapping_figuras.csv para el numero de
# figura vigente) - Especificidad del efecto de la TATA-box sobre el ruido y
# su relacion con la especificidad tisular. Version de figura de los puntos
# p4, p7, p8 y p9 del exploratorio R/44_exploratorio_tata_ruido/tata_ruido.R
# (ver su handoff.md); se recalcula todo desde los datos.
#
# Metodo (igual que el exploratorio):
#   - Ruido = var_rank_sw (add_noise_rank() sobre todos los promotores
#     high-conf de cada replica, como ruido_cgi / ruido_tata); despues se
#     restringe a los promotores FP sin N (los que tienen busqueda del motivo).
#   - Actividad alta = percentil 75 de mean dentro de cada replica.
#   - TATA-box: R/functions/tata_scan.R (matriz de Bucher con pesos W de EPD,
#     corte -8.17; clases canonica / solo no canonica / solo hebra opuesta /
#     sin motivo). "Sin motivo" = sin ninguna coincidencia en el fragmento.
#   - Puntaje = mejor puntaje en -31 a -26 (misma hebra). Puntaje alto / bajo
#     de las canonicas = por encima / por debajo de la mediana de las
#     canonicas (ambas replicas juntas, -5.67). "TATA fuertes" = canonicas de
#     puntaje alto.
#   - AUC = P(ruido de los casos > ruido de la referencia), exacta (pROC), con
#     IC 95% de DeLong; no se estima con menos de 10 casos.
#
# Paneles (cada uno en .jpg/.pdf y con un .tsv de los valores graficados):
#   A panel_a_modelo_AT: coeficientes de clase vs. sin motivo sobre el rango
#     de varianza (lm por replica y estrato de actividad), sin covariable y
#     ajustado por A/T proximal (-50 a -1) y del fragmento (1 - g_c).
#   B panel_b_puntaje_actividad: AUC de ruido de las canonicas (todas, puntaje
#     bajo, puntaje alto) vs. sin motivo del mismo estrato de actividad.
#   C panel_c_gini: indice de Gini de la expresion en FANTOM5
#     (sample_specificity_gini; NA en los no detectados) de las TATA fuertes
#     y del resto, por estrato de actividad. Violin con mediana (linea) y
#     mediana + rango intercuartil (punto y barra).
#   D panel_d_ruido_gini: AUC de ruido de las TATA fuertes vs. sin motivo del
#     mismo tercil de Gini (terciles sobre todos los promotores con Gini de
#     cada replica) y del mismo estrato de actividad. Los no detectados en
#     FANTOM5 quedan afuera. En B y D, las celdas con menos de 20 casos
#     muestran solo el n (la AUC queda en el .tsv).
# Replicas siempre en facetas (columnas). B y D comparten el eje de AUC.
# Figura completa: tata_ruido_especificidad.jpg/.pdf.
#
# Requiere: data/raw/library.bed, BSgenome.Hsapiens.UCSC.hg38,
# data/processed/activity_stats_highconf.tsv y data/processed/prom_df.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
library(pROC)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
source("R/functions/plot_helpers.R")
source("R/functions/tata_scan.R")
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter

out_dir <- fig_dir("tata_ruido_especificidad")

CANON <- TATA_CANON
MIN_N_AUC <- 10
MIN_N_GRAFICO <- 20 # celdas con menos casos: solo la etiqueta de n en B y D (el .tsv trae la AUC)
col_alta <- "#216869"
col_baja <- "#14AFB2"
th <- theme_pubr(base_size = 16) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")
tag_theme <- theme(plot.tag = element_text(size = 24, face = "bold"))
num <- function(x) format(x, decimal.mark = ".", big.mark = "", trim = TRUE)
save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}
save_tsv <- function(x, name) write_tsv(x, file.path(out_dir, name))

auc_vs <- function(x, y) {
  if (length(x) < MIN_N_AUC || length(y) < MIN_N_AUC) {
    return(tibble(n = length(x), n_ref = length(y), AUC = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p = NA_real_))
  }
  r <- pROC::roc(controls = y, cases = x, direction = "<", ci = TRUE, quiet = TRUE)
  auc <- as.numeric(r$auc)
  se <- sqrt(pROC::var(r, method = "delong"))
  tibble(n = length(x), n_ref = length(y), AUC = auc, conf.low = r$ci[1], conf.high = r$ci[3], p = 2 * pnorm(-abs((auc - 0.5) / se)))
}

# --- Datos ----------------------------------------------------------------------

lib <- tata_load_library()
code <- lib$code[!lib$has_n, ] # 1 fragmento con N fuera
sc <- tata_scores(code)
hits <- tata_hits(sc$misma, sc$opuesta, TATA_CORTE_EPD)
at_prox <- tata_at_prox(code, -50, -1)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  filter(seq_id %in% rownames(code)) %>%
  group_by(rep) %>%
  mutate(act_pct = 100 * percent_rank(mean), actividad = factor(if_else(act_pct >= 75, "Alta", "Baja"), levels = c("Alta", "Baja"))) %>%
  ungroup() %>%
  inner_join(tata_classify(hits, rownames(code)), by = "seq_id") %>%
  mutate(
    score_canon = tata_best_in_window(sc$misma, CANON[1], CANON[2])[seq_id],
    at_prox = at_prox[seq_id], at_frag = 1 - g_c
  )
stopifnot(all((data$clase == "Canónica") == data$TATA_EPD))

corte_puntaje <- round(median(data$score_canon[data$clase == "Canónica"]), 2)
data <- data %>%
  mutate(tata_fuerte = clase == "Canónica" & score_canon >= corte_puntaje - TATA_TOL)
message("Corte de puntaje (mediana de las canonicas): ", corte_puntaje)

# --- A: coeficientes de clase, con y sin ajuste por A/T ----------------------------

modelos_lab <- c("Sin ajuste", "Ajustado por A/T")
clases_lab <- setdiff(TATA_CLASES, "Sin motivo")
panel_a_df <- data %>%
  mutate(clase = relevel(clase, ref = "Sin motivo")) %>%
  group_by(actividad, rep) %>%
  group_modify(function(d, k) {
    bind_rows(
      broom::tidy(lm(var_rank_sw ~ clase, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[1]),
      broom::tidy(lm(var_rank_sw ~ clase + at_prox + at_frag, data = d), conf.int = TRUE) %>% mutate(modelo = modelos_lab[2])
    ) %>%
      filter(str_detect(term, "^clase")) %>%
      mutate(clase = str_remove(term, "^clase"), n = map_int(clase, ~ sum(d$clase == .x)))
  }) %>%
  ungroup() %>%
  select(actividad, rep, clase, n, modelo, estimate, conf.low, conf.high, p.value)
save_tsv(panel_a_df, "panel_a_modelo_AT.tsv")

panel_a <- panel_a_df %>%
  mutate(
    clase = factor(clase, levels = rev(clases_lab)),
    modelo = factor(modelo, levels = rev(modelos_lab)),
    actividad = factor(paste("Actividad", tolower(actividad)), levels = c("Actividad alta", "Actividad baja"))
  ) %>%
  ggplot(aes(estimate, clase, color = modelo)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.6), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.6), size = 2.8) +
  geom_text(
    data = . %>% distinct(actividad, rep, clase, n), aes(x = Inf, y = clase, label = paste0("n=", num(n))),
    inherit.aes = FALSE, hjust = 1.05, size = 4, color = "grey30"
  ) +
  scale_color_manual(values = c("Sin ajuste" = "grey55", "Ajustado por A/T" = col_alta), breaks = modelos_lab) +
  scale_x_continuous(labels = num, expand = expansion(mult = c(0.05, 0.25))) +
  facet_grid(actividad ~ rep) +
  th +
  labs(x = "Diferencia de rango de varianza vs. sin motivo", y = NULL, color = NULL)
save_jpg_pdf(panel_a, "panel_a_modelo_AT", 10, 7)

# --- B: puntaje x actividad ---------------------------------------------------------

panel_b_df <- data %>%
  group_by(rep, actividad) %>%
  group_modify(function(d, k) {
    ref <- d$var_rank_sw[d$clase == "Sin motivo"]
    cn <- d[d$clase == "Canónica", ]
    alto <- cn$score_canon >= corte_puntaje - TATA_TOL
    bind_rows(
      auc_vs(cn$var_rank_sw, ref) %>% mutate(puntaje = "Todas"),
      auc_vs(cn$var_rank_sw[!alto], ref) %>% mutate(puntaje = "Bajo"),
      auc_vs(cn$var_rank_sw[alto], ref) %>% mutate(puntaje = "Alto")
    )
  }) %>%
  ungroup() %>%
  mutate(corte_puntaje = corte_puntaje) %>%
  select(rep, actividad, puntaje, everything())
save_tsv(panel_b_df, "panel_b_puntaje_actividad.tsv")

# --- D: ruido de las TATA fuertes por tercil de Gini -------------------------------

espec_lab <- c("1" = "Baja", "2" = "Media", "3" = "Alta")
panel_d_df <- data %>%
  filter(!is.na(sample_specificity_gini)) %>%
  group_by(rep) %>%
  mutate(tercil_gini = factor(espec_lab[as.character(ntile(sample_specificity_gini, 3))], levels = espec_lab)) %>%
  group_by(rep, actividad, tercil_gini) %>%
  group_modify(~ auc_vs(.x$var_rank_sw[.x$tata_fuerte], .x$var_rank_sw[.x$clase == "Sin motivo"]) %>%
    mutate(gini_min = min(.x$sample_specificity_gini), gini_max = max(.x$sample_specificity_gini))) %>%
  ungroup()
save_tsv(panel_d_df, "panel_d_ruido_gini.tsv")

# Eje de AUC comun a B y D, con todos los IC completos
ic_graficados <- bind_rows(panel_b_df, panel_d_df) %>% filter(n >= MIN_N_GRAFICO)
lims_auc <- range(c(0.5, ic_graficados$conf.low, ic_graficados$conf.high), na.rm = TRUE)
lims_auc <- lims_auc + c(-0.02, 0.06) # lugar para las etiquetas de n

forest_auc <- function(df, x, color, valores, titulo_color, xlab) {
  dodge <- position_dodge(width = 0.75)
  df %>%
    mutate(across(c(AUC, conf.low, conf.high), ~ if_else(n < MIN_N_GRAFICO, NA_real_, .x)), y_lab = coalesce(conf.high, 0.5)) %>%
    ggplot(aes(.data[[x]], AUC, color = .data[[color]], group = .data[[color]])) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, position = dodge, linewidth = 1, na.rm = TRUE) +
    geom_point(size = 3, position = dodge, na.rm = TRUE) +
    geom_text(aes(y = y_lab, label = paste0("n=", num(n))), vjust = -0.6, size = 3.6, position = dodge, show.legend = FALSE) +
    scale_color_manual(values = valores) +
    scale_y_continuous(labels = num, limits = lims_auc, expand = expansion(mult = 0)) +
    facet_wrap(~rep) +
    th +
    labs(x = xlab, y = "AUC de ruido vs.\npromotores sin motivo", color = titulo_color)
}

panel_b <- panel_b_df %>%
  mutate(puntaje = factor(puntaje, levels = c("Todas", "Bajo", "Alto")), actividad = factor(actividad, levels = c("Baja", "Alta"))) %>%
  forest_auc("actividad", "puntaje", c(Todas = "grey55", Bajo = col_baja, Alto = col_alta), "Puntaje de la TATA canónica", "Actividad")
save_jpg_pdf(panel_b, "panel_b_puntaje_actividad", 10, 6.5)

panel_d <- panel_d_df %>%
  mutate(actividad = factor(actividad, levels = c("Baja", "Alta"))) %>%
  forest_auc("tercil_gini", "actividad", c(Baja = col_baja, Alta = col_alta), "Actividad", "Especificidad tisular (tercil del índice de Gini)")
save_jpg_pdf(panel_d, "panel_d_ruido_gini", 10, 6.5)

# --- C: indice de Gini por grupo -------------------------------------------------

panel_c_data <- data %>%
  filter(!is.na(sample_specificity_gini)) %>%
  mutate(
    grupo = factor(if_else(tata_fuerte, "TATA fuertes", "Resto"), levels = c("TATA fuertes", "Resto")),
    actividad = factor(actividad, levels = c("Baja", "Alta"))
  )
panel_c_df <- panel_c_data %>%
  group_by(rep, grupo, actividad) %>%
  summarise(
    n = n(), mediana = median(sample_specificity_gini),
    q1 = quantile(sample_specificity_gini, 0.25), q3 = quantile(sample_specificity_gini, 0.75), .groups = "drop"
  ) %>%
  left_join(
    panel_c_data %>%
      group_by(rep, grupo) %>%
      summarise(
        p_wilcoxon_alta_vs_baja = wilcox.test(sample_specificity_gini[actividad == "Alta"], sample_specificity_gini[actividad == "Baja"], exact = FALSE)$p.value,
        .groups = "drop"
      ),
    by = c("rep", "grupo")
  )
save_tsv(panel_c_df, "panel_c_gini.tsv")

dodge_c <- position_dodge(width = 0.85)
panel_c <- ggplot(panel_c_data, aes(actividad, sample_specificity_gini, fill = grupo)) +
  geom_violin(position = dodge_c, alpha = 0.7, color = NA) +
  # convencion de la tesis: linea horizontal = mediana, barra vertical = IQR
  geom_errorbar(data = panel_c_df, aes(y = mediana, ymin = mediana, ymax = mediana, group = grupo), position = dodge_c, width = 0.35, linewidth = 0.8, color = "black", inherit.aes = TRUE, show.legend = FALSE) +
  geom_pointrange(data = panel_c_df, aes(y = mediana, ymin = q1, ymax = q3, group = grupo), position = dodge_c, size = 0.3, color = "black", show.legend = FALSE) +
  geom_text(data = panel_c_df, aes(y = 1.04, label = paste0("n=", num(n)), group = grupo), position = dodge_c, size = 4, color = "grey30") +
  scale_fill_manual(values = c("TATA fuertes" = col_alta, "Resto" = "grey70")) +
  scale_y_continuous(labels = num, limits = c(0.1, 1.08), breaks = seq(0.2, 1, 0.2)) +
  facet_wrap(~rep) +
  th +
  labs(x = "Actividad", y = "Especificidad tisular (índice de Gini)", fill = NULL)
save_jpg_pdf(panel_c, "panel_c_gini", 10, 6.5)

# --- Figura completa --------------------------------------------------------------

figura <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "A") & tag_theme
save_jpg_pdf(figura, "tata_ruido_especificidad", 20, 14)

message("Figura guardada en ", out_dir)
