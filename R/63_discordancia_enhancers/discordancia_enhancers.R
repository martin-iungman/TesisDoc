# discordancia_enhancers (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Discordancia con signo entre el reportero y la actividad endogena maxima
# (FANTOM5) segun el numero de enhancers cercanos, por clase de promotor y
# replica. Pregunta: el efecto de los enhancers en discordancia_features se
# sostiene dentro de housekeeping y de tejido-especificos, o es un proxy de
# especificidad tisular?
#
# Conjunto y discordancia: los de discordancia_features (hek_tpm > 0 y luego
# add_rank_discordance(); dif_signed = rr - re, positivo = el reportero
# sobreestima). Se exige que la binaria reproduzca discordancia_features_lmm.tsv.
#
# Clases: terciles de sample_specificity_gini calculados sobre los detectados y
# aplicados tal cual al subconjunto activo en HEK293 (T1 = housekeeping,
# T3 = tejido-especifico), como en especificidad_tisular_fantom.
#
# Enhancers: enh50kb de prom_df, regiones de FANTOM5
# (data/external/F5.hg38.enhancers.bed) que solapan una ventana de 50 kb de ancho
# total (+-25 kb) anclada en el extremo 5' del oligo.
#
# Grupos: cortes fijos, iguales en las dos replicas. Los puntos son medias crudas
# con IC 95% t. El numero sobre cada eje es el n de ese punto, por clase.
#
# Anotacion en cada panel: coeficiente por enhancer adicional de
# dif_signed ~ avg_rank + enh50kb, por clase y replica (IC 95% y p).
#
# La tabla agrega, por clase y replica: el contraste de cada grupo contra "sin
# enhancer" ajustado por avg_rank, la pendiente solo con promotores con enhancer
# y la binaria. La linea en 0 marca igualdad de rangos, pero el IC de un grupo
# contra 0 no es el test del efecto: lo es la pendiente o el contraste.
#
# Salidas: figura en fig_dir(slug) y tables/discordancia_enhancers.tsv.
# No hay remuestreo, por eso no hay set.seed.
#
# Requiere: data/processed/activity_stats_highconf.tsv, data/processed/prom_df.tsv,
# data/processed/fantom_endo_activity_summary.tsv y
# tables/discordancia_features_lmm.tsv. Run from the TesisDoc repo root.

library(tidyverse)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "discordancia_enhancers"
out_dir <- fig_dir(slug)

clase_colors <- c("Housekeeping" = "#1B8C8E", "Tejido-específico" = "#0D2C54")
grupo_niveles <- c("Sin enhancer", "1", "2", "3", "4-5", "≥6")

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

base <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id")

terciles <- base %>%
  filter(!is.na(sample_specificity_gini)) %>%
  group_by(rep) %>%
  mutate(tercil = as.numeric(cut_number(sample_specificity_gini, n = 3))) %>%
  ungroup() %>%
  select(seq_id, rep, tercil)

data <- base %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  add_rank_discordance() %>%
  left_join(terciles, by = c("seq_id", "rep"))

referencia <- read_tsv("tables/discordancia_features_lmm.tsv", show_col_types = FALSE) %>%
  filter(caracteristica == "Enhancers a 50kb")
fit_ref <- lm(dif_signed ~ avg_rank + rep + I(enh50kb > 0), data = data)
stopifnot(
  nrow(data) == referencia$n_mediciones,
  abs(coef(fit_ref)[["I(enh50kb > 0)TRUE"]] - referencia$estimate) < 1e-3
)

datos_clase <- data %>%
  filter(tercil %in% c(1, 3)) %>%
  mutate(
    clase = if_else(tercil == 1, "Housekeeping", "Tejido-específico"),
    grupo = case_when(
      enh50kb == 0 ~ "Sin enhancer",
      enh50kb %in% 1:3 ~ as.character(enh50kb),
      enh50kb %in% 4:5 ~ "4-5",
      TRUE ~ "≥6"
    ) %>% factor(grupo_niveles)
  )

# --- Modelos -------------------------------------------------------------------

coef_de <- function(fit, termino) {
  ci <- confint(fit)
  s <- summary(fit)$coefficients
  tibble(
    termino = termino, estimado = s[termino, "Estimate"], ic_inf = ci[termino, 1],
    ic_sup = ci[termino, 2], p = s[termino, "Pr(>|t|)"]
  )
}

modelos <- datos_clase %>%
  group_by(rep, clase) %>%
  group_modify(function(d, k) {
    contraste <- lm(dif_signed ~ avg_rank + grupo, data = d)
    terminos <- paste0("grupo", grupo_niveles[-1])
    bind_rows(
      coef_de(lm(dif_signed ~ avg_rank + enh50kb, data = d), "enh50kb") %>%
        mutate(tipo = "pendiente por enhancer", n = nrow(d)),
      coef_de(lm(dif_signed ~ avg_rank + enh50kb, data = filter(d, enh50kb > 0)), "enh50kb") %>%
        mutate(tipo = "pendiente solo con enhancer", n = sum(d$enh50kb > 0)),
      coef_de(lm(dif_signed ~ avg_rank + I(enh50kb > 0), data = d), "I(enh50kb > 0)TRUE") %>%
        mutate(tipo = "con vs. sin enhancer", n = nrow(d)),
      map(terminos, ~ coef_de(contraste, .x)) %>%
        list_rbind() %>%
        mutate(tipo = paste0("contraste contra sin enhancer: ", sub("grupo", "", termino)), n = NA_integer_)
    )
  }) %>%
  ungroup() %>%
  select(-termino)

medias <- datos_clase %>%
  group_by(rep, clase, grupo) %>%
  summarise(
    n = n(), media = mean(dif_signed), se = sd(dif_signed) / sqrt(n()),
    ic_inf = media - qt(0.975, n() - 1) * se, ic_sup = media + qt(0.975, n() - 1) * se,
    .groups = "drop"
  )

write_tsv(
  bind_rows(
    medias %>%
      transmute(rep, clase, tipo = paste0("media cruda: ", grupo), n, estimado = media, ic_inf, ic_sup, p = NA_real_),
    modelos
  ),
  "tables/discordancia_enhancers.tsv"
)

# --- Figura --------------------------------------------------------------------

rango_y <- range(c(medias$ic_inf, medias$ic_sup))
paso <- diff(rango_y) * 0.07

formato_p <- function(p) {
  if (p >= 0.001) sprintf("%.3f", p) else sub("e-0*", "e-", sprintf("%.1e", p))
}

rotulo_beta <- modelos %>%
  filter(tipo == "pendiente por enhancer") %>%
  mutate(
    rotulo = sprintf("β = %.4f [%.4f, %.4f]; p = %s", estimado, ic_inf, ic_sup, map_chr(p, formato_p)),
    y = rango_y[2] + paso * if_else(clase == "Housekeeping", 2, 1)
  )

filas_n <- medias %>%
  mutate(y = rango_y[1] - paso * if_else(clase == "Housekeeping", 1, 2))

p <- ggplot(medias, aes(grupo, media, color = clase, group = clase)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(data = ~ filter(.x, grupo != "Sin enhancer"), linewidth = 0.6) +
  geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0.2, alpha = 0.8) +
  geom_point(size = 2) +
  geom_text(data = filas_n, aes(y = y, label = n), size = 2.8, show.legend = FALSE) +
  geom_text(
    data = rotulo_beta, aes(x = 6.45, y = y, label = rotulo), inherit.aes = FALSE,
    color = clase_colors[rotulo_beta$clase], hjust = 1, size = 3, show.legend = FALSE
  ) +
  facet_wrap(~rep) +
  scale_x_discrete(limits = grupo_niveles) +
  scale_color_manual(values = clase_colors) +
  labs(
    x = "Enhancers en ±25 kb",
    y = "Discordancia media (reportero - endógeno), IC 95%", color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, paste0(slug, ".jpg")), p, width = 10, height = 5.5, units = "in", dpi = 300)
ggsave(file.path(out_dir, paste0(slug, ".pdf")), p, width = 10, height = 5.5, units = "in", device = cairo_pdf)

message("Figura en ", out_dir, "; tabla en tables/discordancia_enhancers.tsv")
