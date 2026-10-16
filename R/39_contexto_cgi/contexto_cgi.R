# EXPLORATORIO - efecto de las caracteristicas de contexto endogeno
# (grupo "endo" de feature_groups, las mismas de summary_features_endogeno)
# sobre la actividad media, por separado dentro de los promotores con y sin
# isla CpG (CGI_frag: G+C >= 50% y CpG o/e > 0.6 en el fragmento de 252 pb).
#
# Mismos criterios que summary_features_endogeno
# (R/15_summary_features/summary_features.R):
#   - Hodges-Lehmann + IC 95% de mean ~ feature, cada replica por separado
#     (coin::wilcox_test(conf.int = TRUE), igual que wilcox_effect_summary()).
#   - p-valor de Wilcoxon estratificado por replica (mean ~ feature | rep).
#   - BH dentro de cada grupo (Global / Con isla CpG / Sin isla CpG; familia =
#     las caracteristicas de contexto de ese grupo). En la figura de tesis la
#     familia son las 34 features (seq + endo), por eso p_corr del grupo Global
#     difiere un poco de activity_summary.tsv (los HL e IC son identicos).
#   - "Sobrevive" = p_corr < 0.05 y estimado e IC (limites inferior y
#     superior) del mismo signo en ambas replicas.
# Las caracteristicas definidas por terciles/medianas (conservacion,
# especificidad tisular, shape, HEK293, DNase-seq) usan los cortes globales de
# build_tidy_features() (por replica, sobre todos los promotores): no se
# recalculan dentro de cada grupo de isla CpG.
#
# Salidas (figures/exploratorio_contexto_cgi/, sin fig_dir: no es figura de
# la tesis):
#   contexto_cgi_efectos.tsv  formato largo: grupo x caracteristica x replica
#                             (n, prevalencia, HL, IC, p, p_corr, sobrevive)
#   contexto_cgi_resumen.tsv  una fila por caracteristica: prevalencias
#                             (Rep 1 / Rep 2, sobre los no-NA),
#                             efectos global / con / sin isla CpG y clase
#                             (ambos grupos / uno solo / ninguno)
#   contexto_cgi_forest.jpg   forest plot de los tres efectos por caracteristica
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(coin)
source("R/functions/plot_helpers.R")

out_dir <- "figures/exploratorio_contexto_cgi"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

GRUPOS <- c("Global", "Con isla CpG", "Sin isla CpG")

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

tidy_data <- build_tidy_features(data, keep = "mean")
endo_features <- feature_groups$feature[feature_groups$group == "endo"]

datos_grupo <- list(
  "Global" = tidy_data,
  "Con isla CpG" = tidy_data %>% filter(`Islas CpG` == "TRUE"),
  "Sin isla CpG" = tidy_data %>% filter(`Islas CpG` == "FALSE")
)

# --- Efecto de una caracteristica dentro de un grupo -------------------------
# Devuelve NA (en vez de error) si en el grupo la caracteristica tiene un solo
# nivel en alguna replica.

efecto_feature <- function(df, feature) {
  d <- df %>%
    select(mean, rep, x = all_of(feature)) %>%
    filter(!is.na(x))
  pval <- tryCatch(pvalue(coin::wilcox_test(mean ~ x | rep, data = d)) %>% as.numeric(), error = function(e) NA_real_)
  map(c("Rep 1", "Rep 2"), function(rep_id) {
    dr <- d %>% filter(rep == rep_id)
    ci <- tryCatch(confint(coin::wilcox_test(mean ~ x, data = dr, conf.int = TRUE)), error = function(e) NULL)
    tibble(
      feature = feature, rep = rep_id,
      n = nrow(dr), n_con = sum(dr$x == "TRUE"),
      prevalencia = mean(dr$x == "TRUE"),
      estimate = if (is.null(ci)) NA_real_ else unname(ci$estimate),
      P2.5 = if (is.null(ci)) NA_real_ else ci$conf.int[1],
      P97.5 = if (is.null(ci)) NA_real_ else ci$conf.int[2],
      pval = pval
    )
  }) %>% list_rbind()
}

efectos <- imap(datos_grupo, function(df, grupo) {
  map(endo_features, ~ efecto_feature(df, .x)) %>%
    list_rbind() %>%
    mutate(grupo = grupo)
}) %>%
  list_rbind() %>%
  mutate(grupo = factor(grupo, levels = GRUPOS)) %>%
  arrange(feature, grupo, rep)

# BH dentro de cada grupo, sobre un p-valor por caracteristica (en el formato
# largo pval se repite una vez por replica).
pcorr <- efectos %>%
  distinct(grupo, feature, pval) %>%
  group_by(grupo) %>%
  mutate(pval_corr = p.adjust(pval, "BH")) %>%
  ungroup()
efectos <- efectos %>%
  left_join(pcorr, by = c("grupo", "feature", "pval")) %>%
  group_by(grupo, feature) %>%
  mutate(
    sobrevive = !is.na(pval_corr[1]) & pval_corr[1] < 0.05 &
      all(!is.na(c(estimate, P2.5, P97.5))) &
      n_distinct(sign(c(estimate, P2.5, P97.5))) == 1 & all(c(estimate, P2.5, P97.5) != 0)
  ) %>%
  ungroup() %>%
  select(grupo, feature, rep, n, n_con, prevalencia, estimate, P2.5, P97.5, pval, pval_corr, sobrevive)

write_tsv(efectos, file.path(out_dir, "contexto_cgi_efectos.tsv"))

# --- Resumen: una fila por caracteristica ------------------------------------

fmt_efecto <- function(est, lo, hi) sprintf("%.3f [%.3f, %.3f]", est, lo, hi)

resumen_grupo <- efectos %>%
  group_by(grupo, feature) %>%
  summarise(
    prevalencia = sprintf("%.1f%% / %.1f%%", 100 * prevalencia[rep == "Rep 1"], 100 * prevalencia[rep == "Rep 2"]),
    efecto = paste0(
      "Rep 1 ", fmt_efecto(estimate[rep == "Rep 1"], P2.5[rep == "Rep 1"], P97.5[rep == "Rep 1"]),
      "; Rep 2 ", fmt_efecto(estimate[rep == "Rep 2"], P2.5[rep == "Rep 2"], P97.5[rep == "Rep 2"])
    ),
    signo = case_when(all(estimate > 0) ~ "+", all(estimate < 0) ~ "-", TRUE ~ "mixto"),
    pval_corr = signif(pval_corr[1], 2),
    sobrevive = sobrevive[1],
    .groups = "drop"
  )

resumen <- resumen_grupo %>%
  mutate(grupo = recode(grupo, "Global" = "global", "Con isla CpG" = "con_cgi", "Sin isla CpG" = "sin_cgi")) %>%
  pivot_wider(names_from = grupo, values_from = c(prevalencia, efecto, signo, pval_corr, sobrevive), names_vary = "slowest") %>%
  mutate(
    sobrevive_en = case_when(
      sobrevive_con_cgi & sobrevive_sin_cgi ~ "ambos grupos",
      sobrevive_con_cgi ~ "solo con isla CpG",
      sobrevive_sin_cgi ~ "solo sin isla CpG",
      TRUE ~ "ninguno"
    ),
    # signo opuesto al global dentro de un grupo donde sobrevive
    inversion = (sobrevive_con_cgi & signo_con_cgi != signo_global) |
      (sobrevive_sin_cgi & signo_sin_cgi != signo_global)
  ) %>%
  left_join(
    efectos %>% filter(grupo == "Global", rep == "Rep 1") %>% select(feature, est_global_rep1 = estimate),
    by = "feature"
  ) %>%
  arrange(desc(est_global_rep1)) %>%
  select(feature, prevalencia_global, prevalencia_con_cgi, prevalencia_sin_cgi, sobrevive_global, sobrevive_en, inversion, everything(), -est_global_rep1)

write_tsv(resumen, file.path(out_dir, "contexto_cgi_resumen.tsv"))

# --- Forest plot --------------------------------------------------------------
# Un punto por grupo (global / con / sin isla CpG) y replica; relleno =
# sobrevive (p_corr < 0.05 y signo consistente), hueco = no.

orden <- efectos %>%
  filter(grupo == "Global", rep == "Rep 1") %>%
  arrange(estimate) %>%
  pull(feature)

# Desplazamiento vertical explicito por grupo (position_dodge no desplaza
# igual a geom_errorbarh y geom_point en un eje y discreto).
desplazamiento <- c("Global" = 0.25, "Con isla CpG" = 0, "Sin isla CpG" = -0.25)

forest <- efectos %>%
  filter(!is.na(estimate)) %>%
  mutate(
    y = match(feature, orden) + desplazamiento[as.character(grupo)],
    relleno = ifelse(sobrevive, as.character(grupo), "no")
  ) %>%
  ggplot(aes(x = estimate, y = y, colour = grupo)) +
  geom_vline(xintercept = 0, colour = "#999999", linewidth = 0.5) +
  geom_segment(aes(x = P2.5, xend = P97.5, yend = y), linewidth = 0.7) +
  geom_point(aes(fill = relleno), shape = 21, size = 2.6, stroke = 1) +
  scale_y_continuous(breaks = seq_along(orden), labels = orden, expand = expansion(add = 0.5)) +
  facet_wrap(~rep, nrow = 1) +
  scale_colour_manual(
    values = c("Global" = "#5f5e5a", "Con isla CpG" = "#2a78d6", "Sin isla CpG" = "#eb6834"),
    breaks = GRUPOS
  ) +
  scale_fill_manual(
    values = c("Global" = "#5f5e5a", "Con isla CpG" = "#2a78d6", "Sin isla CpG" = "#eb6834", "no" = "white"),
    guide = "none"
  ) +
  theme_pubr(base_size = 14) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "#eeeeee"), strip.background = element_blank()) +
  labs(
    x = "Efecto sobre la actividad (Hodges-Lehmann, IC 95%)", y = NULL, colour = NULL,
    caption = "Punto relleno: p BH < 0,05 y mismo signo (estimado e IC) en ambas réplicas; hueco: no."
  )

ggsave(file.path(out_dir, "contexto_cgi_forest.jpg"), forest, width = 13, height = 9, units = "in", dpi = 300)

print(resumen %>% select(feature, prevalencia_global, prevalencia_con_cgi, prevalencia_sin_cgi, sobrevive_global, sobrevive_en, inversion), n = Inf, width = Inf)
message("Salidas guardadas en ", out_dir)
