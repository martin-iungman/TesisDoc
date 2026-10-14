# EXPLORATORIO - edad del elemento transponible (TE) vs. actividad media.
# Pone a prueba la interpretacion de R2.8 (summary_seq: elementos
# transponibles, SINE y LINE asociados a actividad baja) de "origen
# reciente / menor presion selectiva": si fuera cierta, los TE mas
# jovenes (menor divergencia de RepeatMasker respecto de su consenso)
# deberian tener menos actividad que los viejos, o sea correlacion
# positiva entre perc_div y actividad.
#
# Edad = perc_div (milliDiv/10) del TE con mayor superposicion con el
# fragmento (data/processed/te_age.tsv, ver
# R/00_prom_features/build_te_age.R). Ojo: perc_div es divergencia cruda,
# no corregida por CpG (los CpG de TE jovenes mutan rapido y la sesgan
# hacia arriba en elementos ricos en CpG).
#
# Salidas en figures/exploratorio_te_edad/ (no es una figura numerada de
# la tesis todavia, por eso no pasa por fig_dir()).
# Requiere activity_stats_highconf.tsv, prom_df.tsv y te_age.tsv. Run from
# the TesisDoc repo root.

library(tidyverse)
library(ggpubr)

out_dir <- "figures/exploratorio_te_edad"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
te_age <- read_tsv("data/processed/te_age.tsv", show_col_types = FALSE) %>%
  left_join(read_tsv("data/processed/te_age_kimura.tsv", show_col_types = FALSE), by = "seq_id")

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(te_age %>% select(-TE_superclass), by = "seq_id") %>%
  mutate(
    es_te = !is.na(TE_perc_div_dom),
    sc = str_remove(TE_class, "/.+$"),
    TE_frac = TE_overlap_dom / 252
  )

te <- data %>% filter(es_te, sc %in% c("SINE", "LINE", "LTR", "DNA"))
cat("Promotores con TE (SINE/LINE/LTR/DNA) por replica:
")
print(count(te, rep, sc) %>% pivot_wider(names_from = sc, values_from = n))
cat("
Promotores sin TE por replica:
")
print(data %>% filter(!es_te) %>% count(rep))
cat("
% de promotores con TE con mas de 1 TE superpuesto (ahi dominante != ponderado):",
  round(100 * mean(te_age$TE_n > 1), 1), "
")

# Dos medidas de edad: TE dominante (mayor superposicion) y promedio de
# perc_div ponderado por pb superpuestos con el fragmento.
medidas <- c(
  "Cruda, TE dominante" = "TE_perc_div_dom",
  "Cruda, promedio ponderado" = "TE_perc_div_wmean",
  "Kimura sin CpG, dominante" = "kimura_nocpg_dom",
  "Kimura con CpG, dominante" = "kimura_cpg_dom",
  "Kimura con CpG, ponderado" = "kimura_cpg_wmean",
  "Mediana de la subfamilia" = "TE_subfam_med_div"
)

# Spearman con IC95% (Bonett-Wright: SE = sqrt((1 + rho^2/2) / (n - 3)) sobre z de Fisher)
spearman_ci <- function(x, y) {
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  r <- cor(x, y, method = "spearman")
  n <- length(x)
  z <- atanh(r); se <- sqrt((1 + r^2 / 2) / (n - 3))
  tibble(n = n, rho = r, lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se),
    p = cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
}

cor_tab <- imap(medidas, function(col, nom) {
  bind_rows(
    te %>% group_by(rep) %>% reframe(grupo = "Todos", spearman_ci(.data[[col]], mean)),
    te %>% group_by(rep, grupo = sc) %>% reframe(spearman_ci(.data[[col]], mean))
  ) %>% mutate(medida = nom)
}) %>% list_rbind() %>% arrange(medida, grupo, rep)
cat("
Spearman edad vs actividad media (rho>0: mas viejo = mas actividad):
")
print(cor_tab %>% mutate(across(c(rho, lo, hi, p), ~ signif(.x, 3))), n = 40)

mod_tab <- imap(medidas, function(col, nom) {
  map(list("cobertura + G+C + superclase" = c("TE_frac", "g_c", "sc"), "... + CpG o/e del fragmento" = c("TE_frac", "g_c", "sc", "cpg_oe")), function(cov) {
    te %>% group_by(rep) %>% group_modify(~ {
      m <- lm(reformulate(c(col, cov), "mean"), data = .x)
      broom::tidy(m, conf.int = TRUE) %>% filter(term == col)
    }) %>% ungroup()
  }) %>% list_rbind(names_to = "ajuste") %>% mutate(medida = nom)
}) %>% list_rbind()
cat("
lm(mean ~ edad + covariables), coef. de la edad (por 1 punto de % div):
")
print(mod_tab %>% select(medida, ajuste, rep, estimate, conf.low, conf.high, p.value) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 30)

# --- Cuartiles de edad vs. promotores sin TE -----------------------------------
base <- data %>% filter(!es_te) %>% group_by(rep) %>% summarise(med = median(mean), .groups = "drop")
cuart <- imap(medidas[c("Cruda, TE dominante", "Kimura con CpG, dominante")], function(col, nom) {
  te %>% group_by(rep) %>%
    mutate(edad = cut(.data[[col]], breaks = quantile(.data[[col]], 0:4 / 4, na.rm = TRUE), include.lowest = TRUE, labels = c("Q1 (mas joven)", "Q2", "Q3", "Q4 (mas viejo)"))) %>%
    ungroup() %>% mutate(medida = nom)
}) %>% list_rbind()
cat("
Mediana de actividad por cuartil de edad (y sin TE):
")
print(cuart %>% filter(!is.na(edad)) %>% group_by(medida, rep, edad) %>% summarise(n = n(), mediana = round(median(mean), 2), .groups = "drop"), n = 20)
print(base)

# --- Grafico resumen: forest de Spearman -----------------------------------------
orden <- c("Todos", "SINE", "LINE", "LTR", "DNA")
ns <- cor_tab %>% filter(medida == "Cruda, TE dominante") %>% group_by(grupo) %>% summarise(lab = paste0(grupo[1], "
(n = ", paste(n, collapse = " / "), ")"), .groups = "drop")
p_f <- cor_tab %>%
  left_join(ns, by = "grupo") %>%
  mutate(grupo = factor(lab, levels = rev(ns$lab[order(match(ns$grupo, orden))]))) %>%
  ggplot(aes(rho, grupo, col = rep)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey40") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.7) +
  facet_wrap(~medida, nrow = 2) +
  scale_color_manual(values = c("Rep 1" = "#216869", "Rep 2" = "#D6741F")) +
  theme_pubclean(base_size = 18) +
  labs(x = "Spearman: divergencia del TE vs. actividad media (IC 95%)", y = NULL, col = NULL,
    caption = "rho > 0: TE mas viejos con mas actividad (lo esperado si los TE jovenes reprimen)")
ggsave(file.path(out_dir, "te_edad_resumen.jpg"), p_f, width = 15, height = 9, units = "in")

p_q <- ggplot(cuart %>% filter(!is.na(edad)), aes(edad, mean)) +
  geom_hline(data = base, aes(yintercept = med), linetype = "dashed", col = "grey40") +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "#14AFB2", alpha = 0.6) +
  facet_grid(medida ~ rep) +
  coord_cartesian(ylim = c(1, 3.5)) +
  theme_pubclean(base_size = 18) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Divergencia del TE (cuartiles)", y = "Actividad media", caption = "Linea punteada: mediana de promotores sin TE")
ggsave(file.path(out_dir, "te_edad_cuartiles.jpg"), p_q, width = 12, height = 8, units = "in")

write_tsv(cor_tab, file.path(out_dir, "te_edad_spearman.tsv"))
write_tsv(mod_tab, file.path(out_dir, "te_edad_lm.tsv"))
