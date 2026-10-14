# EXPLORATORIO - los promotores con TE, son evolutivamente mas "nuevos"
# que los que no tienen? Dos proxies de la edad del PROMOTOR (no de la
# copia de TE):
#   - turnover (clasificacion de Young et al. 2015: actividad promotora
#     conservada/perdida/ganada entre humano y raton; "human-inserted" =
#     secuencia sin homologo en raton, "mouse-deleted"/"mouse-diminished"
#     = perdida en raton, "expression-turnover" = actividad no conservada);
#   - conservacion (PhyloP 100 vias) en las 3 ventanas del promotor.
# Se usan los 20851 promotores FP de prom_df (no hace falta actividad).
# TE = TE_overlap > 0 (RepeatMasker, sin Low_complexity/Simple_repeat).
#
# Ojo: turnover es NA para mucho mas promotores con TE (61%) que sin TE
# (25%) - probablemente secuencias repetitivas sin ortologo claramente
# asignable - el grafico de barras los excluye (porcentajes sobre los
# clasificados) y el test Fisher de NA se reporta aparte.
#
# Salida en figures/exploratorio_te_edad/. Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)

out_dir <- "figures/exploratorio_te_edad"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>%
  filter(type == "promoter") %>%
  mutate(grupo = factor(ifelse(TE_overlap > 0, "Con TE", "Sin TE"), levels = c("Sin TE", "Con TE")))
cat("Promotores:", table(p$grupo), "\n")

# --- turnover ---------------------------------------------------------------
niv <- c("matched", "ambiguous", "mouse-diminished", "mouse-deleted", "expression-turnover", "human-inserted")
etq <- c(matched = "Conservado (matched)", ambiguous = "Ambiguo", `mouse-diminished` = "Disminuido en ratón",
  `mouse-deleted` = "Delecionado en ratón", `expression-turnover` = "Cambio de expresión", `human-inserted` = "Insertado en humanos")
# Solo promotores con turnover clasificado (se ignoran los NA).
na_pct <- p %>% group_by(grupo) %>% summarise(pct_na = 100 * mean(is.na(turnover)), n_clas = sum(!is.na(turnover)), .groups = "drop")
tab <- p %>%
  filter(!is.na(turnover)) %>%
  mutate(turnover = factor(turnover, levels = niv)) %>%
  count(grupo, turnover) %>%
  group_by(grupo) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
cat("\nTurnover (% sobre los promotores CLASIFICADOS del grupo):\n")
print(tab %>% select(-n) %>% pivot_wider(names_from = grupo, values_from = pct) %>% mutate(across(where(is.numeric), ~ round(.x, 1))))

cl <- p %>% filter(!is.na(turnover))
cat("\nClasificados: ", paste(names(table(cl$grupo)), table(cl$grupo), collapse = " | "), "\n")
cat("Chi-cuadrado turnover (solo clasificados): p =", signif(chisq.test(table(cl$grupo, cl$turnover))$p.value, 3), "\n")
or_tab <- map(c("human-inserted", "mouse-deleted", "mouse-diminished", "expression-turnover", "matched"), function(k) {
  f <- fisher.test(table(cl$grupo, cl$turnover == k))
  tibble(clase = k, pct_sinTE = 100 * mean(cl$turnover[cl$grupo == "Sin TE"] == k), pct_conTE = 100 * mean(cl$turnover[cl$grupo == "Con TE"] == k),
    OR = unname(f$estimate), lo = f$conf.int[1], hi = f$conf.int[2], p = f$p.value)
}) %>% list_rbind()
cat("\nFisher por clase (OR > 1: mas frecuente con TE), solo clasificados:\n")
print(or_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))
cat("\nSin clasificar (NA): Con TE vs Sin TE: Fisher p =", signif(fisher.test(table(p$grupo, is.na(p$turnover)))$p.value, 3), "\n")

# --- conservacion -----------------------------------------------------------
cons <- p %>%
  select(grupo, `16 a -50 pb` = phylop100_close, `-50 a -150 pb` = phylop100_intermediate, `-150 a -235 pb` = phylop100_far) %>%
  pivot_longer(-grupo, names_to = "ventana", values_to = "phylop") %>%
  filter(!is.na(phylop)) %>%
  mutate(ventana = factor(ventana, levels = c("16 a -50 pb", "-50 a -150 pb", "-150 a -235 pb")))
cons_tab <- cons %>% group_by(ventana) %>% summarise(
  mediana_sinTE = median(phylop[grupo == "Sin TE"]), mediana_conTE = median(phylop[grupo == "Con TE"]),
  p = wilcox.test(phylop ~ grupo)$p.value, .groups = "drop")
cat("\nPhyloP (mediana) por ventana:\n")
print(cons_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# --- graficos ---------------------------------------------------------------
pal <- c("Conservado (matched)" = "#216869", "Ambiguo" = "#9DC8C8", "Disminuido en ratón" = "#E8B04B", "Delecionado en ratón" = "#D6741F",
  "Cambio de expresión" = "#A23E48", "Insertado en humanos" = "#5B2A86")
p_t <- tab %>%
  mutate(turnover = factor(etq[as.character(turnover)], levels = etq[niv])) %>%
  ggplot(aes(grupo, pct, fill = turnover)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = pal) +
  theme_pubclean(base_size = 18) +
  theme(legend.position = "right", legend.text = element_text(size = 14), plot.title = element_text(size = 17), plot.caption = element_text(size = 11, hjust = 0), plot.caption.position = "plot") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_x_discrete(labels = setNames(sprintf("%s\n(n = %d)", na_pct$grupo, na_pct$n_clas), as.character(na_pct$grupo))) +
  labs(x = NULL, y = "% de promotores clasificados", fill = NULL, title = "Turnover humano-ratón (Young et al. 2015)",
    caption = sprintf("Se excluyen los promotores sin clasificar: %.0f%% de los Con TE y %.0f%% de los Sin TE", na_pct$pct_na[na_pct$grupo == "Con TE"], na_pct$pct_na[na_pct$grupo == "Sin TE"]))
p_c <- ggplot(cons, aes(ventana, phylop, fill = grupo)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c("Sin TE" = "#9DC8C8", "Con TE" = "#D6741F")) +
  coord_cartesian(ylim = c(-1.5, 3)) +
  theme_pubclean(base_size = 18) +
  theme(plot.title = element_text(size = 17)) +
  labs(x = "Ventana del promotor", y = "PhyloP (100 vías)", fill = NULL, title = "Conservación")
ggsave(file.path(out_dir, "te_vs_sinte_antiguedad.jpg"), ggarrange(p_t, p_c, widths = c(1.25, 1.2)), width = 16, height = 6.5, units = "in")

write_tsv(or_tab, file.path(out_dir, "te_vs_sinte_turnover_fisher.tsv"))
write_tsv(cons_tab, file.path(out_dir, "te_vs_sinte_phylop.tsv"))
