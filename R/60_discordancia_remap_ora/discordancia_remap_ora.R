# discordancia_remap_ora (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Efecto de la union de cada factor de transcripcion (ChIP-seq ReMap2022,
# cualquier linea celular) sobre la discordancia con signo entre el reportero
# y la actividad endogena maxima (FANTOM5), en promotores activos en HEK293
# (hek_tpm > 0), con el mismo modelo que discordancia_features:
#
#   dif_signed ~ avg_rank + rep + TF
#
# con errores estandar robustos agrupados por promotor (ajustar_discordancia_
# robusto() en R/functions/plot_helpers.R) y BH sobre los TF.
#
# dif_signed = rr - re: positivo = el reportero SOBREESTIMA, negativo = la
# SUBESTIMA (ver R/59_discordancia_features/discordancia_features.R para la
# convencion y el diseno). Se testean los TF con mas de 100 promotores unidos
# en ambas replicas.
#
# Criterio de robustez: el mismo modelo con la union restringida a picos no
# redundantes de ReMap a los que contribuye HEK293 (str_detect "HEK293": HEK293,
# HEK293T, HEK293-FT, igual que tf_contexto_endogeno.R). Un TF es concordante si
# tiene FDR < 0.05 y el mismo signo en ambos analisis. El orden es por el menor
# de los dos efectos absolutos. Salida: tables/discordancia_remap_robustez_hek.tsv.
#
# Sobrerrepresentacion de terminos GO (ORA, enrichGO, universo = TF testeados)
# entre los TF asociados a cada direccion. La figura usa FDR < 0.05; la tabla
# trae tambien FDR < 0.001. Los nombres de los terminos GO quedan en ingles,
# como en la ontologia.
#
# Paneles:
# (A) efecto de cada TF; (B) terminos GO por direccion; (C) discordancia segun
# la union de tres correpresores (SETDB1, MGA, HDAC3), con una observacion
# por promotor (promedio de sus replicas) y Wilcoxon sobre esas observaciones.
# Los tres siguen significativos en ambas replicas despues de corregir por
# cobertura (ver R/65_exploratorio_cobertura_discordancia). SETDB1 ademas es
# concordante con la union en HEK293; MGA y HDAC3 no tienen suficientes picos
# de HEK293 para evaluarlo.
#
# Requiere: data/processed/remap_tf_hits.tsv, activity_stats_highconf.tsv,
# prom_df.tsv y fantom_endo_activity_summary.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

select <- dplyr::select
filter <- dplyr::filter

slug <- "discordancia_remap_ora"
out_dir <- fig_dir(slug)

marcadores <- c("SETDB1", "MGA", "HDAC3")
direcciones <- c("sobreestima" = "Reportero sobreestima", "subestima" = "Reportero subestima")

save_panel <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(
      file.path(out_dir, paste0(name, ".", ext)), plot,
      width = width, height = height, units = "in",
      device = if (ext == "pdf") cairo_pdf else "jpeg"
    )
  }
}

remap_todos <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
remap_hits <- remap_todos %>% distinct(seq_id, TF)
remap_hits_hek <- remap_todos %>% filter(str_detect(sample, "HEK293")) %>% distinct(seq_id, TF)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  add_rank_discordance() %>%
  select(seq_id, rep, avg_rank, dif_signed)

message("n = ", n_distinct(data$seq_id), " promotores activos en HEK293")

# --- Modelo por TF ----------------------------------------------------------------

testear_tfs <- function(hits) {
  hits_data <- hits %>% semi_join(data, by = "seq_id")
  tfs <- data %>%
    inner_join(hits_data, by = "seq_id", relationship = "many-to-many") %>%
    count(rep, TF) %>%
    filter(n > 100) %>%
    count(TF, name = "n_rep") %>%
    filter(n_rep == 2, !str_detect(TF, "^H\\d")) %>%
    pull(TF)
  message(length(tfs), " TF con mas de 100 promotores en ambas replicas")

  map(tfs, function(tf) {
    unidos <- hits_data$seq_id[hits_data$TF == tf]
    data %>%
      mutate(presente = factor(seq_id %in% unidos, levels = c(FALSE, TRUE))) %>%
      ajustar_discordancia_robusto() %>%
      transmute(TF = tf, estimate, conf.low, conf.high, std.error, statistic, df, p_value)
  }) %>%
    list_rbind() %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)
}

resultados <- testear_tfs(remap_hits)

write_tsv(resultados, "tables/discordancia_remap_lmm.tsv")
message(
  nrow(resultados), " TF: ", sum(resultados$p_adj < 0.05 & resultados$estimate > 0), " sobreestima / ",
  sum(resultados$p_adj < 0.05 & resultados$estimate < 0), " subestima (FDR < 0.05)"
)
message("Marcadores:")
print(resultados %>% filter(TF %in% marcadores) %>% select(TF, estimate, conf.low, conf.high, p_adj))

# --- Criterio de robustez: union en HEK293 ------------------------------------------

robustez <- resultados %>%
  transmute(TF, estimate_todos = estimate, p_adj_todos = p_adj) %>%
  inner_join(
    testear_tfs(remap_hits_hek) %>% transmute(TF, estimate_hek = estimate, p_adj_hek = p_adj),
    by = "TF"
  ) %>%
  mutate(
    concordante = p_adj_todos < 0.05 & p_adj_hek < 0.05 & sign(estimate_todos) == sign(estimate_hek),
    direccion = case_when(
      !concordante ~ NA_character_,
      estimate_todos > 0 ~ "sobreestima",
      TRUE ~ "subestima"
    ),
    efecto_minimo = pmin(abs(estimate_todos), abs(estimate_hek))
  ) %>%
  group_by(direccion) %>%
  arrange(desc(efecto_minimo), .by_group = TRUE) %>%
  mutate(puesto = if_else(concordante, row_number(), NA_integer_)) %>%
  ungroup()

write_tsv(robustez, "tables/discordancia_remap_robustez_hek.tsv")
message(
  nrow(robustez), " TF con union en HEK293: ", sum(robustez$direccion == "sobreestima", na.rm = TRUE), " sobreestima / ",
  sum(robustez$direccion == "subestima", na.rm = TRUE), " subestima concordantes"
)
print(robustez %>% filter(TF %in% marcadores) %>% select(TF, estimate_todos, estimate_hek, direccion, puesto))

# --- (A) Efecto por TF --------------------------------------------------------------

colores_fdr <- c("FDR < 0.001" = "#c0392b", "FDR < 0.01" = "#e67e22", "FDR < 0.05" = "#2980b9", "ns" = "grey70")

p_tf <- resultados %>%
  mutate(
    sig = case_when(
      p_adj < 0.001 ~ "FDR < 0.001",
      p_adj < 0.01 ~ "FDR < 0.01",
      p_adj < 0.05 ~ "FDR < 0.05",
      TRUE ~ "ns"
    ) %>% factor(levels = names(colores_fdr)),
    TF = fct_reorder(TF, estimate)
  ) %>%
  ggplot(aes(x = estimate, y = TF, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = colores_fdr) +
  labs(
    x = "Efecto (β) (positivo: el reportero sobreestima)",
    y = paste0("Factores de transcripción (ReMap, n = ", nrow(resultados), ")"),
    color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

# --- (B) ORA ----------------------------------------------------------------------

ids <- AnnotationDbi::select(org.Hs.eg.db, keys = resultados$TF, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct() %>%
  filter(!is.na(ENTREZID))

correr_ora <- function(direccion, fdr) {
  genes <- resultados %>%
    filter(p_adj < fdr, if (direccion == "sobreestima") estimate > 0 else estimate < 0) %>%
    inner_join(ids, by = c("TF" = "SYMBOL")) %>%
    pull(ENTREZID)
  ora <- enrichGO(
    gene = genes, universe = ids$ENTREZID, OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID", ont = "all", pAdjustMethod = "BH"
  )
  as_tibble(as.data.frame(ora)) %>%
    arrange(p.adjust) %>%
    mutate(direccion = direccion, fdr_tf = fdr, n_tf = length(genes))
}

ora <- expand_grid(direccion = names(direcciones), fdr = c(0.05, 0.001)) %>%
  mutate(res = map2(direccion, fdr, correr_ora))

print(ora %>% transmute(direccion, fdr, n_tf = map_dbl(res, ~ .x$n_tf[1]), n_terminos = map_int(res, nrow)))

write_xlsx(
  set_names(ora$res, paste0(ora$direccion, "_fdr", ora$fdr)),
  "tables/discordancia_remap_ora.xlsx"
)

terminos <- ora %>%
  filter(fdr == 0.05) %>%
  pull(res) %>%
  map(~ slice_head(.x, n = 12)) %>%
  list_rbind() %>%
  mutate(
    direccion = factor(direcciones[direccion], levels = direcciones),
    neglog10_p = -log10(p.adjust),
    termino = paste0(str_wrap(Description, width = 45), "  ")
  ) %>%
  group_by(direccion) %>%
  mutate(termino = fct_reorder(termino, neglog10_p)) %>%
  ungroup()

p_ora <- ggplot(terminos, aes(x = neglog10_p, y = termino, fill = Count)) +
  geom_col() +
  facet_grid(direccion ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(width = 12)) +
  scale_fill_gradient(low = "#7FB800", high = "#0D2C54") +
  labs(x = "-log10(p ajustado)", y = NULL, fill = "Número de TF") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = rel(0.8), lineheight = 0.85),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )

# --- (C) Correpresores --------------------------------------------------------------

por_marcador <- map(marcadores, function(m) {
  data %>% mutate(marcador = m, union = ifelse(seq_id %in% remap_hits$seq_id[remap_hits$TF == m], "Presente", "Ausente"))
}) %>%
  list_rbind() %>%
  mutate(marcador = factor(marcador, levels = marcadores))

una_por_promotor <- por_marcador %>%
  group_by(marcador, seq_id, union) %>%
  summarise(dif_signed = mean(dif_signed), .groups = "drop")

etiquetas_p <- una_por_promotor %>%
  group_by(marcador) %>%
  summarise(p = wilcox.test(dif_signed ~ union)$p.value, .groups = "drop") %>%
  mutate(etiqueta = ifelse(p < 2.2e-16, "p < 2.2e-16", sprintf("p = %.1e", p)))

p_violin <- una_por_promotor %>%
  ggplot(aes(union, dif_signed, fill = union)) +
  geom_violin(draw_quantiles = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = etiquetas_p, aes(x = 1.5, y = 1.2, label = etiqueta), inherit.aes = FALSE, size = 4) +
  coord_cartesian(ylim = c(-1, 1.3)) +
  facet_wrap(~marcador, ncol = 3) +
  scale_fill_manual(values = c("Presente" = "#c0392b", "Ausente" = "grey70")) +
  labs(x = "Unión del factor", y = "Discordancia\n(reportero - endógeno)", fill = NULL) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# --- Figura completa y paneles sueltos ----------------------------------------------

save_panel(p_tf, "discordancia_remap_tf", 9, 6.75)
save_panel(p_ora, "discordancia_remap_ora_terminos", 11, 9.5)
save_panel(p_violin, "discordancia_remap_represores", 11, 4.5)

figura <- ((p_tf | p_ora) + plot_layout(widths = c(1, 1.15))) / p_violin +
  plot_layout(heights = c(1.5, 0.7)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))
save_panel(figura, slug, 15, 13)

message("Figuras guardadas en ", out_dir)
