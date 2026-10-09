# cgi_actividad, tata_actividad y ccaat_gcbox_actividad (ver
# docs/mapping_figuras.csv para el numero de figura vigente)
# Efecto de la presencia de islas CpG, TATA-box, CCAAT-box y GC-box sobre la
# actividad media, agrupando promotores en bins de 100 segun su rank de
# actividad. Una figura por motivo (CCAAT-box y GC-box juntos) - antes era
# una sola figura (tata_cgi_actividad, 12 paneles). Las islas CpG se definen
# por la composicion del fragmento (CGI_frag), y cgi_actividad suma un panel
# con la razon CpG o/e como variable continua (cgi_oe_deciles.jpg).
# Tablas de soporte del texto: tata_estratificado_cgi.tsv (tata_actividad),
# ccaat_gcbox_estratificado_cgi.tsv, motivos_prevalencia_cgi.tsv y
# ccaat_gcbox_combinacion.tsv (ccaat_gcbox_actividad).
# Ademas guarda cgi_correlaciones.tsv (Spearman de la actividad con cpg_oe y
# g_c, en total y dentro de cada grupo de CGI_frag, por replica).
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y R/01_activity_stats).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

out_dir_cgi <- fig_dir("cgi_actividad")
out_dir_tata <- fig_dir("tata_actividad")
out_dir_ccaat_gcbox <- fig_dir("ccaat_gcbox_actividad")

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))
# copia sin recortar (los bins descartan los promotores sobrantes de menor
# rank); las tablas de efectos usan todos los promotores, igual que summary_features_secuencia.
data_full <- data

# bins de 100 promotores por actividad media (rank), por separado en cada
# replica; se descartan los de menor rank sobrantes para que el total sea
# multiplo de 100.
excess <- data %>%
  group_by(rep) %>%
  summarise(excess = n() %% 100)
data <- data %>%
  ungroup() %>%
  group_split(rep) %>%
  map2(excess$excess, ~ .x %>%
    mutate(mean_rank_sw = row_number(mean)) %>%
    filter(mean_rank_sw > .y) %>%
    mutate(mean_rank_sw = row_number(mean), mean_sw = factor(ceiling(mean_rank_sw / 100))) %>%
    group_by(rep, mean_sw) %>%
    mutate(var_rank_sw = row_number(var)) %>%
    ungroup()) %>%
  bind_rows()

presence_labels <- c("FALSE" = "Ausencia", "TRUE" = "Presencia")

motif_scatter <- function(data, motif, y_label, title) {
  data %>%
    group_by(mean_sw, rep) %>%
    summarise(prop = sum(.data[[motif]]) / n(), .groups = "drop") %>%
    ggplot(aes(as.numeric(mean_sw), prop)) +
    geom_point(col = "#14AFB2") +
    geom_smooth(col = "#216869") +
    theme_pubclean() +
    theme(text = element_text(size = 30), panel.spacing = unit(3, "lines")) +
    facet_wrap(~rep) +
    ylim(0, 1) +
    scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = c(0.03, 0.1))) +
    labs(x = "Bins de actividad media", y = y_label) +
    ggtitle(title)
}

violin_by_motif <- function(data, rep_id, motif, fill_label) {
  p <- ggviolin(data %>% filter(rep == rep_id),
    x = motif, y = "mean", fill = motif, draw_quantiles = 0.5,
    add = "median_q1q3", palette = c("#14AFB2", "#216869"), alpha = 0.7
  )
  p + stat_compare_means() +
    theme_pubclean() +
    scale_x_discrete(labels = presence_labels) +
    labs(fill = fill_label, x = NULL, y = "Actividad media") +
    theme(legend.position = "none", axis.ticks.x = element_blank(), text = element_text(size = 40)) +
    ylim(c(0.5, 6.2))
}

# --- TATA-box ---------------------------------------------------------

panel_tata_scatter <- motif_scatter(data, "TATA_EPD", "Proporción de promotores\ncon TATA-box", "TATA-box")
ggsave(file.path(out_dir_tata, "tata_scatter.jpg"), panel_tata_scatter, width = 12, height = 6.75, units = "in")

ggsave(file.path(out_dir_tata, "tata_violin_rep1.jpg"), violin_by_motif(data, "Rep 1", "TATA_EPD", "TATA-box"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir_tata, "tata_violin_rep2.jpg"), violin_by_motif(data, "Rep 2", "TATA_EPD", "TATA-box"), width = 9, height = 6.75, units = "in")

# --- Isla CpG -----------------------------------------------------------
# Usa CGI_frag (el fragmento de 252 pb cumple G+C >= 50% y CpG o/e > 0.6),
# no la anotacion genomica de UCSC (CGI); ver R/00_prom_features/build_prom_features.R.

panel_cgi_scatter <- motif_scatter(data, "CGI_frag", "Proporción de promotores\ncon islas CpG", "Isla CpG")
ggsave(file.path(out_dir_cgi, "cgi_scatter.jpg"), panel_cgi_scatter, width = 12, height = 6.75, units = "in")

ggsave(file.path(out_dir_cgi, "cgi_violin_rep1.jpg"), violin_by_motif(data, "Rep 1", "CGI_frag", "Isla CpG"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir_cgi, "cgi_violin_rep2.jpg"), violin_by_motif(data, "Rep 2", "CGI_frag", "Isla CpG"), width = 9, height = 6.75, units = "in")

# Razon CpG o/e como variable continua: percentil de actividad media
# (dentro de cada replica) por decil de cpg_oe. Punto = media del decil,
# barra = IC 95%; linea punteada vertical = umbral 0.6 de CGI_frag.
cgi_oe_data <- data %>%
  filter(!is.na(cpg_oe)) %>%
  group_by(rep) %>%
  mutate(act_pct = 100 * percent_rank(mean), oe_decil = ntile(cpg_oe, 10)) %>%
  ungroup()
cgi_oe_rho <- cgi_oe_data %>%
  group_by(rep) %>%
  summarise(rho = cor(cpg_oe, mean, method = "spearman"), .groups = "drop") %>%
  mutate(label = sprintf("rho == %.2f", rho))
panel_cgi_oe_deciles <- cgi_oe_data %>%
  group_by(rep, oe_decil) %>%
  summarise(oe = median(cpg_oe), m = mean(act_pct), se = sd(act_pct) / sqrt(n()), .groups = "drop") %>%
  ggplot(aes(oe, m)) +
  geom_hline(yintercept = 50, linetype = "dotted", col = "grey60") +
  geom_vline(xintercept = 0.6, linetype = "dashed", col = "grey40") +
  geom_errorbar(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), width = 0, col = "#216869") +
  geom_line(col = "#216869") +
  geom_point(col = "#14AFB2", size = 3) +
  geom_text(data = cgi_oe_rho, aes(x = -Inf, y = Inf, label = label), parse = TRUE,
    hjust = -0.2, vjust = 1.5, size = 9, inherit.aes = FALSE) +
  theme_pubclean() +
  theme(text = element_text(size = 30), panel.spacing = unit(3, "lines")) +
  facet_wrap(~rep) +
  labs(x = "Razón CpG o/e del fragmento (mediana del decil)", y = "Percentil de\nactividad media") +
  ggtitle("Razón CpG observado/esperado")
ggsave(file.path(out_dir_cgi, "cgi_oe_deciles.jpg"), panel_cgi_oe_deciles, width = 12, height = 6.75, units = "in")

# Correlaciones de Spearman de la actividad media con la razon CpG o/e y con
# el contenido de G+C del fragmento, por replica: en todos los promotores y
# por separado dentro de los que cumplen / no cumplen CGI_frag. Sostienen
# el texto de Resultados (efecto gradual de CpG o/e incluso dentro de las
# islas; G+C con la mitad de correlacion y ~0 dentro de las islas). Se
# calculan sobre el mismo `data` que los paneles (bins recortados).
cor_spearman <- function(df, var) {
  ct <- suppressWarnings(cor.test(df[[var]], df$mean, method = "spearman", exact = FALSE))
  tibble(rho = unname(ct$estimate), pval = ct$p.value)
}
cgi_cor <- data %>%
  filter(!is.na(cpg_oe)) %>%
  mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")) %>%
  bind_rows(mutate(., grupo = "Todos")) %>%
  group_by(rep, grupo) %>%
  reframe(
    n = n(),
    bind_rows(
      mutate(cor_spearman(pick(everything()), "cpg_oe"), variable = "cpg_oe"),
      mutate(cor_spearman(pick(everything()), "g_c"), variable = "g_c")
    )
  ) %>%
  mutate(grupo = factor(grupo, levels = c("Todos", "Con isla CpG", "Sin isla CpG"))) %>%
  arrange(rep, variable, grupo) %>%
  select(rep, variable, grupo, n, rho, pval)
write_tsv(cgi_cor, file.path(out_dir_cgi, "cgi_correlaciones.tsv"))
print(cgi_cor, n = Inf)

# --- CCAAT-box ------------------------------------------------------------

panel_ccaat_scatter <- motif_scatter(data, "CCAAT_EPD", "Proporción de promotores\ncon CCAAT-box", "CCAAT-box")
ggsave(file.path(out_dir_ccaat_gcbox, "ccaat_scatter.jpg"), panel_ccaat_scatter, width = 12, height = 6.75, units = "in")

ggsave(file.path(out_dir_ccaat_gcbox, "ccaat_violin_rep1.jpg"), violin_by_motif(data, "Rep 1", "CCAAT_EPD", "CCAAT-box"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir_ccaat_gcbox, "ccaat_violin_rep2.jpg"), violin_by_motif(data, "Rep 2", "CCAAT_EPD", "CCAAT-box"), width = 9, height = 6.75, units = "in")

# --- GC-box ---------------------------------------------------------------

panel_gcbox_scatter <- motif_scatter(data, "GCbox_EPD", "Proporción de promotores\ncon GC-box", "GC-box")
ggsave(file.path(out_dir_ccaat_gcbox, "gcbox_scatter.jpg"), panel_gcbox_scatter, width = 12, height = 6.75, units = "in")

ggsave(file.path(out_dir_ccaat_gcbox, "gcbox_violin_rep1.jpg"), violin_by_motif(data, "Rep 1", "GCbox_EPD", "GC-box"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir_ccaat_gcbox, "gcbox_violin_rep2.jpg"), violin_by_motif(data, "Rep 2", "GCbox_EPD", "GC-box"), width = 9, height = 6.75, units = "in")

# --- Tablas de soporte para el texto de Resultados -------------------------
# Efecto (Wilcoxon + estimador de Hodges-Lehmann con IC 95%) de TATA-box,
# CCAAT-box y GC-box sobre la actividad media, en todos los promotores y por
# separado con / sin isla CpG (CGI_frag), por replica; prevalencia de cada
# motivo segun CGI_frag en la library completa; y efecto de las
# combinaciones CCAAT-box x GC-box frente a promotores sin ninguno de los dos.
# Usan data_full (todos los promotores, sin el recorte de los bins).

hl_test <- function(x, y) {
  w <- wilcox.test(x, y, conf.int = TRUE)
  tibble(HL = unname(w$estimate), IC_inf = w$conf.int[1], IC_sup = w$conf.int[2], pval = w$p.value)
}

motif_cols <- c("TATA-box" = "TATA_EPD", "CCAAT-box" = "CCAAT_EPD", "GC-box" = "GCbox_EPD")

motivos_estratificado_cgi <- imap(motif_cols, function(col, motivo) {
  data_full %>%
    mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")) %>%
    bind_rows(mutate(., grupo = "Todos")) %>%
    group_by(rep, grupo) %>%
    reframe(
      n = n(), n_motivo = sum(.data[[col]]), pct_motivo = 100 * mean(.data[[col]]),
      hl_test(mean[.data[[col]]], mean[!.data[[col]]])
    ) %>%
    mutate(motivo = motivo)
}) %>%
  list_rbind() %>%
  mutate(grupo = factor(grupo, levels = c("Todos", "Con isla CpG", "Sin isla CpG"))) %>%
  arrange(motivo, rep, grupo) %>%
  select(motivo, rep, grupo, n, n_motivo, pct_motivo, HL, IC_inf, IC_sup, pval)

motivos_prevalencia_cgi <- prom_df %>%
  filter(!is.na(CGI_frag)) %>%
  mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")) %>%
  group_by(grupo) %>%
  summarise(n = n(), across(all_of(unname(motif_cols)), ~ 100 * mean(.x), .names = "pct_{.col}"))

ccaat_gcbox_combinacion <- data_full %>%
  mutate(combinacion = case_when(
    CCAAT_EPD & GCbox_EPD ~ "Ambos",
    CCAAT_EPD ~ "Solo CCAAT-box",
    GCbox_EPD ~ "Solo GC-box",
    TRUE ~ "Ninguno"
  )) %>%
  group_by(rep) %>%
  group_modify(function(df, ...) {
    map(c("Solo CCAAT-box", "Solo GC-box", "Ambos"), function(k) {
      tibble(combinacion = k, n = sum(df$combinacion == k), hl_test(df$mean[df$combinacion == k], df$mean[df$combinacion == "Ninguno"]))
    }) %>% list_rbind()
  }) %>%
  ungroup()

write_tsv(filter(motivos_estratificado_cgi, motivo == "TATA-box"), file.path(out_dir_tata, "tata_estratificado_cgi.tsv"))
write_tsv(filter(motivos_estratificado_cgi, motivo != "TATA-box"), file.path(out_dir_ccaat_gcbox, "ccaat_gcbox_estratificado_cgi.tsv"))
write_tsv(motivos_prevalencia_cgi, file.path(out_dir_ccaat_gcbox, "motivos_prevalencia_cgi.tsv"))
write_tsv(ccaat_gcbox_combinacion, file.path(out_dir_ccaat_gcbox, "ccaat_gcbox_combinacion.tsv"))
print(motivos_estratificado_cgi, n = Inf)
print(motivos_prevalencia_cgi)
print(ccaat_gcbox_combinacion)

message("Figuras guardadas en ", paste(c(out_dir_cgi, out_dir_tata, out_dir_ccaat_gcbox), collapse = ", "))
