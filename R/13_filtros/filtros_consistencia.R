# filtros_consistencia_replicas (ver docs/mapping_figuras.csv para el
# numero de figura vigente)
# Panel A: correlacion de actividad media entre replicas, antes de
#   excluir promotores bimodales/inconsistentes (coloreado por si el
#   promotor termina en el subconjunto de alta confianza).
# Paneles B-D: solo secuencias de tipo promotor (seq_data$type ==
#   "promoter"; se excluyen enhancers y cancer_wt) - el subconjunto de alta
#   confianza que usa la tesis (12284 promotores).
# Panel B: la misma correlacion, ya solo sobre el subconjunto filtrado.
# Panel C: promotores de alta confianza por replica. Dos representaciones
#   (n por replica en barras / venn contra la library de promotores),
#   portadas de transcriptional_library/Tesis/tesis.R (n_per_rep.jpg y
#   venn_rep_postfilter.jpg).
# Panel D: correlacion de la VARIANZA (ruido) entre replicas, sobre el
#   subconjunto de alta confianza - analogo al panel B pero para var en
#   vez de mean. Portado de transcriptional_library/Analysis/scripts/
#   final_github.R ("Variance correlation" / var_replicates.pdf).
#
# Requiere: data/processed/activity_stats_full.tsv,
# data/processed/activity_stats_highconf.tsv y data/processed/data_long.tsv
# (ver R/01_activity_stats/build_activity_stats.R) y data/raw/seq_data.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggExtra)
library(patchwork)
source("R/functions/fig_paths.R")
source("R/functions/venn_library.R")

slug <- "filtros_consistencia_replicas"
out_dir <- fig_dir(slug)

full_stats <- read_tsv("data/processed/activity_stats_full.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
spike_id <- c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")
all_seq <- read_tsv("data/raw/seq_data.tsv", show_col_types = FALSE) %>% filter(!seq_id %in% dupl_id)
promoter_ids <- all_seq %>% filter(type == "promoter") %>% pull(seq_id)
highconf_prom <- stats_highconf %>% filter(seq_id %in% promoter_ids)

save_jpg_pdf <- function(name, plot, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in")
  }
}

# --- Panel A: correlacion entre replicas, pre-filtro (coloreada) ---------

panel_a_base <- full_stats %>%
  filter(sum_norm_counts > 1000) %>%
  select(seq_id, mean, rep) %>%
  pivot_wider(values_from = "mean", names_from = "rep") %>%
  mutate(Confianza = ifelse(seq_id %in% stats_highconf$seq_id, "Alta", "Baja") %>% fct_rev()) %>%
  arrange(desc(Confianza))

p_a <- panel_a_base %>%
  ggplot(aes(`Rep 1`, `Rep 2`, col = Confianza)) +
  geom_point(size = 0.5, alpha = 0.7) +
  labs(x = "Actividad media (Rep 1)", y = "Actividad media (Rep 2)") +
  scale_color_manual(values = c("Alta" = "#1B8C8E", "Baja" = "#0D2C54")) +
  theme_bw() +
  theme(text = element_text(size = 25)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
p_a <- ggExtra::ggMarginal(p_a, type = "density", margins = "both", col = "#0D2C54")
save_jpg_pdf("panel_a_correlacion_prefiltro", p_a, 9, 6.75)

cor_prefiltro <- cor.test(panel_a_base$`Rep 1`, panel_a_base$`Rep 2`, method = "pearson")
message("Pearson pre-filtro: ", round(cor_prefiltro$estimate, 3), ", p = ", format.pval(cor_prefiltro$p.value, digits = 3))

# --- Panel B: correlacion entre replicas, post-filtro ----------------------

panel_b_base <- highconf_prom %>%
  filter(sum_norm_counts > 1000) %>%
  select(seq_id, mean, rep) %>%
  pivot_wider(values_from = "mean", names_from = "rep")

p_b <- panel_b_base %>%
  ggplot(aes(`Rep 1`, `Rep 2`)) +
  geom_point(size = 0.5, alpha = 0.5, col = "#1B8C8E") +
  labs(x = "Actividad media (Rep 1)", y = "Actividad media (Rep 2)") +
  theme_bw() +
  theme(text = element_text(size = 25))
p_b <- ggExtra::ggMarginal(p_b, type = "density", margins = "both", col = "#0D2C54")
save_jpg_pdf("panel_b_correlacion_postfiltro", p_b, 7, 6.75)

cor_postfiltro <- cor.test(panel_b_base$`Rep 1`, panel_b_base$`Rep 2`, method = "pearson")
message("Pearson post-filtro: ", round(cor_postfiltro$estimate, 3), ", p = ", format.pval(cor_postfiltro$p.value, digits = 3))

# --- Panel C candidato 1: N de promotores detectados por replica ---------

data_long <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE)
n_lib <- length(promoter_ids)

count_by_rep_combo <- function(seq_id_rep_df, total_n) {
  tmp <- seq_id_rep_df %>%
    mutate(rep = factor("1|2")) %>%
    bind_rows(seq_id_rep_df) %>%
    distinct(seq_id, rep)
  tmp %>%
    count(seq_id) %>%
    filter(n == 3) %>%
    mutate(rep = factor("1&2")) %>%
    bind_rows(tmp %>% distinct(seq_id, rep)) %>%
    count(rep) %>%
    mutate(total_n = total_n)
}

p1 <- count_by_rep_combo(data_long %>% filter(counts > 0, seq_id %in% promoter_ids) %>% distinct(seq_id, rep), n_lib) %>%
  ggplot(aes(rep, n, fill = rep, label = n)) +
  geom_col(aes(y = total_n), fill = "#dfedf6") +
  geom_col(aes(y = n), fill = "#358AAA") +
  geom_label(nudge_y = 1800, aes(y = n), fill = "#AADAD4") +
  ggtitle("Totales observadas") +
  xlab("Réplica") +
  ylab("Frecuencia") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15))

p2 <- count_by_rep_combo(highconf_prom %>% distinct(seq_id, rep), n_lib) %>%
  ggplot(aes(rep, n, fill = rep, label = n)) +
  geom_col(aes(y = total_n), fill = "#dfedf6") +
  geom_col(aes(y = n), fill = "#358AAA") +
  geom_label(nudge_y = 1800, aes(y = n), fill = "#AADAD4") +
  ggtitle("Post-filtros") +
  xlab("Réplica") +
  ylab("Frecuencia") +
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 15))

save_jpg_pdf("panel_c_candidato1_n_por_replica", p1 + p2, 9, 6.75)

# --- Panel C candidato 2: Venn post-filtro ---------------------------------

highconf_rep <- function(r) {
  highconf_prom %>% filter(rep == r) %>% pull(seq_id) %>% unique()
}
p_c <- venn_library(highconf_rep("Rep 1"), highconf_rep("Rep 2"), promoter_ids, "Réplica 1", "Réplica 2")
message("Promotores de alta confianza: ", n_distinct(highconf_prom$seq_id))
for (ext in c("png", "pdf")) {
  ggsave(file.path(out_dir, paste0("panel_c_candidato2_venn_postfiltro.", ext)), p_c, width = 4.3, height = 4.3, units = "in", dpi = 600, bg = "white")
}

# --- Panel D: correlacion de varianza (ruido) entre replicas -------------

panel_d_base <- highconf_prom %>%
  select(seq_id, var, rep) %>%
  pivot_wider(values_from = "var", names_from = "rep")

p_d <- panel_d_base %>%
  ggplot(aes(`Rep 1`, `Rep 2`)) +
  geom_point(size = 0.5, alpha = 0.3, col = "#1B8C8E") +
  labs(x = "Varianza (Rep 1)", y = "Varianza (Rep 2)") +
  theme_bw() +
  theme(text = element_text(size = 25))
p_d <- ggExtra::ggMarginal(p_d, type = "density", margins = "both", col = "#0D2C54")
save_jpg_pdf("panel_d_var_replicates", p_d, 7, 6.75)

cor_var <- cor.test(panel_d_base$`Rep 1`, panel_d_base$`Rep 2`, method = "pearson")
message("Pearson varianza (post-filtro): ", round(cor_var$estimate, 3), ", p = ", format.pval(cor_var$p.value, digits = 3))

# --- Figura completa: A-D en una sola imagen (panel C = venn) ------------

p_c_grande <- venn_library(highconf_rep("Rep 1"), highconf_rep("Rep 2"), promoter_ids, "Réplica 1", "Réplica 2", text_scale = 1.6)
# fila vacia intermedia para separar las dos filas de paneles
figura <- wrap_elements(full = p_a) + wrap_elements(full = p_b) +
  plot_spacer() + plot_spacer() +
  wrap_elements(full = p_c_grande) + wrap_elements(full = p_d) +
  plot_layout(ncol = 2, widths = c(9, 7), heights = c(1, 0.06, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 32, face = "bold"))
save_jpg_pdf(slug, figura, 16, 14.2)

message("Paneles guardados en ", out_dir, ". Panel C tiene 2 representaciones (barras/venn); la figura completa usa el venn.")
