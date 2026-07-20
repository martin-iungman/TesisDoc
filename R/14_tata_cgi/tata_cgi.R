# tata_cgi_actividad (ver docs/mapping_figuras.csv para el numero de
# figura vigente - a la fecha sin resolver, ver notas del CSV)
# Efecto de la presencia de TATA-box e islas CpG sobre la actividad media,
# agrupando promotores en bins de 100 segun su rank de actividad.
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y R/01_activity_stats).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "tata_cgi_actividad"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

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

violin_by_motif <- function(data, rep_id, motif, fill_label) {
  p <- ggviolin(data %>% filter(rep == rep_id),
    x = motif, y = "mean", fill = motif, draw_quantiles = 0.5,
    add = "median_q1q3", palette = c("#14AFB2", "#216869"), alpha = 0.7
  )
  p + stat_compare_means() +
    theme_pubclean() +
    labs(fill = fill_label, x = NULL, y = "Actividad media") +
    theme(legend.position = "none", axis.ticks.x = element_blank(), text = element_text(size = 40)) +
    ylim(c(0.5, 6.2))
}

# --- TATA-box ---------------------------------------------------------

panel_tata_scatter <- data %>%
  group_by(mean_sw, rep) %>%
  summarise(TATA_EPD = sum(TATA_EPD) / n(), .groups = "drop") %>%
  ggplot(aes(as.numeric(mean_sw), TATA_EPD)) +
  geom_point(col = "#14AFB2") +
  geom_smooth(col = "#216869") +
  theme_pubclean() +
  theme(text = element_text(size = 30)) +
  facet_wrap(~rep) +
  ylim(0, 1) +
  labs(x = "Actividad media (agrupadas por rank)", y = "Proporción de promotores\ncon TATA-box") +
  ggtitle("TATA-box")
ggsave(file.path(out_dir, "tata_scatter.jpg"), panel_tata_scatter, width = 9, height = 6.75, units = "in")

ggsave(file.path(out_dir, "tata_violin_rep1.jpg"), violin_by_motif(data, "Rep 1", "TATA_EPD", "TATA-box"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "tata_violin_rep2.jpg"), violin_by_motif(data, "Rep 2", "TATA_EPD", "TATA-box"), width = 9, height = 6.75, units = "in")

# --- Isla CpG -----------------------------------------------------------

panel_cgi_scatter <- data %>%
  group_by(mean_sw, rep) %>%
  summarise(CGI = sum(CGI) / n(), .groups = "drop") %>%
  ggplot(aes(as.numeric(mean_sw), CGI)) +
  geom_point(col = "#14AFB2") +
  geom_smooth(col = "#216869") +
  theme_pubclean() +
  theme(text = element_text(size = 30)) +
  facet_wrap(~rep) +
  ylim(0, 1) +
  labs(x = "Bins de actividad media", y = "Proporción de promotores\ncon islas CpG") +
  ggtitle("Isla CpG")
ggsave(file.path(out_dir, "cgi_scatter.jpg"), panel_cgi_scatter, width = 9, height = 6.75, units = "in")

ggsave(file.path(out_dir, "cgi_violin_rep1.jpg"), violin_by_motif(data, "Rep 1", "CGI", "Isla CpG"), width = 9, height = 6.75, units = "in")
ggsave(file.path(out_dir, "cgi_violin_rep2.jpg"), violin_by_motif(data, "Rep 2", "CGI", "Isla CpG"), width = 9, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
