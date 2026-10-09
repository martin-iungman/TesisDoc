# Comparacion del efecto sobre la actividad media del reportero entre el
# MOTIVO (EPD) y la UNION DEL TF correspondiente (ChIP-seq de ReMap2022):
# CCAAT-box vs. NFYA, GC-box vs. SP1/SP2 (union en cualquier linea celular
# y restringida a HEK293 para SP1/SP2). Mismo estadistico que summary_features_secuencia: efecto
# de Wilcoxon (estimador de Hodges-Lehmann, mean ~ feature, IC95% por
# replica) via wilcox_effect_summary() - no se filtra por significancia,
# para ver los 5-7 features lado a lado.
#
# Sale en la carpeta de tf_nfya_sp_contexto_endogeno, no es una
# figura nueva en mapping_figuras.csv.
#
# Requiere: data/processed/remap_tf_hits.tsv, activity_stats_highconf.tsv
# y prom_df.tsv. Run from the TesisDoc repo root.

library(tidyverse)
library(coin)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

out_dir <- fig_dir("tf_nfya_sp_contexto_endogeno")

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

tf_flag <- function(tf, hek_only = FALSE) {
  hits <- remap_hits %>% filter(TF == tf)
  if (hek_only) hits <- hits %>% filter(str_detect(sample, "HEK293"))
  data$seq_id %in% hits$seq_id
}

features <- tibble(
  seq_id = data$seq_id, rep = data$rep, mean = data$mean,
  `CCAAT-box` = data$CCAAT_EPD,
  `NFYA (ReMap)` = tf_flag("NFYA"),
  `GC-box` = data$GCbox_EPD,
  `SP1 (ReMap)` = tf_flag("SP1"),
  `SP2 (ReMap)` = tf_flag("SP2"),
  `SP1 (ReMap, HEK293)` = tf_flag("SP1", hek_only = TRUE),
  `SP2 (ReMap, HEK293)` = tf_flag("SP2", hek_only = TRUE)
)
feature_levels <- setdiff(names(features), c("seq_id", "rep", "mean"))

# n por feature y solapamiento motivo-TF (para interpretar el grafico)
overlap <- features %>%
  filter(rep == "Rep 1") %>%
  summarise(
    n = n(),
    CCAAT = sum(`CCAAT-box`), NFYA = sum(`NFYA (ReMap)`), CCAAT_y_NFYA = sum(`CCAAT-box` & `NFYA (ReMap)`),
    GCbox = sum(`GC-box`), SP1 = sum(`SP1 (ReMap)`), SP2 = sum(`SP2 (ReMap)`),
    GCbox_y_SP1 = sum(`GC-box` & `SP1 (ReMap)`), GCbox_y_SP2 = sum(`GC-box` & `SP2 (ReMap)`)
  )
print(overlap)

tidy_data <- features %>%
  select(-seq_id) %>%
  mutate(
    rep = as.factor(rep),
    across(all_of(feature_levels), ~ factor(.x, levels = c("TRUE", "FALSE")))
  )

wilcox_df <- wilcox_effect_summary(tidy_data)
write_tsv(wilcox_df, file.path(out_dir, "motivo_vs_tf_efecto.tsv"))

p <- wilcox_df %>%
  pivot_wider(names_from = val, values_from = estimate) %>%
  mutate(
    feature = factor(as.character(feature), levels = rev(feature_levels)),
    grupo = ifelse(str_detect(feature, "CCAAT|NFYA"), "CCAAT / NFYA", "GC-box / SP1 / SP2"),
    sig = ifelse(pval_corr < 0.05, "FDR < 0.05", "ns")
  ) %>%
  ggplot(aes(x = estimate, y = feature, group = fct_rev(rep))) +
  geom_col(orientation = "y", position = "dodge", aes(fill = grupo, alpha = rep)) +
  geom_errorbarh(aes(xmax = P2.5, xmin = P97.5), position = position_dodge(0.9), height = 0.15, col = "#777777", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.6)) +
  scale_fill_manual(values = c("CCAAT / NFYA" = "#D6741F", "GC-box / SP1 / SP2" = "#216869")) +
  labs(x = "Efecto sobre la actividad (Wilcoxon)", y = NULL, fill = NULL, alpha = NULL) +
  theme_pubr(base_size = 18) +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "motivo_vs_tf_efecto.jpg"), p, width = 10, height = 6, units = "in")

message("Figura guardada en ", out_dir)
