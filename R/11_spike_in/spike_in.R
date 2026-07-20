# spike_in_robustez (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# Robustez de los controles internos de spike-in (30/300/3000 celulas)
# a lo largo de las muestras (gates) y replicas, relativizado al Spike-In 2.
#
# Requiere: data/processed/data_long.tsv (generado por
# R/01_activity_stats/build_activity_stats.R). Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "spike_in_robustez"
out_dir <- fig_dir(slug)
clrs <- ghibli::ghibli_palette(name = "PonyoMedium")

data_long <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE)

spike_id <- c("spikeIn_SV40", "spikeIn_CMVe", "spikeIn_CMVeMut")
spikes <- data_long %>%
  filter(seq_id %in% spike_id) %>%
  mutate(spike = ifelse(seq_id == "spikeIn_CMVeMut", "Spike-In 1",
    ifelse(seq_id == "spikeIn_SV40", "Spike-In 2", "Spike-In 3")
  ))

p <- spikes %>%
  filter(spike == "Spike-In 2") %>%
  pivot_wider(names_from = "spike", values_from = "counts") %>%
  select(`Spike-In 2`, sample, rep) %>%
  left_join(spikes, ., by = c("sample", "rep")) %>%
  mutate(rel_spike = counts / `Spike-In 2`) %>%
  ggplot(aes(as_factor(as.numeric(sample)), rel_spike, col = spike, group = spike)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  scale_y_log10(labels = scales::label_number(), breaks = 10^(-2:2)) +
  scale_color_manual(values = clrs[c(3, 4, 6)]) +
  facet_wrap(~rep) +
  ylab("Counts relativizados (log)") +
  xlab("Muestra") +
  labs(col = "Spike-In") +
  theme_pubclean() +
  theme(text = element_text(size = 15))

ggsave(file.path(out_dir, paste0(slug, ".jpg")), p, width = 9, height = 6.75, units = "in")

message("Figura guardada en ", out_dir)
