# validacion_rtqpcr_egfp (ver docs/mapping_figuras.csv para el numero de
# figura vigente)
# RT-qPCR de ARNm de EGFP por gate y replica, relativizado al gate 1.
#
# Requiere: data/raw/qPCR_Gates.xlsx. Run from the TesisDoc repo root.
# NOTA: los nombres de replica del archivo fuente (rep_1/2/3) no
# coinciden con el orden final ("Rep 3"/"Rep 1"/"Rep 2" respectivamente)
# - reasignacion intencional heredada de transcriptional_library/Tesis/tesis.R.

library(tidyverse)
library(ggpubr)
source("R/functions/fig_paths.R")

slug <- "validacion_rtqpcr_egfp"
out_dir <- fig_dir(slug)

qpcr <- readxl::read_xlsx("data/raw/qPCR_Gates.xlsx") %>%
  filter(Task == "Unknown", Sample != "NTC") %>%
  mutate(
    Cq = as.numeric(Cq),
    Gate = str_remove(Sample, "_rep_.") %>% str_remove("Gate_") %>% as.numeric() %>% as_factor(),
    Rep = str_remove(Sample, "Gate_._")
  )
qpcr$tech <- rep(c(1, 2, 3), times = 7 * 3 * 2)

qpcr_wide <- qpcr %>%
  select(Cq, Target, Sample, Rep, tech, Gate) %>%
  pivot_wider(values_from = Cq, names_from = Target) %>%
  mutate(delta = EGFP - GAPDH)

df <- qpcr_wide %>%
  group_by(Rep, Gate) %>%
  mutate(mean_GAPDH = mean(GAPDH), delta = 2^-(EGFP - mean_GAPDH))
mean_g2 <- df %>%
  filter(Gate == 1) %>%
  select(Rep, delta) %>%
  group_by(Rep) %>%
  summarise(mean_g2 = mean(delta))
df <- df %>%
  left_join(mean_g2, by = "Rep") %>%
  mutate(FC = delta / mean_g2, Gate = as.factor(as.numeric(Gate) - 1))

qpcr_panel <- function(data, rep_id, title, show_ylab) {
  data %>%
    filter(Rep == rep_id) %>%
    ggplot(aes(Gate, FC, fill = Gate, col = Gate)) +
    geom_line(linewidth = 2) +
    geom_point(size = 5) +
    geom_point(shape = 1, size = 5, colour = "black") +
    ggtitle(title) +
    scale_y_log10() +
    ylab(if (show_ylab) "ARN EGFP (Veces de cambio respecto a la Población 0)" else "") +
    rcartocolor::scale_color_carto_d(palette = "Emrld") +
    theme_pubr() +
    theme(legend.position = "none")
}

p1 <- qpcr_panel(df, "rep_1", "Rep 3", FALSE)
p2 <- qpcr_panel(df, "rep_2", "Rep 1", TRUE)
p3 <- qpcr_panel(df, "rep_3", "Rep 2", FALSE)

combined <- cowplot::plot_grid(p2, p3, p1, nrow = 1)
ggsave(file.path(out_dir, paste0(slug, ".jpg")), combined, width = 12, height = 6.75, units = "in")

message("Figura guardada en ", out_dir)
