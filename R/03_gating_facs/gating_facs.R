# gating_facs (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Estrategia de gating por fluorescencia EGFP.
# Panel A: definicion del umbral inferior (densidad de fluorescencia EGFP
#   de la library vs. el control no transfectado).
# Panel B: discretizacion final en 6 gates (mas la fraccion sin clasificar).
#
# Procesa los .fcs crudos (excluye el archivo "library 1 sort.fcs", de
# 264MB, que no hace falta para estos dos paneles) - no usa ningun
# intermedio de transcriptional_library/Analysis/.
#
# Requiere: data/raw/Citometry_18-09-2023/{library 1,non transfected
# control}.fcs. Run from the TesisDoc repo root.

library(tidyverse)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
source("R/functions/fig_paths.R")

# flowCore/BiocGenerics define S4 generics that shadow dplyr's
# filter/select; pin them to the dplyr versions (see also
# R/00_prom_features/build_prom_features.R for the same issue).
filter <- dplyr::filter
select <- dplyr::select

slug <- "gating_facs"
out_dir <- fig_dir(slug)
data_path <- "data/raw/Citometry_18-09-2023"

fs <- read.flowSet(
  files = list.files(data_path, pattern = "fcs$", full.names = TRUE),
  transformation = FALSE, emptyValue = FALSE, alter.names = TRUE, truncate_max_range = FALSE
)
sampleNames(fs) <- str_remove(sampleNames(fs), ".fcs$")

# --- Gating: outliers, luego singlets -------------------------------------

outlier_gate <- rectangleGate(filterId = "-outlier", "FSC1.Height" = c(.4e9, 3.5e9), "SSC1.Height" = c(5e7, 1e9))
gs <- GatingSet(fs)
gs_pop_add(gs, outlier_gate, parent = "root")
recompute(gs)

singlet_gate <- rectangleGate(filterId = "singlet", "FSC1.Height" = c(0, 1e10), "FSC1.Width" = c(0, 2e9))
gs_pop_add(gs, singlet_gate, parent = "-outlier")
recompute(gs)

data <- gs_pop_get_data(gs, "singlet")
data_df <- map2(seq_along(data), sampleNames(data), ~ exprs(data[[.x]]) %>%
  as_tibble() %>%
  select(FL18.Height) %>%
  mutate(sample = .y)) %>% list_rbind()

# --- Panel A: umbral inferior (Library vs. control no transfectado) -------

panel_a <- data_df %>%
  mutate(sample = factor(ifelse(sample == "library 1", "Library", "WT"), levels = c("WT", "Library"))) %>%
  ggplot(aes(FL18.Height, fill = sample)) +
  scale_x_log10(lim = c(1E4, 1E7)) +
  geom_density(alpha = 0.8) +
  xlab("log EGFP") +
  scale_fill_manual(values = c(Library = "#97E196", WT = "#f5fdff")) +
  theme_bw() +
  theme(
    legend.key.size = unit(4, "line"), legend.text = element_text(size = 35),
    legend.title = element_text(size = 40), axis.title = element_text(size = 45),
    axis.text = element_text(size = 35), title = element_text(size = 55)
  )
ggsave(file.path(out_dir, "panel_a_umbral_inferior.jpg"), panel_a, width = 12, height = 9, units = "in")

# --- Panel B: discretizacion en 6 gates -----------------------------------
# breaks: limites de gate calibrados en el software del sorter (EGFP.plo
# en transcriptional_library/Experimental_data/Citometry/18-09-2023), no
# se derivan de estos datos.

breaks <- c(5.56, 8.35, 12.53, 17.98, 27, 40.54, 60.86)
df <- data_df %>% filter(sample == "library 1", FL18.Height > 0)

dens <- density(df$FL18.Height, n = 100000)
dens <- data.frame(x = dens$x, y = dens$y) %>% filter(x > 0)
dens <- dens %>%
  rowwise() %>%
  mutate(bin = findInterval(x, breaks * 15000) - 1)
dens <- do.call("rbind", lapply(split(dens, dens$bin), function(d) {
  d <- rbind(d[1, ], d, d[nrow(d), ])
  d$y[c(1, nrow(d))] <- 0
  d
}))
dens <- dens %>%
  ungroup() %>%
  mutate(Gate = factor(ifelse(bin == -1, "Sin\nclasificar", bin)))
gate_colors <- setNames(
  c("grey75", rcartocolor::carto_pal(7, "Emrld")),
  c("Sin\nclasificar", as.character(0:6))
)

panel_b <- dens %>%
  ggplot(aes(x, y)) +
  geom_polygon(aes(fill = Gate)) +
  scale_fill_manual(values = gate_colors, name = "Fracción") +
  xlab("log10 EGFP") +
  ylab("Densidad") +
  scale_x_log10() +
  theme_bw() +
  theme(legend.key.size = unit(1, "line"), text = element_text(size = 20))
ggsave(file.path(out_dir, "panel_b_discretizacion_gates.jpg"), panel_b, width = 9, height = 6.75, units = "in")

message("Paneles A y B guardados en ", out_dir)
