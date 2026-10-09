# ruido_summary (ver docs/mapping_figuras.csv para el numero de figura
# vigente)
# Resumen del efecto de cada feature booleana del promotor sobre el
# ruido transcripcional: AUC de la feature prediciendo var_rank_sw alto
# (rango de varianza dentro de cada bin de actividad media - ver
# R/19_ruido_cgi_tata), un solo panel combinado (sin separar seq/endo,
# a diferencia de R7).
#
# Significancia: p-valor de Wald usando la varianza de DeLong del AUC
# (z = (AUC-0.5)/SE_delong), corregido por multiples tests con BH -
# resync 2026-09-30, reemplaza al criterio anterior de "IC bootstrap que
# no cruza 0.5". Se verifico que el metodo anterior en realidad ya
# calculaba el IC via DeLong (pROC::ci.auc usa DeLong por default; el
# argumento boot.n se ignoraba silenciosamente sin method="bootstrap"),
# asi que el cambio no es metodologico sino que agrega el p-valor
# explicito para poder corregir por tests multiples, cosa que un IC no
# permite. Se comparo explicitamente contra bootstrap real (boot.n=2000,
# resampleo estratificado) para las 33 features x 2 replicas: acuerdo
# practicamente perfecto (Spearman en -log10(p) = 0.99), con la unica
# discrepancia de significancia post-BH en un caso limite (p~0.05), y el
# bootstrap con piso de resolucion en 1/boot.n=0.0005 (da p=0 exacto en
# 25/66 tests donde el efecto es muy fuerte, cosa que DeLong no tiene).
#
# Requiere: data/processed/activity_stats_highconf.tsv y
# data/processed/prom_df.tsv (ver R/00_prom_features y
# R/01_activity_stats). Corre en segundos (DeLong es una formula cerrada,
# no hace falta resamplear). Run from the TesisDoc repo root.
#
# Ported from transcriptional_library/Analysis/scripts/final_github.R
# (bloque de resumen de la seccion de ruido, plot_auc2). build_tidy_
# features()/feature_groups (identicas a las de R/15_summary_features,
# el resumen de actividad) viven en R/functions/plot_helpers.R. Incluye
# "Promotor unidireccional" (ver orientacion_promotor) como cualquier otra feature de
# build_tidy_features(), sin tratamiento especial.

library(tidyverse)
library(ggpubr)
library(fastDummies)
library(pROC)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

slug <- "ruido_summary"
out_dir <- fig_dir(slug)

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name"))

data <- add_noise_rank(data)

tidy_data <- build_tidy_features(data, keep = "var_rank_sw")

# --- AUC (p-valor DeLong, corregido por BH) de cada feature prediciendo --
# --- ruido alto ------------------------------------------------------------

vbles <- tidy_data %>% select(-rep, -var_rank_sw) %>% names()
repname <- unique(tidy_data$rep)

auc_df <- map(repname, function(r) {
  df <- tidy_data %>% filter(rep == r)
  map(vbles, function(f) {
    resp <- df[[f]]
    ok <- !is.na(resp)
    roc_obj <- pROC::roc(resp[ok], df[["var_rank_sw"]][ok], ci = TRUE, direction = ">", quiet = TRUE)
    auc_val <- as.numeric(roc_obj$ci[2])
    se_delong <- sqrt(pROC::var(roc_obj, method = "delong"))
    pval <- 2 * pnorm(-abs((auc_val - 0.5) / se_delong))
    tibble(feature = f, rep = r, AUC = auc_val, ci2.5 = as.numeric(roc_obj$ci[1]), ci97.5 = as.numeric(roc_obj$ci[3]), pval = pval)
  }) %>% list_rbind()
}) %>% list_rbind() %>%
  mutate(pval_corr = p.adjust(pval, method = "BH"))

write_tsv(auc_df, file.path(out_dir, "ruido_summary_auc.tsv"))

# --- Grafico resumen --------------------------------------------------
# Features con pval_corr<0.05 Y mismo sentido (ruido alto/bajo) en ambas
# replicas (mismo criterio de consistencia que R7). plot_noise_auc_
# summary() en plot_helpers.R (compartida con R3.3, la version de esto
# para TFs de ReMap, que todavia usa el criterio de IC ya que no calcula
# pval_corr).

panel <- plot_noise_auc_summary(auc_df)
ggsave(file.path(out_dir, "summary_noise.jpg"), panel, width = 13.5, height = 6.75, units = "in")

message("Figuras guardadas en ", out_dir)
