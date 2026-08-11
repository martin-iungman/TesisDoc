# Shared intermediate: clasificacion de promotores alternativos, usada
# por toda la seccion Resultados 5 (R5.x). Portado de
# transcriptional_library/Analysis/scripts/final_github.R (seccion
# "### Alternative promoters", lineas 1295-1521) - variantes mas
# elaboradas del mismo analisis viven en alt_prom.qmd/alt_prom.R, pero
# final_github.R es la version que efectivamente genera los plots del
# paper.
#
# Clasificacion primaria (por gen): "non_detected" (ningun promotor del
# gen supera 1 TPM en ninguna muestra), "unique" (perc_tpm_gene>0.999:
# un solo promotor concentra >99.9% de la actividad del gen - esto
# incluye, pero no se limita a, genes con un solo promotor anotado en
# EPD; ver nota mas abajo sobre por que epd_unique del original es
# codigo muerto), "Main promoter" (cumple las dos condiciones: mayor
# %tpm del gen Y mayor cantidad de muestras donde es el mas activo) o
# "unclassified" (cumple solo una de esas dos condiciones - xor).
#
# Clasificacion secundaria (por par Main-Secondary, dentro de genes con
# multiples promotores activos): "switch" (el secundario supera al
# principal en al menos una muestra: tpm>5 y, o bien ratio_alt_main>1.5
# en escala log, o bien el principal esta practicamente inactivo,
# main_tpm<1), "correlated" (sin switch, pero correlacion de Pearson o
# Spearman entre ambos >0.5 a lo largo de las muestras) o "independent"
# (ninguna de las anteriores).
#
# NOTA epd_unique: el original tambien calcula epd_unique (genes con un
# solo promotor ANOTADO en EPD, sin importar su actividad) pero nunca lo
# usa - es codigo muerto. Se puede demostrar que es redundante: si un
# gen tiene un solo promotor en sample_activity_prom, tpm_gene (la suma
# agrupada por gen) es necesariamente igual a su propio tpm, entonces
# perc_tpm_gene=tpm/tpm=1 siempre (los "no detectados" ya se excluyen
# aparte) - es decir epd_unique (detectados) es subconjunto de
# unique_prom por construccion. No se porta.
#
# CORRECCION 2026-08-05: N_switch=sum(switch) en el original NO usa
# na.rm=TRUE. switch es NA cuando tpm=main_tpm=0 en la misma muestra
# (log10(0)/log10(0)=NaN), lo cual puede pasar dentro del subconjunto de
# muestras usado aca (solo se filtro tpm_gene_sample>1 a nivel del GEN,
# no de cada promotor individual). Sin na.rm, un solo NA en cualquier
# muestra volvia NA todo N_switch de ese promotor, que en el ifelse
# posterior (prom_alt=ifelse(N_switch>0,...)) le asignaba prom_alt=NA y
# lo hacia desaparecer silenciosamente de todos los plots downstream -
# aun si en el resto de sus muestras mostraba evidencia clara de switch/
# correlated/independent. Corregido con na.rm=TRUE (una muestra
# ambigua 0-vs-0 se trata como "no es evidencia de switch", no como dato
# faltante para todo el promotor).
#
# Requiere: data/external/EPD/human38_epdnew.bed, data/raw/library.bed y
# la excepcion sample_CAGE_activity.tsv (ver
# R/00_prom_features/analysis_tables_exceptions.R). Run from the
# TesisDoc repo root. Output:
# data/processed/prom_alt_classification.tsv (name, prom_alt - primaria
#   y secundaria combinadas en una sola columna, como promalt_df del
#   original)
# data/processed/prom_alt_pairs.tsv (pares Secondary-Main: name,
#   main_name, sample, tpm, main_tpm, cor_pearson, cor_spearman,
#   ratio_alt_main, switch, N_switch, prom_alt - como wide_df)

library(tidyverse)
library(rtracklayer)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

source("R/00_prom_features/analysis_tables_exceptions.R")

# --- Cargar datos crudos --------------------------------------------------

epd <- import.bed("data/external/EPD/human38_epdnew.bed") %>%
  IRanges::promoters(downstream = 16, upstream = 236)
lib <- import.bed("data/raw/library.bed")
lib$name <- str_remove(lib$name, "^FP.{6}_")

sample_activity <- read_tsv(path_sample_cage_activity, show_col_types = FALSE) %>%
  mutate(tpm = counts * 1e6 / libsize)
sample_activity_prom <- sample_activity %>%
  filter(name %in% epd$name) %>%
  mutate(gene_sym = str_remove(name, "_.{1,2}$"))

# --- Clasificacion primaria ------------------------------------------------

detected_names <- sample_activity_prom %>% filter(tpm > 1) %>% pull(name) %>% unique()
non_detected <- sample_activity_prom %>% pull(name) %>% unique() %>% setdiff(detected_names)

tot_counts <- sample_activity_prom %>%
  group_by(name, gene_sym) %>%
  summarise(tpm = sum(tpm), .groups = "drop") %>%
  group_by(gene_sym) %>%
  mutate(tpm_gene = sum(tpm), perc_tpm_gene = tpm / tpm_gene) %>%
  ungroup()

unique_prom <- tot_counts %>%
  filter(perc_tpm_gene > 0.999, !name %in% non_detected) %>%
  pull(name) %>%
  unique()

sample_activity_alt <- sample_activity_prom %>%
  filter(!name %in% c(non_detected, unique_prom)) %>%
  group_by(gene_sym, sample) %>%
  mutate(tpm_gene_sample = sum(tpm), is_highest_in_sample = tpm == max(tpm)) %>%
  ungroup()

highest_counts <- sample_activity_alt %>%
  filter(tpm_gene_sample > 1) %>%
  group_by(name) %>%
  summarise(highest_in_N_samples = sum(is_highest_in_sample), .groups = "drop")

sample_activity_alt <- highest_counts %>%
  inner_join(sample_activity_alt, by = "name") %>%
  group_by(gene_sym, sample) %>%
  mutate(main_N_samples = highest_in_N_samples / sum(highest_in_N_samples)) %>%
  ungroup()

main_prom <- sample_activity_alt %>%
  left_join(tot_counts %>% select(-tpm), by = c("name", "gene_sym")) %>%
  group_by(gene_sym) %>%
  filter((perc_tpm_gene == max(perc_tpm_gene)) & (main_N_samples == max(main_N_samples))) %>%
  ungroup()

unclassified <- sample_activity_alt %>%
  left_join(tot_counts %>% select(-tpm), by = c("name", "gene_sym")) %>%
  group_by(gene_sym) %>%
  filter(xor(perc_tpm_gene == max(perc_tpm_gene), main_N_samples == max(main_N_samples))) %>%
  ungroup()

# --- Clasificacion secundaria (pares Secondary-Main) -----------------------

# dplyr detecta un join many-to-many aca (multiples filas de un lado y
# del otro comparten la misma clave sample+gene_sym dentro de un gen con
# varios secundarios) - verificado que no genera duplicados espurios: 0
# de los 9628 promotores secundarios terminan con mas de un main_name
# distinto tras el join + filter(name != main_name).
wide_df <- sample_activity_alt %>%
  select(name, gene_sym, tpm, sample) %>%
  inner_join(
    main_prom %>% select(name, gene_sym, tpm, sample) %>% unique() %>%
      rename(main_name = name, main_tpm = tpm),
    by = c("sample", "gene_sym"), relationship = "many-to-many"
  ) %>%
  filter(name != main_name)

cor_data_pearson <- wide_df %>%
  ungroup() %>%
  filter(tpm + main_tpm > 0) %>%
  group_split(name) %>%
  map(~ summarise(.x,
    cor = cor(tpm, main_tpm, method = "pearson"),
    n_samples = n(), name = unique(name), main_name = unique(main_name)
  )) %>%
  list_rbind()

cor_data_spearman <- wide_df %>%
  ungroup() %>%
  filter(tpm + main_tpm > 0) %>%
  group_split(name) %>%
  map(~ summarise(.x,
    cor = cor(tpm, main_tpm, method = "spearman"),
    n_samples = n(), name = unique(name), main_name = unique(main_name)
  )) %>%
  list_rbind()

cor_data <- left_join(cor_data_pearson, cor_data_spearman,
  by = c("n_samples", "name", "main_name"),
  suffix = c("_pearson", "_spearman")
)

wide_df <- wide_df %>%
  left_join(cor_data, by = c("name", "main_name")) %>%
  mutate(
    ratio_alt_main = log10(tpm) / log10(main_tpm),
    switch = (ratio_alt_main > 1.5) & (tpm > 5),
    switch = ifelse(tpm > 5 & main_tpm < 1, TRUE, switch)
  ) %>%
  group_by(name) %>%
  mutate(N_switch = sum(switch, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(N_switch))

wide_df <- wide_df %>%
  mutate(prom_alt = ifelse(N_switch > 0, "switch",
    ifelse(cor_pearson > 0.5 | cor_spearman > 0.5, "correlated", "independent")
  ))

# --- Combinar clasificacion primaria + secundaria por promotor -------------

promalt_df <- bind_rows(
  tibble(name = unique_prom, prom_alt = "unique"),
  tibble(name = non_detected, prom_alt = "non_detected"),
  tibble(name = unique(main_prom$name), prom_alt = "Main promoter"),
  wide_df %>% select(name, prom_alt) %>% unique(),
  tibble(name = unique(unclassified$name), prom_alt = "unclassified")
)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(promalt_df, "data/processed/prom_alt_classification.tsv")
message("Wrote data/processed/prom_alt_classification.tsv with ", nrow(promalt_df), " rows")

pairs_df <- wide_df %>%
  select(name, main_name, sample, tpm, main_tpm, cor_pearson, cor_spearman, ratio_alt_main, switch, N_switch, prom_alt)
write_tsv(pairs_df, "data/processed/prom_alt_pairs.tsv")
message("Wrote data/processed/prom_alt_pairs.tsv with ", nrow(pairs_df), " rows")
