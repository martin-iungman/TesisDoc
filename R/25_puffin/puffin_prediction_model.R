# PUFFIN prediction pipeline sobre la library propia - REFERENCIA, NO SE
# EJECUTA.
#
# Portado (sin correr) de transcriptional_library/Analysis/scripts/
# puffin_processing.qmd. Corre el modelo de secuencia "Puffin" (Dudnyk
# et al. 2024) via puffin.py (herramienta externa Python/TensorFlow +
# pesos del modelo, no incluidos en este repo) sobre cada promotor de la
# library para predecir un score de expresion en la posicion del TSS.
# Necesita ~24000 secuencias x 5 replicas de relleno aleatorio upstream
# (~120000 predicciones) y se corrio una unica vez, fuera de R, generando
# el archivo que este script escribe al final via write_tsv().
#
# Los CSV de prediccion por-secuencia que este pipeline escribe
# (puffin/pred/*.csv) ya no existen en disco - solo sobrevivio el
# resumen agregado (transcriptional_library/Analysis/Tables/
# Dudnyk_puffin_prediction_summ.tsv), referenciado como path_puffin_pred
# en R/00_prom_features/analysis_tables_exceptions.R.
# R/25_puffin/puffin.R usa ese resumen directamente y NO llama a este
# script.
#
# Se mantiene aca para poder rehacer la prediccion desde las secuencias
# crudas si puffin.py y los pesos del modelo vuelven a estar disponibles.

library(tidyverse)
library(furrr)

# --- 1. Armar un FASTA por promotor de la library, con n_random_seq
# variantes independientes de relleno aleatorio upstream (se necesitan
# 325pb de cada lado para ser independiente del contexto genomico de la
# insercion; el relleno downstream es una secuencia real conocida:
# primers/LoxP/parte del CDS de dsRed-EGFP) --------------------------------

n_random_seq <- 5
up <- "gtacatcaagtgtatcactagcgtttaaacttaagcttccatggattacaaggatgacgatgacaagggggtacctgccccaaaaaaaaaacgcaaagtggaggacccagtaccaggatctagaggtaggtgatcctcctgctgctttggttcagggttttgcttgaggggggggggtggtgatttccttgccatgggcagactgagcagaaaaggccattgggaccatgttctgaatgcctccacctcaaccaccggccggtaggaccaaagccaccccgtgttttctcaggatctcttttcccagggagatccctcggcccaa"
down <- tolower("ATAACTTCGTATAGCATACATTATACGAAGTTATATGGATCCATACTAGgacattgattattg")

fasta <- readLines("Library_data/data/promoters_wo_dupl.fa")
ids <- fasta[grep("^>", fasta)]
seqs <- fasta[-grep("^>", fasta)]

for (i in 1:n_random_seq) {
  random <- sample(c("A", "C", "G", "T"), size = 325 - nchar(down), replace = TRUE) %>% paste0(collapse = "")
  for (j in seq_along(seqs)) {
    id_ <- paste0(ids[j], "_random", i)
    cat(id_, file = "puffin/library_puffin.fa", append = TRUE, sep = "\n")
    cat(paste0(random, down, seqs[j], up), file = "puffin/library_puffin.fa", append = TRUE, sep = "\n")
  }
}

# --- 2. Correr el modelo (herramienta Python externa, modificada para
# escribir toda la salida en puffin/pred/) ----------------------------------

# cd puffin && mkdir pred pred/processed
# python puffin.py sequence puffin/library_puffin.fa
system("cd puffin && mkdir -p pred pred/processed && python puffin.py sequence puffin/library_puffin.fa")

# --- 3. Resumir las 5 "replicas" de relleno aleatorio por promotor
# (mediana por posicion) y calcular estadisticos por promotor: score de
# prediccion en la posicion del TSS de EPD (posicion 16), score medio en
# la ventana +-10pb del TSS, y posicion/valor del score maximo en toda
# la secuencia. Procesado en chunks (FILE_BINS ids a la vez) y en
# paralelo - son ~24000*5 archivos. -----------------------------------------

process_sample <- function(id) {
  if (is.na(id)) {
    return(NULL)
  }
  id <- str_replace(id, "\\+", "_") %>% str_replace(">", "_")
  files <- paste0("puffin/pred/", id, "_random", 1:5, ".csv")
  obj <- map(files, ~ read_csv(.x, show_col_types = FALSE) %>% t() %>% as_tibble()) %>% keep(~ nrow(.x) > 0)
  if (length(obj) != 5) message("WARNING: ", id, " with only ", length(obj), " elements")
  for (i in seq_along(obj)) {
    names(obj[[i]]) <- obj[[i]][1, ]
    obj[[i]] <- obj[[i]][-1, ] %>% mutate(across(-Sequence, as.numeric))
  }
  obj <- obj %>%
    list_rbind() %>%
    group_by(Coordinate) %>%
    summarise(across(where(is.numeric), median), sd_prediction = sd(Prediction), Sequence = unique(Sequence), rep_random = length(obj)) %>%
    ungroup() %>%
    mutate(seq_id = id)
  write_csv(obj, paste0("puffin/pred/processed/", id, ".csv"))
  obj
}

df_pred <- data.frame()
FILE_BINS <- 200
plan(multisession, workers = 10)
for (i in 0:(length(ids) / FILE_BINS)) {
  message("starting ID ", i * FILE_BINS)
  m <- future_map(ids[(i * FILE_BINS + 1):((i + 1) * FILE_BINS)] %>% str_remove("^>"), process_sample)
  g1 <- m %>%
    list_rbind() %>%
    group_by(seq_id) %>%
    filter(Prediction == max(Prediction)) %>%
    rename(max_score_pred_puffin = Prediction, pos_max_score_pred_puffin = Coordinate) %>%
    select(seq_id, max_score_pred_puffin, pos_max_score_pred_puffin)
  g2 <- m %>%
    list_rbind() %>%
    filter(Coordinate == 252 - 16) %>%
    rename(
      teoTSS_score = Prediction, teoTSS_inr_eff = `Sum of initiator effect`,
      teoTSS_motif_eff = `Sum of motif effect`, teoTSS_trinucl_eff = `Sum of trinucleotide effect`
    ) %>%
    select(teoTSS_score, teoTSS_inr_eff, teoTSS_motif_eff, teoTSS_trinucl_eff, seq_id)
  g3 <- m %>%
    list_rbind() %>%
    filter(Coordinate > 252 - 26, Coordinate < 252 - 6) %>%
    group_by(seq_id) %>%
    summarise(TSS20_mean_score_pred_puffin = mean(Prediction))
  df_pred <- rbind(df_pred, full_join(g2, g1) %>% full_join(g3))
}

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(df_pred, "data/processed/dudnyk_puffin_prediction.tsv")
