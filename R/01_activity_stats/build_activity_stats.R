# Builds the shared per-promoter activity statistics (mean/median/mode/
# variance of expression across sorting gates, from the amplicon-seq
# gate counts) consumed by the figure scripts for M8
# (correccion_sesgo_muestreo), R1 (representatividad_library_venn), R2
# (sesgo_representatividad), R3 (spike_in_robustez), R4
# (distribuciones_expresion), R5 (filtros_consistencia_replicas) and R6
# (tata_cgi_actividad).
# Run from the TesisDoc repo root. Output: data/processed/data_long.tsv,
# data/processed/activity_stats_full.tsv, data/processed/activity_stats.tsv
# (post sum_counts>=100 filter only - NOT the bimodal/consistency filters,
# those are figure-specific and applied in R/05_correccion_sesgo and
# R/13_filtros).

library(tidyverse)
source("R/functions/perc_cells.R")

# --- Load and normalize gate counts ---------------------------------------

files <- list.files("data/raw/Count_table", "counts_ampliconseq2023_.-..tsv", full.names = TRUE)
data <- map(files, ~ read_tsv(.x, col_names = c("seq_id", "counts"), show_col_types = FALSE))
names(data) <- files %>% str_remove("^.+2023_") %>% str_remove(".tsv$")
data <- imap(data, ~ mutate(.x, counts = counts / 2, sample_rep = .y))

data_long <- data %>%
  list_rbind() %>%
  mutate(
    sample = factor(as.numeric(str_remove(sample_rep, "-.$")) - 1),
    rep = factor(str_replace(sample_rep, "^.-", "Rep ") %>% str_replace("3", "2")),
    name = str_remove(seq_id, "^FP.{6}_")
  )

# library-size normalization (per replicate) then sampling-effort correction
counts_tot <- data_long %>%
  group_by(sample, rep) %>%
  summarise(sumCounts = sum(counts), .groups = "drop") %>%
  group_by(rep) %>%
  mutate(mean_size_rep = mean(sumCounts)) %>%
  ungroup()
data_long <- data_long %>%
  left_join(counts_tot, by = c("sample", "rep")) %>%
  mutate(counts_norm_libsize = counts / (sumCounts / mean_size_rep)) %>%
  select(-mean_size_rep)
data_long <- data_long %>%
  left_join(perc_cells, by = "sample") %>%
  mutate(counts_norm = counts_norm_libsize * rel_factor) %>%
  select(-c(sumCounts, perc, sample_rep, rel_factor))

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_tsv(data_long, "data/processed/data_long.tsv")
message("Wrote data/processed/data_long.tsv with ", nrow(data_long), " rows")

# --- Per-promoter activity stats across gates (weighted by counts) -------
# Treats each normalized count as one "observation" of its gate (1-6) and
# computes the mean/median/mode/variance of that distribution per
# seq_id x replicate - i.e. where in the sorted population a promoter's
# activity is centered, and how spread out (noisy) it is.

get_stats <- function(df, col_counts) {
  if (!all(df[[col_counts]] == 0)) {
    k <- sum(df[[col_counts]] != 0)
    df$sample <- as.numeric(df$sample)
    vctr <- map(seq_along(unique(df$sample)), ~ rep(df$sample[.x], times = df[[col_counts]][.x])) %>% list_c()
    return(data.frame(
      seq_id = unique(df$seq_id), name = unique(df$name), positive_gates = k,
      mean = mean(vctr), median = median(vctr),
      mode = df$sample[which.max(df[[col_counts]])], var = var(vctr)
    ))
  }
}

t0 <- Sys.time()
stats <- map(unique(data_long$rep), ~ data_long %>%
  filter(rep == .x) %>%
  group_split(seq_id) %>%
  map(~ get_stats(.x, "counts_norm")) %>%
  list_rbind() %>%
  mutate(rep = .x)) %>%
  list_rbind()
message("get_stats took ", round(as.numeric(Sys.time() - t0, units = "secs")), "s")

stats <- data_long %>%
  group_by(seq_id, rep) %>%
  summarise(
    sum_counts = sum(counts), sum_norm_counts_lib = sum(counts_norm_libsize),
    sum_norm_counts = sum(counts_norm), .groups = "drop"
  ) %>%
  left_join(stats, ., by = c("seq_id", "rep"))

write_tsv(stats, "data/processed/activity_stats_full.tsv")
message("Wrote data/processed/activity_stats_full.tsv with ", nrow(stats), " rows")

# minimum-counts filter
stats_filt <- stats %>% filter(sum_counts >= 100)
write_tsv(stats_filt, "data/processed/activity_stats.tsv")
message("Wrote data/processed/activity_stats.tsv with ", nrow(stats_filt), " rows")

# --- Staged QC filters (shared by M8, R1, R2, R4, R5, R6) -----------------
# consistent: drop promoters whose median gate differs by >1 between reps
# (this is what R5/filtros_consistencia_replicas actually plots).
# highconf: consistent AND not bimodal (var > 1.2x the max variance
# expected for that many positive gates) - the final high-confidence set.

inconsistent_median <- stats_filt %>%
  select(seq_id, median, rep) %>%
  pivot_wider(names_from = rep, values_from = median) %>%
  mutate(median_dif = abs(as.numeric(`Rep 1`) - as.numeric(`Rep 2`))) %>%
  filter(median_dif > 1) %>%
  pull(seq_id)
stats_consistent <- stats_filt %>% filter(!seq_id %in% inconsistent_median)

max_var <- function(k) {
  if (k %% 2 == 1) return((sum(seq(0, (k - 1) / 2)^2) * 2) / k)
  if (k %% 2 == 0) return((sum((seq(1, k / 2) - 0.5)^2) * 2) / k)
}
stats_consistent <- stats_consistent %>%
  mutate(max_var = sapply(positive_gates, max_var), bimodal = var > max_var * 1.2)
write_tsv(stats_consistent, "data/processed/activity_stats_consistent.tsv")
message("Wrote data/processed/activity_stats_consistent.tsv with ", nrow(stats_consistent), " rows")

stats_highconf <- stats_consistent %>% filter(!bimodal)
write_tsv(stats_highconf, "data/processed/activity_stats_highconf.tsv")
message("Wrote data/processed/activity_stats_highconf.tsv with ", nrow(stats_highconf), " rows")
