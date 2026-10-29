# Exploratorio, continuacion de remap_metilacion.R (misma infraestructura):
# p7 maquinaria de pausa / P-TEFb y p8 TFIIA (GTF2A2) y TATA-box, mas el caso
# lateral de SNAPC1.
#
# Estadisticos (los de remap_metilacion.R): coeficiente de `unido` del modelo
# ajustado (var_rank_sw ~ unido + log10(hek_tpm + 1) + mean_dnase); residuo de
# un loess robusto del coeficiente sobre log10(n unidos); percentil en la
# ventana de cobertura +-25% con p por permutacion exacto y BH sobre los 982
# factores. Cada variante de modelo o estrato recalcula residuo y ventana
# contra los 982 factores bajo esa misma variante.
#
# p7 Pausa y elongacion:
#   a) Wilcoxon de grupo (pausa, preiniciacion, mediador) vs. el resto y
#      pausa vs. preiniciacion, con y sin CDK9 / CCNT2 (y preiniciacion con y
#      sin GTF2A2), sobre residuo y percentil;
#   b) distinguir P-TEFb de "estar ocupado por la maquinaria": modelo
#      + log10(n de los otros 981 factores unidos + 1) para los 982 factores;
#      modelo conjunto a nivel de promotor con los indicadores de union de
#      pausa y preiniciacion; CDK9 dentro de los promotores unidos por
#      TAF1, TBP y GTF2B;
#   c) estratos por isla CpG: se leen de p4_residuo_ventana_por_factor_cgi.tsv;
#   d) procedencia (lineas de ReMap, experimentos de ChIP-Atlas) y efecto por
#      linea celular de CDK9 y CCNT2.
# p8 GTF2A2:
#   a) enriquecimiento en TATA-box canonica (TATA_EPD; identica a la clase
#      "Canonica" de tata_scan.R, ver tata_ruido_especificidad) de los
#      promotores unidos, odds ratio por replica, ubicado entre los 982;
#   b) modelo + TATA_EPD (los 982, con ventana) y efecto dentro de promotores
#      sin / con TATA, para GTF2A2 y los controles TAF1, TBP, GTF2B;
#   c) estratos de actividad (alta = percentil 75 de mean por replica, como
#      ROC_TATA_expr), con y sin TATA como covariable, e interaccion
#      unido x actividad alta;
#   d) procedencia de GTF2A2.
#   SNAPC1: genes de snRNA en la library (GENETYPE de org.Hs.eg.db) y genes de
#   snRNA a menos de 2 kb de los promotores unidos (CHRLOC de org.Hs.eg.db).
#
# Salida en figures/exploratorio_remap_metilacion/ con prefijos p7 y p8.
# Requiere lo mismo que remap_metilacion.R mas sus salidas
# p1_residuo_ventana_por_factor.tsv y p4_residuo_ventana_por_factor_cgi.tsv,
# y data/raw/library.bed. Tarda unos minutos. Run from the TesisDoc repo root.

library(tidyverse)
library(org.Hs.eg.db)
library(ggrepel)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
first <- dplyr::first
desc <- dplyr::desc
slice <- dplyr::slice
mad <- stats::mad

out_dir <- "figures/exploratorio_remap_metilacion"
hs <- org.Hs.eg.db

cov_window <- 0.25
min_estrato <- 50

# grupos a priori. Pausa: P-TEFb, SEC, 7SK snRNP, DSIF, NELF, BRD4 (los del
# pedido mas AFF1, ELL/ELL2/ELL3, EAF1/2, MLLT1, HEXIM1/2, LARP7, MEPCE).
# Preiniciacion: TFIID (TBP y TAF de TFIID; TAF15 no: es una proteina FET, no
# una subunidad de TFIID), TFIIA, TFIIB, TFIIE, TFIIF, TFIIH, Pol II.
grupos_def <- list(
  pausa = c("CDK9", "CCNT1", "CCNT2", "NELFA", "NELFB", "NELFCD", "NELFE", "SUPT4H1", "SUPT5H", "AFF4", "BRD4", "CDK12", "CDK13", "MLLT3",
    "AFF1", "ELL", "ELL2", "ELL3", "EAF1", "EAF2", "MLLT1", "HEXIM1", "HEXIM2", "LARP7", "MEPCE"),
  preiniciacion = c("TBP", paste0("TAF", c(1:8, 9, "9B", 10:13)), "GTF2A1", "GTF2A2", "GTF2B", "GTF2E1", "GTF2E2", "GTF2F1", "GTF2F2",
    "GTF2H1", "GTF2H2", "GTF2H3", "GTF2H4", "ERCC2", "ERCC3", "CDK7", "CCNH", "MNAT1", "POLR2A"),
  mediador = c("MED", paste0("MED", c(1, 4, 6:31)))
)
grupo_lab <- c(pausa = "Pausa / elongación", preiniciacion = "Preiniciación", mediador = "Mediador", resto = "Resto")
grupo_clr <- c(pausa = "#AD343E", preiniciacion = "#1B8C8E", mediador = "#D6741F", resto = "grey70")
extremos_pausa <- c("CDK9", "CCNT2")
foco_p8 <- c("GTF2A2", "TAF1", "TBP", "GTF2B", "SNAPC1")
procedencia_tfs <- c("CDK9", "CCNT2", "GTF2A2", "SNAPC1")

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}

# --- Datos ---------------------------------------------------------------

cob <- read_tsv(file.path(out_dir, "p1_residuo_ventana_por_factor.tsv"), show_col_types = FALSE) %>% select(-grupo)
cob_cgi <- read_tsv(file.path(out_dir, "p4_residuo_ventana_por_factor_cgi.tsv"), show_col_types = FALSE) %>% select(-grupo)
tfs <- unique(cob$TF)

remap_raw <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
remap_hits <- remap_raw %>% distinct(seq_id, TF)
sets <- remap_hits %>% filter(TF %in% tfs) %>% with(split(seq_id, TF))
n_ocup <- remap_hits %>% filter(TF %in% tfs) %>% count(seq_id, name = "n_ocup")

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  group_by(rep) %>%
  mutate(act_alta = percent_rank(mean) >= 0.75) %>%
  ungroup() %>%
  mutate(log_hek = log10(hek_tpm + 1)) %>%
  left_join(n_ocup, by = "seq_id") %>%
  mutate(n_ocup = replace_na(n_ocup, 0L)) %>%
  select(seq_id, rep, gene_sym, var_rank_sw, log_hek, mean_dnase, CGI_frag, TATA_EPD, act_alta, n_ocup)
stopifnot(!anyNA(select(data, -gene_sym)))
repname <- unique(data$rep)

n_lib <- remap_hits %>%
  filter(TF %in% unlist(grupos_def)) %>%
  inner_join(distinct(data, seq_id, rep), by = "seq_id", relationship = "many-to-many") %>%
  count(TF, rep) %>%
  pivot_wider(names_from = rep, values_from = n, names_prefix = "n_")
grupos_tab <- enframe(grupos_def, "grupo", "TF") %>%
  unnest(TF) %>%
  left_join(n_lib, by = "TF") %>%
  mutate(
    incluido = TF %in% tfs,
    motivo = case_when(
      incluido ~ "incluido",
      !TF %in% remap_hits$TF ~ "sin picos de ReMap sobre promotores de la library",
      TRUE ~ "100 promotores unidos o menos en alguna réplica (umbral de remap_tf_noise_gsea)"
    )
  )
write_tsv(grupos_tab, file.path(out_dir, "p7_grupos_definicion.tsv"))
print(grupos_tab %>% filter(incluido | !is.na(`n_Rep 1`)), n = 60)
grupos <- map(grupos_def, ~ intersect(.x, tfs))
grupo_de <- function(tf) case_when(tf %in% grupos$pausa ~ "pausa", tf %in% grupos$preiniciacion ~ "preiniciacion", tf %in% grupos$mediador ~ "mediador", TRUE ~ "resto")

# --- Funciones (las de remap_metilacion.R; fit_coefs con `prep` para
# covariables que dependen del factor) ------------------------------------

fit_coefs <- function(d, fo, min_n = 0, prep = NULL) {
  map(tfs, function(tf) {
    d$unido <- as.numeric(d$seq_id %in% sets[[tf]])
    if (!is.null(prep)) d <- prep(d)
    n1 <- sum(d$unido)
    n0 <- nrow(d) - n1
    if (n1 < min_n || n0 < min_n) {
      return(tibble(TF = tf, n_unidos = n1, n_no_unidos = n0, estimate = NA_real_, se = NA_real_, pval = NA_real_))
    }
    co <- summary(lm(fo, data = d))$coefficients["unido", ]
    tibble(TF = tf, n_unidos = n1, n_no_unidos = n0, estimate = co[1], se = co[2], pval = co[4])
  }) %>% list_rbind()
}

cobertura_stats <- function(df) {
  df <- df %>% filter(!is.na(estimate))
  df$log_n <- log10(df$n_unidos)
  fit_lo <- loess(estimate ~ log_n, data = df, family = "symmetric", span = 0.5)
  df$ajuste_cobertura <- predict(fit_lo, df)
  df$residuo <- df$estimate - df$ajuste_cobertura
  win <- map(seq_len(nrow(df)), function(i) {
    ref <- df$estimate[-i][abs(df$n_unidos[-i] - df$n_unidos[i]) <= cov_window * df$n_unidos[i]]
    x <- df$estimate[i]
    n_le <- sum(ref <= x)
    n_ge <- sum(ref >= x)
    tibble(n_ventana = length(ref), percentil_ventana = n_le / length(ref), p_perm = min(1, 2 * (min(n_le, n_ge) + 1) / (length(ref) + 1)))
  }) %>% list_rbind()
  bind_cols(df, win) %>% mutate(p_perm_bh = p.adjust(p_perm, method = "BH"))
}

wilcox_tab <- function(df, var, comps) {
  imap(comps, function(cmp, nm) {
    xa <- df[[var]][df$TF %in% cmp$a & !is.na(df[[var]])]
    xb <- df[[var]][df$TF %in% cmp$b & !is.na(df[[var]])]
    w <- wilcox.test(xa, xb, conf.int = TRUE, exact = FALSE)
    tibble(
      comparacion = nm, variable = var, n_a = length(xa), n_b = length(xb),
      mediana_a = median(xa), mediana_b = median(xb),
      # rank-biserial: > 0 = el primer grupo tiene valores mas altos (mas ruido)
      rank_biserial = 2 * w$statistic / (length(xa) * length(xb)) - 1,
      HL = unname(w$estimate), HL_ci_low = w$conf.int[1], HL_ci_high = w$conf.int[2], pval = w$p.value
    )
  }) %>% list_rbind()
}

comps_p7 <- function(all_tf) {
  pausa_sin <- setdiff(grupos$pausa, extremos_pausa)
  pic_sin <- setdiff(grupos$preiniciacion, "GTF2A2")
  maq <- unlist(grupos)
  list(
    pausa_vs_resto = list(a = grupos$pausa, b = setdiff(all_tf, grupos$pausa)),
    pausa_sin_CDK9_CCNT2_vs_resto = list(a = pausa_sin, b = setdiff(all_tf, grupos$pausa)),
    preiniciacion_vs_resto = list(a = grupos$preiniciacion, b = setdiff(all_tf, grupos$preiniciacion)),
    preiniciacion_sin_GTF2A2_vs_resto = list(a = pic_sin, b = setdiff(all_tf, grupos$preiniciacion)),
    mediador_vs_resto = list(a = grupos$mediador, b = setdiff(all_tf, grupos$mediador)),
    pausa_vs_preiniciacion = list(a = grupos$pausa, b = grupos$preiniciacion),
    pausa_vs_preiniciacion_sin_GTF2A2 = list(a = grupos$pausa, b = pic_sin),
    pausa_sin_CDK9_CCNT2_vs_preiniciacion_sin_GTF2A2 = list(a = pausa_sin, b = pic_sin),
    maquinaria_vs_resto = list(a = maq, b = setdiff(all_tf, maq))
  )
}
run_tests <- function(cobdf, by = "rep") {
  cobdf %>%
    group_split(across(all_of(by))) %>%
    map(function(df) {
      map(c("residuo", "percentil_ventana"), ~ wilcox_tab(df, .x, comps_p7(df$TF))) %>%
        list_rbind() %>%
        bind_cols(distinct(df[1, by]))
    }) %>%
    list_rbind() %>%
    relocate(all_of(by))
}

# variantes de modelo para los 982 factores; estratos: solo factores
# evaluables (>= min_estrato unidos y no unidos) en ambas replicas
fo_aj <- var_rank_sw ~ unido + log_hek + mean_dnase
fo_tata <- update(fo_aj, ~ . + TATA_EPD)
fo_ocup <- update(fo_aj, ~ . + ocup_otros)
prep_ocup <- function(d) mutate(d, ocup_otros = log10(n_ocup - unido + 1))
variantes <- list(
  ajustado_mas_ocupacion = list(sub = identity, fo = fo_ocup, min_n = 0, prep = prep_ocup),
  ajustado_mas_tata = list(sub = identity, fo = fo_tata, min_n = 0, prep = NULL),
  actividad_alta = list(sub = function(d) filter(d, act_alta), fo = fo_aj, min_n = min_estrato, prep = NULL),
  actividad_baja = list(sub = function(d) filter(d, !act_alta), fo = fo_aj, min_n = min_estrato, prep = NULL),
  actividad_alta_mas_tata = list(sub = function(d) filter(d, act_alta), fo = fo_tata, min_n = min_estrato, prep = NULL),
  actividad_baja_mas_tata = list(sub = function(d) filter(d, !act_alta), fo = fo_tata, min_n = min_estrato, prep = NULL)
)
coef_var <- map(repname, function(r) {
  imap(variantes, function(v, m) {
    message("Ajustando ", m, " ", r)
    fit_coefs(v$sub(filter(data, rep == r)), v$fo, v$min_n, v$prep) %>% mutate(rep = r, modelo = m)
  }) %>% list_rbind()
}) %>%
  list_rbind() %>%
  group_by(modelo, TF) %>%
  mutate(estimate = if (any(is.na(estimate))) NA_real_ else estimate) %>%
  ungroup()

cob_var <- bind_rows(
  cob %>% select(TF, rep, n_unidos, estimate) %>% mutate(modelo = "ajustado"),
  cob_cgi %>% filter(modelo %in% c("estrato_con_isla", "estrato_sin_isla")) %>% select(TF, rep, n_unidos, estimate, modelo),
  coef_var %>% select(TF, rep, n_unidos, estimate, modelo)
) %>%
  group_split(modelo, rep) %>%
  map(~ cobertura_stats(.x) %>% mutate(modelo = .x$modelo[1])) %>%
  list_rbind() %>%
  mutate(grupo = grupo_de(TF))
modelos_lab <- c(
  ajustado = "Ajustado", ajustado_mas_ocupacion = "Ajustado + ocupación por otros factores", ajustado_mas_tata = "Ajustado + TATA-box",
  estrato_con_isla = "Solo promotores con isla CpG", estrato_sin_isla = "Solo promotores sin isla CpG",
  actividad_alta = "Actividad alta (cuartil superior)", actividad_baja = "Actividad baja (resto)",
  actividad_alta_mas_tata = "Actividad alta + TATA-box", actividad_baja_mas_tata = "Actividad baja + TATA-box"
)
cob_var <- cob_var %>% mutate(modelo = factor(modelo, levels = names(modelos_lab))) %>% arrange(modelo, rep, TF)
write_tsv(cob_var, file.path(out_dir, "p7_p8_residuo_ventana_por_factor_variantes.tsv"))

n_var <- coef_var %>%
  group_by(modelo, rep) %>%
  summarise(n_factores_evaluables = sum(!is.na(estimate)), n_promotores = first(n_unidos + n_no_unidos), .groups = "drop")
print(n_var)

# --- p7: pausa y elongacion ------------------------------------------------

miembros_p7 <- cob_var %>%
  filter(grupo %in% c("pausa", "preiniciacion", "mediador"), modelo %in% c("ajustado", "ajustado_mas_ocupacion", "estrato_con_isla", "estrato_sin_isla")) %>%
  select(modelo, grupo, TF, rep, n_unidos, estimate, residuo, n_ventana, percentil_ventana, p_perm, p_perm_bh) %>%
  arrange(grupo, TF, modelo, rep)
write_tsv(miembros_p7, file.path(out_dir, "p7_miembros_grupos.tsv"))
print(as.data.frame(miembros_p7 %>% filter(modelo == "ajustado")), digits = 3)

tests_p7 <- run_tests(cob_var %>% filter(modelo %in% c("ajustado", "ajustado_mas_ocupacion", "estrato_con_isla", "estrato_sin_isla")), by = c("modelo", "rep"))
write_tsv(tests_p7, file.path(out_dir, "p7_test_grupos.tsv"))
print(as.data.frame(tests_p7 %>% select(modelo, rep, comparacion, variable, n_a, n_b, mediana_a, mediana_b, HL, pval)), digits = 3)

# modelo conjunto a nivel de promotor: indicadores de union de pausa y
# preiniciacion en el mismo lm (efectos parciales)
conj_tfs <- intersect(c("CDK9", "CCNT2", "NELFA", "NELFE", "SUPT5H", "AFF4", "BRD4", "TAF1", "TBP", "GTF2B", "GTF2F1"), tfs)
pic_nucleo <- c("TAF1", "TBP", "GTF2B")
conjunto <- map(repname, function(r) {
  d <- data %>% filter(rep == r)
  for (tf in conj_tfs) d[[tf]] <- as.numeric(d$seq_id %in% sets[[tf]])
  d$pic_nucleo <- as.numeric(rowSums(d[pic_nucleo]) == length(pic_nucleo))
  fo <- reformulate(c(conj_tfs, "log_hek", "mean_dnase"), "var_rank_sw")
  fit <- lm(fo, data = d)
  ci <- confint(fit)
  out <- tibble(analisis = "modelo_conjunto", term = conj_tfs, estimate = coef(fit)[conj_tfs], ci_low = ci[conj_tfs, 1], ci_high = ci[conj_tfs, 2],
    pval = summary(fit)$coefficients[conj_tfs, 4], n_unidos = map_dbl(conj_tfs, ~ sum(d[[.x]])))
  # cada factor por separado, en el subconjunto con y sin el nucleo de
  # preiniciacion (TAF1 + TBP + GTF2B) unido
  dentro <- map(c("CDK9", "CCNT2", "NELFA", "TAF1"), function(tf) {
    map(c(1, 0), function(k) {
      dd <- d %>% filter(pic_nucleo == k) %>% mutate(unido = .data[[tf]])
      if (tf %in% pic_nucleo && k == 1) return(NULL)
      fit <- lm(fo_aj, data = dd)
      tibble(analisis = if (k == 1) "dentro_de_nucleo_preiniciacion" else "fuera_de_nucleo_preiniciacion", term = tf,
        estimate = coef(fit)["unido"], ci_low = confint(fit)["unido", 1], ci_high = confint(fit)["unido", 2],
        pval = summary(fit)$coefficients["unido", 4], n_unidos = sum(dd$unido), n_subconjunto = nrow(dd))
    }) %>% list_rbind()
  }) %>% list_rbind()
  # solapamiento de CDK9 con el resto del modelo conjunto
  a <- d$CDK9 == 1
  solap <- map(setdiff(conj_tfs, "CDK9"), ~ tibble(analisis = "solapamiento_con_CDK9", term = .x, jaccard = sum(a & d[[.x]] == 1) / sum(a | d[[.x]] == 1),
    frac_de_CDK9 = sum(a & d[[.x]] == 1) / sum(a))) %>% list_rbind()
  bind_rows(out, dentro, solap) %>% mutate(rep = r, n_promotores = nrow(d), n_nucleo_preiniciacion = sum(d$pic_nucleo))
}) %>% list_rbind() %>% relocate(rep, analisis)
write_tsv(conjunto, file.path(out_dir, "p7_modelo_conjunto_maquinaria.tsv"))
print(as.data.frame(conjunto), digits = 3)

# --- Procedencia (p7: CDK9, CCNT2; p8: GTF2A2, SNAPC1) ------------------------

remap_lineas <- remap_raw %>%
  filter(TF %in% procedencia_tfs) %>%
  separate_rows(sample, sep = ",") %>%
  distinct(seq_id, TF, sample)
proc_lineas <- remap_lineas %>%
  group_by(TF) %>%
  mutate(n_lineas_tf = n_distinct(sample)) %>%
  group_by(TF, seq_id) %>%
  mutate(n_lineas_promotor = n_distinct(sample)) %>%
  group_by(TF, sample, n_lineas_tf) %>%
  summarise(promotores_library = n_distinct(seq_id), exclusivos_de_la_linea = sum(n_lineas_promotor == 1), .groups = "drop") %>%
  group_by(TF) %>%
  mutate(total_promotores_tf = n_distinct(remap_hits$seq_id[remap_hits$TF == TF[1]]), frac_del_total = promotores_library / total_promotores_tf) %>%
  ungroup() %>%
  arrange(TF, desc(promotores_library))
chipatlas <- read_tsv(path_chipatlas_explist, col_names = FALSE, col_select = 1:6, show_col_types = FALSE) %>%
  filter(X2 == "hg38", X4 %in% procedencia_tfs) %>%
  count(TF = X4, cell_type_chipatlas = X6, name = "experimentos_chipatlas")
efecto_linea <- remap_lineas %>%
  distinct(TF, sample) %>%
  pmap(function(TF, sample) {
    s <- remap_lineas$seq_id[remap_lineas$TF == TF & remap_lineas$sample == sample]
    map(repname, function(r) {
      d <- data %>% filter(rep == r)
      d$unido <- as.numeric(d$seq_id %in% s)
      if (sum(d$unido) < min_estrato) return(NULL)
      fit <- summary(lm(fo_aj, data = d))$coefficients["unido", ]
      ref <- cob %>% filter(rep == r, TF != !!TF, abs(n_unidos - sum(d$unido)) <= cov_window * sum(d$unido))
      tibble(TF, sample, rep = r, n_unidos = sum(d$unido), estimate = fit[1], pval = fit[4], n_ventana = nrow(ref), percentil_ventana = mean(ref$estimate <= fit[1]))
    }) %>% list_rbind()
  }) %>%
  list_rbind()
write_tsv(proc_lineas, file.path(out_dir, "p7_p8_procedencia_lineas_remap.tsv"))
write_tsv(chipatlas, file.path(out_dir, "p7_p8_experimentos_chipatlas.tsv"))
write_tsv(efecto_linea, file.path(out_dir, "p7_p8_efecto_por_linea_celular.tsv"))
print(proc_lineas, n = 80)
print(chipatlas %>% group_by(TF) %>% summarise(experimentos = sum(experimentos_chipatlas), tipos_celulares = n(), detalle = paste(cell_type_chipatlas, experimentos_chipatlas, collapse = "; ")))
print(as.data.frame(efecto_linea), digits = 3)

# --- p8: GTF2A2 y TATA-box ---------------------------------------------------

# a) TATA en los unidos, para los 982 (ubica a GTF2A2 en la distribucion)
tata_or <- map(repname, function(r) {
  d <- data %>% filter(rep == r)
  map(tfs, function(tf) {
    u <- d$seq_id %in% sets[[tf]]
    ft <- fisher.test(table(factor(u, c(FALSE, TRUE)), factor(d$TATA_EPD, c(FALSE, TRUE))))
    tibble(TF = tf, rep = r, n_unidos = sum(u), n_tata_unidos = sum(d$TATA_EPD[u]), frac_tata_unidos = mean(d$TATA_EPD[u]),
      frac_tata_no_unidos = mean(d$TATA_EPD[!u]), odds_ratio = unname(ft$estimate), or_ci_low = ft$conf.int[1], or_ci_high = ft$conf.int[2], pval_fisher = ft$p.value)
  }) %>% list_rbind()
}) %>%
  list_rbind() %>%
  group_by(rep) %>%
  mutate(rango_or_entre_982 = rank(-odds_ratio), pval_fisher_bh = p.adjust(pval_fisher, method = "BH")) %>%
  ungroup()
write_tsv(tata_or, file.path(out_dir, "p8_tata_por_factor.tsv"))
print(as.data.frame(tata_or %>% filter(TF %in% foco_p8)), digits = 3)
print(as.data.frame(tata_or %>% group_by(rep) %>% slice_min(rango_or_entre_982, n = 10)), digits = 3)

# b, c) efecto de los factores foco en cada modelo / subconjunto, con IC
focal_modelos <- list(
  ajustado = list(sub = identity, fo = fo_aj, term = "unido"),
  ajustado_mas_tata = list(sub = identity, fo = fo_tata, term = "unido"),
  ajustado_mas_ocupacion = list(sub = identity, fo = fo_ocup, term = "unido"),
  sin_tata = list(sub = function(d) filter(d, !TATA_EPD), fo = fo_aj, term = "unido"),
  con_tata = list(sub = function(d) filter(d, TATA_EPD), fo = fo_aj, term = "unido"),
  actividad_alta = list(sub = function(d) filter(d, act_alta), fo = fo_aj, term = "unido"),
  actividad_baja = list(sub = function(d) filter(d, !act_alta), fo = fo_aj, term = "unido"),
  actividad_alta_mas_tata = list(sub = function(d) filter(d, act_alta), fo = fo_tata, term = "unido"),
  interaccion_actividad = list(sub = identity, fo = var_rank_sw ~ unido * act_alta + log_hek + mean_dnase, term = "unido:act_altaTRUE")
)
focal <- map(repname, function(r) {
  map(union(foco_p8, extremos_pausa), function(tf) {
    imap(focal_modelos, function(m, nm) {
      d <- m$sub(filter(data, rep == r)) %>% mutate(unido = as.numeric(seq_id %in% sets[[tf]])) %>% prep_ocup()
      fit <- lm(m$fo, data = d)
      ci <- confint(fit)
      tibble(rep = r, TF = tf, modelo = nm, term = m$term, n_unidos = sum(d$unido), n_subconjunto = nrow(d),
        estimate = coef(fit)[m$term], ci_low = ci[m$term, 1], ci_high = ci[m$term, 2], pval = summary(fit)$coefficients[m$term, 4])
    }) %>% list_rbind()
  }) %>% list_rbind()
}) %>% list_rbind()
focal <- focal %>%
  left_join(cob_var %>% select(TF, rep, modelo, percentil_ventana, n_ventana, p_perm) %>% mutate(modelo = as.character(modelo)), by = c("TF", "rep", "modelo"))
write_tsv(focal, file.path(out_dir, "p8_efecto_focal_por_modelo.tsv"))
print(as.data.frame(focal), digits = 3)

# SNAPC1: genes de snRNA en la library y cerca de los promotores unidos
gt <- suppressMessages(AnnotationDbi::select(hs, keys = keys(hs, "ENTREZID"), columns = c("SYMBOL", "GENETYPE"), keytype = "ENTREZID"))
snrna <- gt %>% filter(GENETYPE == "snRNA")
loc <- suppressMessages(AnnotationDbi::select(hs, keys = snrna$ENTREZID, columns = "CHRLOC", keytype = "ENTREZID")) %>%
  filter(!is.na(CHRLOC), !str_detect(CHRLOCCHR, "_")) %>%
  transmute(SYMBOL = snrna$SYMBOL[match(ENTREZID, snrna$ENTREZID)], chr = paste0("chr", CHRLOCCHR), pos = abs(CHRLOC)) %>%
  distinct()
bed <- read_tsv("data/raw/library.bed", col_names = c("chr", "start", "end", "seq_id", "score", "strand"), show_col_types = FALSE)
lib_prom <- bed %>% filter(seq_id %in% data$seq_id)
dist_snrna <- lib_prom %>%
  inner_join(loc, by = "chr", relationship = "many-to-many") %>%
  mutate(dist = pmax(0, pmax(start - pos, pos - end))) %>%
  group_by(seq_id) %>%
  slice_min(dist, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(seq_id, snrna_cercano = SYMBOL, dist_snrna = dist)
snapc <- map(repname, function(r) {
  d <- data %>% filter(rep == r) %>% left_join(dist_snrna, by = "seq_id") %>% mutate(unido = seq_id %in% sets[["SNAPC1"]], cerca = replace_na(dist_snrna <= 2000, FALSE))
  excl <- d %>% filter(!cerca) %>% mutate(unido = as.numeric(unido))
  fit <- lm(fo_aj, data = excl)
  tibble(rep = r, n_promotores = nrow(d), n_genetype_snrna_library = sum(d$gene_sym %in% snrna$SYMBOL), n_snrna_hg38_org = n_distinct(loc$SYMBOL),
    n_promotores_a_2kb_de_snrna = sum(d$cerca), n_unidos_SNAPC1 = sum(d$unido), n_unidos_a_2kb_de_snrna = sum(d$unido & d$cerca),
    snrna_cercanos_unidos = paste(sort(unique(d$snrna_cercano[d$unido & d$cerca])), collapse = ", "),
    frac_tata_unidos = mean(d$TATA_EPD[d$unido]), frac_cgi_unidos = mean(d$CGI_frag[d$unido]), frac_cgi_library = mean(d$CGI_frag),
    estimate_sin_cercanos = coef(fit)["unido"], pval_sin_cercanos = summary(fit)$coefficients["unido", 4], n_unidos_sin_cercanos = sum(excl$unido))
}) %>% list_rbind()
write_tsv(snapc, file.path(out_dir, "p8_snapc1_snrna.tsv"))
print(as.data.frame(snapc), digits = 3)

# --- Figuras ----------------------------------------------------------------

tema <- ggpubr::theme_pubr(base_size = 14) + theme(strip.background = element_blank(), legend.position = "top")
etiquetas <- c(unlist(grupos[c("pausa", "preiniciacion")]), "SNAPC1")

plot_cob <- function(df, y, ylab, linea = TRUE) {
  dg <- df %>% filter(grupo != "resto")
  ggplot(df, aes(n_unidos, .data[[y]])) +
    geom_hline(yintercept = 0, linetype = "dotted", col = "grey40") +
    geom_point(data = ~ filter(.x, grupo == "resto"), col = "grey70", size = 0.9, alpha = 0.6) +
    { if (linea) geom_line(data = ~ arrange(.x, n_unidos), aes(y = ajuste_cobertura), col = "grey20", linewidth = 0.7) } +
    geom_point(data = dg, aes(col = grupo), size = 2.4) +
    geom_text_repel(data = ~ filter(.x, TF %in% etiquetas), aes(label = TF, col = grupo), size = 3.2, min.segment.length = 0, max.overlaps = 50, show.legend = FALSE) +
    scale_x_log10() +
    scale_color_manual(values = grupo_clr, labels = grupo_lab) +
    facet_wrap(~rep) +
    tema +
    labs(x = "Promotores unidos (ReMap, escala log10)", y = ylab, col = NULL)
}
p_cob <- plot_cob(cob_var %>% filter(modelo == "ajustado"), "estimate", "Efecto ajustado de la unión sobre\nel rango de varianza (1-100)") +
  labs(caption = "Línea: loess robusto sobre log10(n unidos). Mediador sin etiqueta.")
save_jpg_pdf(p_cob, "p7_coef_vs_cobertura", 13, 6.5)

p_grp <- cob_var %>%
  filter(modelo %in% c("ajustado", "ajustado_mas_ocupacion")) %>%
  mutate(modelo = factor(modelos_lab[as.character(modelo)], levels = modelos_lab), grupo = factor(grupo, levels = names(grupo_lab))) %>%
  ggplot(aes(grupo, percentil_ventana, col = grupo)) +
  geom_boxplot(data = ~ filter(.x, grupo == "resto"), outlier.shape = NA, col = "grey50") +
  geom_jitter(data = ~ filter(.x, grupo == "resto"), width = 0.2, size = 0.5, alpha = 0.3) +
  geom_point(data = ~ filter(.x, grupo != "resto"), size = 2.2) +
  geom_text_repel(data = ~ filter(.x, TF %in% etiquetas), aes(label = TF), size = 2.8, direction = "y", nudge_x = 0.35, min.segment.length = 0, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = grupo_clr, labels = grupo_lab) +
  scale_x_discrete(labels = c(pausa = "Pausa /\nelongación", preiniciacion = "Preiniciación", mediador = "Mediador", resto = "Resto")) +
  facet_grid(modelo ~ rep) +
  tema +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Percentil en la ventana de cobertura ±25%\n(0 = mayor reducción del ruido)")
save_jpg_pdf(p_grp, "p7_percentil_por_grupo", 12, 10)

orden_tf <- cob_var %>% filter(modelo == "ajustado", grupo %in% c("pausa", "preiniciacion")) %>% group_by(TF, grupo) %>% summarise(m = mean(percentil_ventana), .groups = "drop") %>% arrange(grupo, m)
mods_pct <- c("ajustado", "ajustado_mas_ocupacion", "estrato_con_isla", "estrato_sin_isla")
p_pct <- cob_var %>%
  filter(grupo %in% c("pausa", "preiniciacion"), modelo %in% mods_pct) %>%
  mutate(modelo = factor(modelos_lab[as.character(modelo)], levels = modelos_lab[mods_pct]), TF = factor(TF, levels = orden_tf$TF)) %>%
  ggplot(aes(percentil_ventana, TF, col = modelo)) +
  geom_vline(xintercept = 0.5, col = "grey60") +
  geom_point(size = 2.4, position = position_dodge(0.6)) +
  scale_color_manual(values = c("black", "#7A5195", "#1B8C8E", "#AD343E")) +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(~rep) +
  tema +
  theme(legend.direction = "vertical", axis.text.y = element_text(color = ifelse(orden_tf$grupo == "pausa", grupo_clr["pausa"], grupo_clr["preiniciacion"]))) +
  labs(x = "Percentil en la ventana de cobertura ±25% (0 = mayor reducción del ruido)", y = NULL, col = NULL,
    caption = "Rojo: pausa / elongación. Verde azulado: preiniciación. Sin punto = factor no evaluable en el estrato.")
save_jpg_pdf(p_pct, "p7_percentil_ventana_miembros", 12, 10)

mods_foc <- c(ajustado = "Ajustado", ajustado_mas_tata = "+ TATA-box", sin_tata = "Solo sin TATA", con_tata = "Solo con TATA",
  actividad_alta = "Actividad alta", actividad_baja = "Actividad baja", actividad_alta_mas_tata = "Actividad alta + TATA-box")
p_foc <- focal %>%
  filter(TF %in% foco_p8, modelo %in% names(mods_foc)) %>%
  mutate(modelo = factor(mods_foc[modelo], levels = rev(mods_foc)), TF = factor(TF, levels = foco_p8)) %>%
  ggplot(aes(estimate, modelo)) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0) +
  geom_point(aes(col = TF == "GTF2A2"), size = 2.4) +
  scale_color_manual(values = c(`TRUE` = "#AD343E", `FALSE` = "grey20"), guide = "none") +
  facet_grid(rep ~ TF) +
  tema +
  labs(x = "Efecto de la unión sobre el rango de varianza (IC 95%)", y = NULL)
save_jpg_pdf(p_foc, "p8_efecto_focal_por_modelo", 15, 7.5)

p_or <- tata_or %>%
  mutate(grupo = grupo_de(TF)) %>%
  ggplot(aes(n_unidos, odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dotted", col = "grey40") +
  geom_point(col = "grey70", size = 0.9, alpha = 0.6) +
  geom_point(data = ~ filter(.x, TF %in% foco_p8), col = "#AD343E", size = 2.4) +
  geom_text_repel(data = ~ filter(.x, TF %in% foco_p8 | rango_or_entre_982 <= 5), aes(label = TF), size = 3.2, min.segment.length = 0) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~rep) +
  tema +
  labs(x = "Promotores unidos (ReMap, escala log10)", y = "Odds ratio de TATA-box canónica\nen los unidos (escala log10)")
save_jpg_pdf(p_or, "p8_tata_odds_ratio_vs_cobertura", 13, 6.5)

message("Exploratorio p7/p8 guardado en ", out_dir)
