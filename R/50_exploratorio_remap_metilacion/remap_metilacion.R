# Exploratorio sobre R/48_exploratorio_remap_ruido_ajustado: ¿hay un eje
# ligado a la metilacion del ADN en el efecto de la union de cada factor
# (ReMap) sobre el ruido, una vez controlada la cobertura (n de promotores
# unidos)?
#
# Parte del coeficiente de `unido` del modelo ajustado de p1
# (var_rank_sw ~ unido + log10(hek_tpm + 1) + mean_dnase) para los 982
# factores de remap_tf_noise_auc.tsv. Puntos del pedido:
#
# p1 Estadistico controlado por cobertura, por replica:
#   - residuo de un loess robusto de coef sobre log10(n_unidos) (span 0.5,
#     family = "symmetric"). La relacion es plana hasta ~3000 promotores y
#     cae cada vez mas rapido por encima de ~4000 (el grupo no unido se
#     reduce a los promotores mas atipicos): una recta en log10(n) deja un
#     sesgo de ~+1 en 1000-5000 y de -1 a -2 en 6000-9000, justo donde estan
#     los dos grupos. El residuo de MASS::rlm lineal se guarda solo como
#     referencia (residuo_rlm);
#   - posicion en la ventana de cobertura +-25% (mismo procedimiento que
#     remap_ruido_cobertura.R / nfy_especificidad). Percentil = fraccion de
#     los factores de la ventana con coeficiente <= al del focal (0 = mayor
#     reduccion del ruido). p por permutacion exacto: bajo la hipotesis de
#     que el focal es intercambiable con los factores de su ventana, su rango
#     es uniforme; p bilateral = 2 * min(rango inferior, rango superior) /
#     (n_ventana + 1). BH sobre los 982 factores en cada replica.
# p2 Test de grupo (Wilcoxon de rangos, grupo vs. el resto) sobre el residuo
#   y el percentil de ventana, en cada replica. Grupos definidos a priori.
# p3 GSEA (mismos parametros que remap_tf_noise_gsea) sobre el ranking por
#   residuo, comparado con los rankings AUC original y ajustado de p3.
# p4 Control por isla CpG (CGI_frag, la definicion por composicion del
#   fragmento que usa todo R3):
#   a) modelo ajustado + CGI_frag;
#   b) modelo ajustado dentro de los promotores con y sin isla (var_rank_sw
#      global, restringido al estrato; se exige >= 50 unidos y >= 50 no
#      unidos dentro del estrato en ambas replicas).
#   En cada version se recalculan residuo, ventana y test de grupo.
# p5 Procedencia: lineas celulares de los picos en ReMap (nombre TF:biotipo
#   del bed no redundante) y efecto por linea celular; experimentos de
#   ChIP-Atlas como referencia cruzada (el bed no redundante de ReMap no
#   guarda el numero de experimentos).
# p6 Solapamiento DNMT3B / ZBTB33: conjuntos unidos y no unidos, y efecto
#   por categoria (solo uno, ambos, ninguno).
#
# Salida en figures/exploratorio_remap_metilacion/ (exploratorio, jpg y pdf),
# con prefijo por punto. Requiere data/processed/{remap_tf_hits,
# activity_stats_highconf,prom_df}.tsv, remap_tf_noise_auc.tsv y, de
# figures/exploratorio_remap_ruido_ajustado/, p1_coeficientes_por_factor.tsv
# y p3_gsea_todos_los_terminos.tsv. Tarda unos minutos. Run from the
# TesisDoc repo root.

library(tidyverse)
library(clusterProfiler)
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

BiocParallel::register(BiocParallel::SerialParam())

in_dir <- "figures/exploratorio_remap_ruido_ajustado"
out_dir <- "figures/exploratorio_remap_metilacion"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hs <- org.Hs.eg.db

cov_window <- 0.25
min_estrato <- 50

grupos_def <- list(
  metilacion = c("DNMT1", "DNMT3A", "DNMT3B", "ZBTB33", "MBD1", "MBD2", "MBD3", "MBD4", "MECP2", "UHRF1", "TET1", "TET2", "TET3"),
  cpg_no_metilado = c("MAZ", "SP1", "SP2", "SP3", "ZBTB14", "E2F1", "E2F4", "NRF1", "GABPA", "KLF16", "CXXC1", "KDM2B")
)
grupo_lab <- c(metilacion = "Asociados a metilación", cpg_no_metilado = "Asociados a CpG no metilado", resto = "Resto")
grupo_clr <- c(metilacion = "#AD343E", cpg_no_metilado = "#1B8C8E", resto = "grey70")
procedencia_tfs <- c("DNMT3B", "ZBTB33", "MAZ", "SP2", "ZBTB14")

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = "")

# --- Datos ---------------------------------------------------------------

auc_orig <- read_tsv(file.path(fig_dir("remap_tf_noise_gsea"), "remap_tf_noise_auc.tsv"), show_col_types = FALSE)
tfs <- unique(auc_orig$feature)
coef_p1 <- read_tsv(file.path(in_dir, "p1_coeficientes_por_factor.tsv"), show_col_types = FALSE)

remap_raw <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
remap_hits <- remap_raw %>% distinct(seq_id, TF)
sets <- remap_hits %>% filter(TF %in% tfs) %>% with(split(seq_id, TF))

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  mutate(log_hek = log10(hek_tpm + 1)) %>%
  select(seq_id, rep, var_rank_sw, log_hek, mean_dnase, CGI_frag)
stopifnot(!anyNA(data))
repname <- unique(data$rep)

# grupos: solo factores entre los 982 (>100 promotores unidos en ambas
# replicas, criterio de remap_tf_noise_gsea); los demas se reportan con motivo
n_lib <- remap_hits %>%
  filter(TF %in% unlist(grupos_def)) %>%
  inner_join(distinct(data, seq_id, rep), by = "seq_id", relationship = "many-to-many") %>%
  count(TF, rep) %>%
  pivot_wider(names_from = rep, values_from = n, names_prefix = "n_")
grupos_tab <- enframe(grupos_def, "grupo", "TF") %>%
  unnest(TF) %>%
  left_join(n_lib, by = "TF") %>%
  mutate(
    en_remap_library = TF %in% remap_hits$TF,
    incluido = TF %in% tfs,
    motivo = case_when(
      incluido ~ "incluido",
      !en_remap_library ~ "sin picos de ReMap sobre promotores de la library",
      TRUE ~ "100 promotores unidos o menos en alguna réplica (umbral de remap_tf_noise_gsea)"
    )
  )
write_tsv(grupos_tab, file.path(out_dir, "p2_grupos_definicion.tsv"))
print(grupos_tab, n = 30)
grupos <- map(grupos_def, ~ intersect(.x, tfs))
grupo_de <- function(tf) case_when(tf %in% grupos$metilacion ~ "metilacion", tf %in% grupos$cpg_no_metilado ~ "cpg_no_metilado", TRUE ~ "resto")

# --- Funciones -------------------------------------------------------------

# coeficiente de `unido` para todos los factores en un subconjunto de
# promotores; NA si hay menos de `min_n` unidos o no unidos
fit_coefs <- function(d, fo, min_n = 0) {
  map(tfs, function(tf) {
    d$unido <- as.numeric(d$seq_id %in% sets[[tf]])
    n1 <- sum(d$unido)
    n0 <- nrow(d) - n1
    if (n1 < min_n || n0 < min_n) {
      return(tibble(TF = tf, n_unidos = n1, n_no_unidos = n0, estimate = NA_real_, se = NA_real_, pval = NA_real_))
    }
    co <- summary(lm(fo, data = d))$coefficients["unido", ]
    tibble(TF = tf, n_unidos = n1, n_no_unidos = n0, estimate = co[1], se = co[2], pval = co[4])
  }) %>% list_rbind()
}

# residuo por cobertura (loess robusto sobre log10(n); rlm lineal como
# referencia) y posicion en la ventana +-25%, dentro de una replica
cobertura_stats <- function(df) {
  df <- df %>% filter(!is.na(estimate))
  df$log_n <- log10(df$n_unidos)
  fit_rlm <- MASS::rlm(estimate ~ log_n, data = df, maxit = 100)
  fit_lo <- loess(estimate ~ log_n, data = df, family = "symmetric", span = 0.5)
  df$ajuste_cobertura <- predict(fit_lo, df)
  df$residuo <- df$estimate - df$ajuste_cobertura
  df$residuo_z <- df$residuo / mad(df$residuo)
  df$residuo_rlm <- df$estimate - predict(fit_rlm, df)
  win <- map(seq_len(nrow(df)), function(i) {
    ref <- df$estimate[-i][abs(df$n_unidos[-i] - df$n_unidos[i]) <= cov_window * df$n_unidos[i]]
    x <- df$estimate[i]
    n_le <- sum(ref <= x)
    n_ge <- sum(ref >= x)
    tibble(
      n_ventana = length(ref), percentil_ventana = n_le / length(ref),
      mediana_ventana = median(ref),
      p_perm = min(1, 2 * (min(n_le, n_ge) + 1) / (length(ref) + 1))
    )
  }) %>% list_rbind()
  bind_cols(df, win) %>%
    mutate(
      p_perm_bh = p.adjust(p_perm, method = "BH"),
      rlm_intercepto = coef(fit_rlm)[1], rlm_pendiente = coef(fit_rlm)[2]
    )
}

# Wilcoxon de cada grupo contra el resto (y metilacion vs. CpG no metilado)
test_grupos <- function(df, var) {
  comps <- list(
    metilacion_vs_resto = list(a = grupos$metilacion, b = setdiff(df$TF, grupos$metilacion)),
    metilacion_sin_DNMT3B_ZBTB33_vs_resto = list(a = setdiff(grupos$metilacion, c("DNMT3B", "ZBTB33")), b = setdiff(df$TF, grupos$metilacion)),
    cpg_no_metilado_vs_resto = list(a = grupos$cpg_no_metilado, b = setdiff(df$TF, grupos$cpg_no_metilado)),
    metilacion_vs_cpg_no_metilado = list(a = grupos$metilacion, b = grupos$cpg_no_metilado)
  )
  imap(comps, function(cmp, nm) {
    xa <- df[[var]][df$TF %in% cmp$a & !is.na(df[[var]])]
    xb <- df[[var]][df$TF %in% cmp$b & !is.na(df[[var]])]
    w <- wilcox.test(xa, xb, conf.int = TRUE, exact = FALSE)
    tibble(
      comparacion = nm, variable = var, n_a = length(xa), n_b = length(xb),
      mediana_a = median(xa), mediana_b = median(xb),
      # rank-biserial: > 0 = el primer grupo tiene valores mas altos (mas ruido)
      rank_biserial = 2 * w$statistic / (length(xa) * length(xb)) - 1,
      HL = unname(w$estimate), pval = w$p.value
    )
  }) %>% list_rbind()
}

# --- p1: estadistico controlado por cobertura -------------------------------

coef_aj <- coef_p1 %>% transmute(TF, rep, n_unidos, estimate = estimate_ajustado, pval_corr_ajustado)
cob <- coef_aj %>% group_split(rep) %>% map(cobertura_stats) %>% list_rbind() %>% mutate(grupo = grupo_de(TF))
write_tsv(cob, file.path(out_dir, "p1_residuo_ventana_por_factor.tsv"))

rho_cob <- cob %>%
  group_by(rep) %>%
  summarise(
    rho_coef_n = cor(estimate, n_unidos, method = "spearman"),
    rho_residuo_n = cor(residuo, n_unidos, method = "spearman"),
    rho_residuo_rlm_n = cor(residuo_rlm, n_unidos, method = "spearman"),
    rho_percentil_n = cor(percentil_ventana, n_unidos, method = "spearman"),
    rho_residuo_percentil = cor(residuo, percentil_ventana, method = "spearman"),
    rlm_pendiente = first(rlm_pendiente), mad_residuo = mad(residuo),
    n_bh_005 = sum(p_perm_bh < 0.05), n_bh_01 = sum(p_perm_bh < 0.1), p_perm_min = min(p_perm)
  )
write_tsv(rho_cob, file.path(out_dir, "p1_resumen_cobertura.tsv"))
print(as.data.frame(rho_cob), digits = 3)

miembros_p1 <- cob %>%
  filter(grupo != "resto") %>%
  select(grupo, TF, rep, n_unidos, estimate, residuo, residuo_z, residuo_rlm, n_ventana, percentil_ventana, mediana_ventana, p_perm, p_perm_bh) %>%
  arrange(grupo, TF, rep)
write_tsv(miembros_p1, file.path(out_dir, "p1_miembros_grupos.tsv"))
print(as.data.frame(miembros_p1), digits = 3)

# extremos de cada replica por percentil de ventana (para ver quienes acompañan
# a DNMT3B/ZBTB33 y a MAZ/SP2)
extremos <- cob %>%
  select(TF, rep, n_unidos, estimate, residuo, percentil_ventana, p_perm) %>%
  pivot_wider(names_from = rep, values_from = c(n_unidos, estimate, residuo, percentil_ventana, p_perm), names_glue = "{.value}_{str_remove(rep, ' ')}") %>%
  mutate(percentil_medio = (percentil_ventana_Rep1 + percentil_ventana_Rep2) / 2, residuo_medio = (residuo_Rep1 + residuo_Rep2) / 2, grupo = grupo_de(TF)) %>%
  arrange(desc(residuo_medio))
write_tsv(extremos, file.path(out_dir, "p1_ranking_residuo_ambas_replicas.tsv"))

# --- p2: test de grupo -------------------------------------------------------

tests_p2 <- cob %>%
  group_split(rep) %>%
  map(function(df) map(c("residuo", "residuo_rlm", "percentil_ventana"), ~ test_grupos(df, .x)) %>% list_rbind() %>% mutate(rep = df$rep[1])) %>%
  list_rbind() %>%
  relocate(rep)
write_tsv(tests_p2, file.path(out_dir, "p2_test_grupos.tsv"))
print(as.data.frame(tests_p2), digits = 3)

# --- p3: GSEA sobre el ranking por residuo ----------------------------------

ids <- AnnotationDbi::select(hs, keys = tfs, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>% distinct()
run_gsea <- function(df) {
  pre <- df %>%
    arrange(desc(residuo)) %>%
    left_join(ids, by = c("TF" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  gseGO(setNames(pre$residuo, pre$ENTREZID), ont = "CC", OrgDb = "org.Hs.eg.db", pvalueCutoff = 1, seed = TRUE)
}
set.seed(2024)
gsea_res <- map(set_names(repname), ~ run_gsea(cob %>% filter(rep == .x)))
gsea_res_df <- imap(gsea_res, ~ as.data.frame(.x) %>% mutate(ranking = "residuo_cobertura", rep = .y)) %>%
  list_rbind() %>%
  as_tibble() %>%
  select(ranking, rep, ID, Description, setSize, NES, pvalue, p.adjust, core_enrichment)
gsea_p3 <- read_tsv(file.path(in_dir, "p3_gsea_todos_los_terminos.tsv"), show_col_types = FALSE) %>% filter(ranking %in% c("AUC_original", "ajustado"))
gsea_all <- bind_rows(gsea_p3, gsea_res_df)
write_tsv(gsea_res_df, file.path(out_dir, "p3_gsea_residuo_todos_los_terminos.tsv"))

patron_tema <- "methyl|chromatin|heterochromatin|histone|CpG|nucleosome|DNA repair|replication fork"
terminos_ref <- c("GO:0044665", "GO:0071339", "GO:0140513", "GO:1902494", "GO:0032991", "GO:0000785", "GO:1990234")
terminos <- union(
  terminos_ref,
  gsea_all %>% filter(p.adjust < 0.1 | (ranking == "residuo_cobertura" & str_detect(Description, regex(patron_tema, ignore_case = TRUE)) & pvalue < 0.05)) %>% pull(ID)
)
gsea_comp <- gsea_all %>%
  filter(ID %in% terminos) %>%
  mutate(rep = str_remove(rep, " ")) %>%
  select(ID, Description, ranking, rep, NES, p.adjust) %>%
  pivot_wider(id_cols = c(ID, Description), names_from = c(ranking, rep), values_from = c(NES, p.adjust), names_glue = "{.value}_{ranking}_{rep}") %>%
  mutate(tema_metilacion_cromatina = str_detect(Description, regex(patron_tema, ignore_case = TRUE))) %>%
  arrange(p.adjust_residuo_cobertura_Rep1)
col_orden <- c("ID", "Description", "tema_metilacion_cromatina", unlist(map(c("AUC_original", "ajustado", "residuo_cobertura"), function(rk) map(c("Rep1", "Rep2"), ~ paste0(c("NES_", "p.adjust_"), rk, "_", .x)))))
gsea_comp <- gsea_comp %>% select(any_of(col_orden))
write_tsv(gsea_comp, file.path(out_dir, "p3_gsea_comparacion.tsv"))

# chequeo: el NES de "chromatin" (GO:0000785) sobre el residuo, ¿lo arrastran
# los factores de baja cobertura (residuos mas dispersos)? Wilcoxon anotados
# vs. no anotados dentro de cada tercio de cobertura
chrom_tfs <- suppressMessages(AnnotationDbi::select(hs, keys = "GO:0000785", columns = "SYMBOL", keytype = "GOALL"))$SYMBOL %>% unique()
chrom_tercil <- cob %>%
  mutate(chromatin = TF %in% chrom_tfs, tercil_cobertura = ntile(n_unidos, 3), .by = rep) %>%
  group_by(rep, tercil_cobertura) %>%
  summarise(
    n_min = min(n_unidos), n_max = max(n_unidos), n_chromatin = sum(chromatin), n_otros = sum(!chromatin),
    mediana_residuo_chromatin = median(residuo[chromatin]), mediana_residuo_otros = median(residuo[!chromatin]),
    p_residuo = wilcox.test(residuo[chromatin], residuo[!chromatin])$p.value,
    mediana_percentil_chromatin = median(percentil_ventana[chromatin]), mediana_percentil_otros = median(percentil_ventana[!chromatin]),
    p_percentil = wilcox.test(percentil_ventana[chromatin], percentil_ventana[!chromatin])$p.value,
    .groups = "drop"
  )
write_tsv(chrom_tercil, file.path(out_dir, "p3_chromatin_por_tercil_cobertura.tsv"))
print(as.data.frame(chrom_tercil), digits = 3)
print(as.data.frame(gsea_comp %>% head(40)), digits = 3)

# --- p4: control por isla CpG ----------------------------------------------

fo_aj <- var_rank_sw ~ unido + log_hek + mean_dnase
fo_cgi <- var_rank_sw ~ unido + log_hek + mean_dnase + CGI_frag

frac_cgi <- map(repname, function(r) {
  d <- data %>% filter(rep == r)
  tibble(TF = tfs, rep = r, frac_cgi_unidos = map_dbl(tfs, ~ mean(d$CGI_frag[d$seq_id %in% sets[[.x]]])), frac_cgi_no_unidos = map_dbl(tfs, ~ mean(d$CGI_frag[!d$seq_id %in% sets[[.x]]])))
}) %>% list_rbind()
frac_cgi_lib <- data %>% group_by(rep) %>% summarise(frac_cgi_library = mean(CGI_frag), n = n())

variantes_cgi <- list(
  ajustado_mas_cgi = list(sub = function(d) d, fo = fo_cgi, min_n = 0),
  estrato_con_isla = list(sub = function(d) filter(d, CGI_frag), fo = fo_aj, min_n = min_estrato),
  estrato_sin_isla = list(sub = function(d) filter(d, !CGI_frag), fo = fo_aj, min_n = min_estrato)
)
coef_cgi <- map(repname, function(r) {
  imap(variantes_cgi, function(v, m) fit_coefs(v$sub(filter(data, rep == r)), v$fo, v$min_n) %>% mutate(rep = r, modelo = m)) %>% list_rbind()
}) %>% list_rbind()
# en los estratos, solo factores evaluables en ambas replicas
coef_cgi <- coef_cgi %>%
  group_by(modelo, TF) %>%
  mutate(estimate = if (any(is.na(estimate))) NA_real_ else estimate) %>%
  ungroup()

cob_cgi <- bind_rows(
  cob %>% select(TF, rep, n_unidos, estimate) %>% mutate(modelo = "ajustado"),
  coef_cgi %>% select(TF, rep, n_unidos, estimate, modelo)
) %>%
  group_split(modelo, rep) %>%
  map(~ cobertura_stats(.x) %>% mutate(modelo = .x$modelo[1])) %>%
  list_rbind() %>%
  mutate(grupo = grupo_de(TF)) %>%
  left_join(frac_cgi, by = c("TF", "rep"))
modelos_cgi <- c(ajustado = "Ajustado", ajustado_mas_cgi = "Ajustado + isla CpG", estrato_con_isla = "Solo promotores con isla CpG", estrato_sin_isla = "Solo promotores sin isla CpG")
cob_cgi <- cob_cgi %>% mutate(modelo = factor(modelo, levels = names(modelos_cgi))) %>% arrange(modelo, rep, TF)
write_tsv(cob_cgi, file.path(out_dir, "p4_residuo_ventana_por_factor_cgi.tsv"))

n_estrato <- coef_cgi %>%
  group_by(modelo, rep) %>%
  summarise(n_factores_evaluables = sum(!is.na(estimate)), n_promotores = first(n_unidos + n_no_unidos), .groups = "drop") %>%
  left_join(frac_cgi_lib, by = "rep")
write_tsv(n_estrato, file.path(out_dir, "p4_n_por_estrato.tsv"))
print(n_estrato)

miembros_p4 <- cob_cgi %>%
  filter(grupo != "resto") %>%
  select(modelo, grupo, TF, rep, n_unidos, frac_cgi_unidos, frac_cgi_no_unidos, estimate, residuo, percentil_ventana, n_ventana, p_perm) %>%
  arrange(grupo, TF, modelo, rep)
write_tsv(miembros_p4, file.path(out_dir, "p4_miembros_grupos_cgi.tsv"))
print(as.data.frame(miembros_p4 %>% filter(TF %in% c("DNMT3B", "ZBTB33", "MAZ", "SP2", "ZBTB14", "E2F1"))), digits = 3)

tests_p4 <- cob_cgi %>%
  group_split(modelo, rep) %>%
  map(function(df) map(c("residuo", "percentil_ventana"), ~ test_grupos(df, .x)) %>% list_rbind() %>% mutate(modelo = df$modelo[1], rep = df$rep[1])) %>%
  list_rbind() %>%
  relocate(modelo, rep)
write_tsv(tests_p4, file.path(out_dir, "p4_test_grupos_cgi.tsv"))
print(as.data.frame(tests_p4 %>% filter(variable == "residuo")), digits = 3)

# composicion CGI de los conjuntos unidos de los miembros, frente a la library
print(frac_cgi_lib)

# --- p5: procedencia ---------------------------------------------------------

# el nombre de cada pico del bed no redundante es TF:linea1,linea2,... (un
# pico fusionado puede venir de varias lineas)
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
write_tsv(proc_lineas, file.path(out_dir, "p5_procedencia_lineas_remap.tsv"))
write_tsv(chipatlas, file.path(out_dir, "p5_experimentos_chipatlas.tsv"))
print(proc_lineas, n = 60)
print(chipatlas %>% group_by(TF) %>% summarise(experimentos = sum(experimentos_chipatlas), tipos_celulares = n()))

# efecto por linea celular: union = pico en esa linea (>100 unidos en ambas
# replicas), percentil en la ventana de cobertura contra los 982 factores
efecto_linea <- remap_lineas %>%
  distinct(TF, sample) %>%
  pmap(function(TF, sample) {
    s <- remap_lineas$seq_id[remap_lineas$TF == TF & remap_lineas$sample == sample]
    map(repname, function(r) {
      d <- data %>% filter(rep == r)
      d$unido <- as.numeric(d$seq_id %in% s)
      if (sum(d$unido) <= 100) return(NULL)
      fit <- summary(lm(fo_aj, data = d))$coefficients["unido", ]
      ref <- cob %>% filter(rep == r, TF != !!TF, abs(n_unidos - sum(d$unido)) <= cov_window * sum(d$unido))
      tibble(TF, sample, rep = r, n_unidos = sum(d$unido), frac_cgi_unidos = mean(d$CGI_frag[d$unido == 1]), estimate = fit[1], pval = fit[4], n_ventana = nrow(ref), percentil_ventana = mean(ref$estimate <= fit[1]))
    }) %>% list_rbind()
  }) %>%
  list_rbind()
write_tsv(efecto_linea, file.path(out_dir, "p5_efecto_por_linea_celular.tsv"))
print(as.data.frame(efecto_linea), digits = 3)

# --- p6: solapamiento DNMT3B / ZBTB33 ----------------------------------------

par <- c("DNMT3B", "ZBTB33")
solap <- map(repname, function(r) {
  d <- data %>% filter(rep == r)
  a <- d$seq_id %in% sets[[par[1]]]
  b <- d$seq_id %in% sets[[par[2]]]
  ft <- fisher.test(table(a, b))
  tibble(
    rep = r, n_promotores = nrow(d), n_DNMT3B = sum(a), n_ZBTB33 = sum(b), n_ambos = sum(a & b), esperado_ambos = sum(a) * sum(b) / nrow(d),
    jaccard_unidos = sum(a & b) / sum(a | b), frac_DNMT3B_en_ZBTB33 = sum(a & b) / sum(a), frac_ZBTB33_en_DNMT3B = sum(a & b) / sum(b),
    n_no_unidos_DNMT3B = sum(!a), n_no_unidos_ZBTB33 = sum(!b), n_no_unidos_ambos = sum(!a & !b),
    jaccard_no_unidos = sum(!a & !b) / sum(!a | !b), odds_ratio = unname(ft$estimate), pval_fisher = ft$p.value
  )
}) %>% list_rbind()

# Jaccard de DNMT3B y de ZBTB33 con todos los otros factores (unidos y no
# unidos) para ubicar el solapamiento entre ellos
jacc_ref <- map(repname, function(r) {
  ids_r <- unique(data$seq_id[data$rep == r])
  map(par, function(f) {
    a <- ids_r %in% sets[[f]]
    map(setdiff(tfs, f), function(o) {
      b <- ids_r %in% sets[[o]]
      tibble(rep = r, focal = f, otro = o, jaccard_unidos = sum(a & b) / sum(a | b), jaccard_no_unidos = sum(!a & !b) / sum(!a | !b))
    }) %>% list_rbind()
  }) %>% list_rbind()
}) %>%
  list_rbind() %>%
  group_by(rep, focal) %>%
  mutate(rango_jaccard_unidos = rank(-jaccard_unidos), rango_jaccard_no_unidos = rank(-jaccard_no_unidos)) %>%
  ungroup()
write_tsv(jacc_ref %>% group_by(rep, focal) %>% slice_min(rango_jaccard_unidos, n = 15) %>% ungroup(), file.path(out_dir, "p6_jaccard_top15.tsv"))
solap <- solap %>%
  left_join(jacc_ref %>% filter(focal == "DNMT3B", otro == "ZBTB33") %>% select(rep, rango_ZBTB33_entre_981_para_DNMT3B = rango_jaccard_unidos), by = "rep") %>%
  left_join(jacc_ref %>% filter(focal == "ZBTB33", otro == "DNMT3B") %>% select(rep, rango_DNMT3B_entre_981_para_ZBTB33 = rango_jaccard_unidos), by = "rep")
write_tsv(solap, file.path(out_dir, "p6_solapamiento.tsv"))
print(as.data.frame(solap), digits = 3)

# efecto por categoria (referencia: ninguno de los dos) y efecto de cada uno
# sin los promotores unidos por el otro
cat_lab <- c(ninguno = "Ninguno", solo_DNMT3B = "Solo DNMT3B", solo_ZBTB33 = "Solo ZBTB33", ambos = "Ambos")
efecto_cat <- map(repname, function(r) {
  d <- data %>%
    filter(rep == r) %>%
    mutate(
      a = seq_id %in% sets[[par[1]]], b = seq_id %in% sets[[par[2]]],
      categoria = factor(case_when(a & b ~ "ambos", a ~ "solo_DNMT3B", b ~ "solo_ZBTB33", TRUE ~ "ninguno"), levels = names(cat_lab))
    )
  fit <- lm(var_rank_sw ~ categoria + log_hek + mean_dnase, data = d)
  ci <- confint(fit)
  cat_tab <- tibble(term = names(coef(fit)), estimate = coef(fit), ci_low = ci[, 1], ci_high = ci[, 2], pval = summary(fit)$coefficients[, 4]) %>%
    filter(str_detect(term, "^categoria")) %>%
    mutate(categoria = str_remove(term, "categoria"), analisis = "categoria_vs_ninguno") %>%
    left_join(d %>% count(categoria = as.character(categoria), name = "n") , by = "categoria") %>%
    left_join(d %>% group_by(categoria = as.character(categoria)) %>% summarise(frac_cgi = mean(CGI_frag), media_var_rank = mean(var_rank_sw)), by = "categoria")
  sin_otro <- map2(par, rev(par), function(f, o) {
    dd <- d %>% filter(!seq_id %in% sets[[o]]) %>% mutate(unido = as.numeric(seq_id %in% sets[[f]]))
    fit <- lm(fo_aj, data = dd)
    tibble(categoria = paste0(f, "_sin_promotores_de_", o), analisis = "sin_el_otro", estimate = coef(fit)["unido"], ci_low = confint(fit)["unido", 1], ci_high = confint(fit)["unido", 2], pval = summary(fit)$coefficients["unido", 4], n = sum(dd$unido))
  }) %>% list_rbind()
  bind_rows(cat_tab %>% select(-term), sin_otro) %>% mutate(rep = r)
}) %>% list_rbind() %>% relocate(rep, analisis, categoria)
write_tsv(efecto_cat, file.path(out_dir, "p6_efecto_por_categoria.tsv"))
print(as.data.frame(efecto_cat), digits = 3)

# --- Figuras ----------------------------------------------------------------

tema <- ggpubr::theme_pubr(base_size = 14) + theme(strip.background = element_blank(), legend.position = "top")

# p1: coeficiente ajustado vs. n unidos (log10), grupos destacados, por replica
plot_cob <- function(df, y, ylab, linea = TRUE) {
  dg <- df %>% filter(grupo != "resto")
  ggplot(df, aes(n_unidos, .data[[y]])) +
    geom_hline(yintercept = 0, linetype = "dotted", col = "grey40") +
    geom_point(data = ~ filter(.x, grupo == "resto"), col = "grey70", size = 0.9, alpha = 0.6) +
    { if (linea) geom_line(data = ~ arrange(.x, n_unidos), aes(y = ajuste_cobertura), col = "grey20", linewidth = 0.7) } +
    geom_point(data = dg, aes(col = grupo), size = 2.6) +
    geom_text_repel(data = dg, aes(label = TF, col = grupo), size = 3.4, min.segment.length = 0, max.overlaps = 50, show.legend = FALSE) +
    scale_x_log10() +
    scale_color_manual(values = grupo_clr, labels = grupo_lab) +
    facet_wrap(~rep) +
    tema +
    labs(x = "Promotores unidos (ReMap, escala log10)", y = ylab, col = NULL)
}
p_cob <- plot_cob(cob, "estimate", "Efecto ajustado de la unión sobre\nel rango de varianza (1-100)") +
  labs(caption = "Línea: loess robusto sobre log10(n unidos).")
save_jpg_pdf(p_cob, "p1_coef_vs_cobertura", 13, 6.5)
p_res <- plot_cob(cob, "residuo", "Residuo del efecto ajustado\nrespecto de la cobertura", linea = FALSE)
save_jpg_pdf(p_res, "p1_residuo_vs_cobertura", 13, 6.5)

# p2: distribucion del residuo y del percentil por grupo
p_grp <- cob %>%
  select(TF, rep, grupo, Residuo = residuo, `Percentil en la ventana ±25%` = percentil_ventana) %>%
  pivot_longer(c(Residuo, `Percentil en la ventana ±25%`), names_to = "variable") %>%
  mutate(grupo = factor(grupo, levels = names(grupo_lab))) %>%
  ggplot(aes(grupo, value, col = grupo)) +
  geom_boxplot(data = ~ filter(.x, grupo == "resto"), outlier.shape = NA, col = "grey50") +
  geom_jitter(data = ~ filter(.x, grupo == "resto"), width = 0.2, size = 0.5, alpha = 0.3) +
  geom_point(data = ~ filter(.x, grupo != "resto"), size = 2.4) +
  geom_text_repel(data = ~ filter(.x, grupo != "resto"), aes(label = TF), size = 3, direction = "y", nudge_x = 0.3, min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual(values = grupo_clr, labels = grupo_lab) +
  scale_x_discrete(labels = c(metilacion = "Metilación", cpg_no_metilado = "CpG no\nmetilado", resto = "Resto")) +
  facet_grid(variable ~ rep, scales = "free_y", switch = "y") +
  tema +
  theme(legend.position = "none", strip.placement = "outside") +
  labs(x = NULL, y = NULL)
save_jpg_pdf(p_grp, "p2_residuo_por_grupo", 11, 9)

# p3: NES de los terminos mas relevantes en los tres rankings
top_res <- gsea_res_df %>% group_by(ID) %>% summarise(pmin = min(p.adjust)) %>% slice_min(pmin, n = 15, with_ties = FALSE) %>% pull(ID)
p_nes <- gsea_all %>%
  filter(ID %in% union(terminos_ref, top_res)) %>%
  mutate(
    ranking = factor(ranking, levels = c("AUC_original", "ajustado", "residuo_cobertura"), labels = c("AUC original", "Coeficiente ajustado", "Residuo por cobertura")),
    Description = fct_reorder(paste0(Description, " (", ID, ")"), NES),
    sig = p.adjust < 0.1
  ) %>%
  ggplot(aes(NES, Description, col = ranking, shape = sig)) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_point(size = 3, position = position_dodge(0.6)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p ajustado < 0.1", `FALSE` = "p ajustado ≥ 0.1")) +
  scale_color_manual(values = c("black", "#AD343E", "#1B8C8E")) +
  facet_wrap(~rep) +
  tema +
  theme(legend.direction = "vertical") +
  labs(x = "NES (negativo = factores asociados a menor ruido)", y = NULL, col = "Ranking", shape = NULL)
save_jpg_pdf(p_nes, "p3_gsea_nes_comparacion", 15, 9)

# p4: coeficiente vs. cobertura dentro de cada estrato
p_estr <- cob_cgi %>%
  filter(modelo %in% c("estrato_con_isla", "estrato_sin_isla")) %>%
  mutate(modelo = factor(modelos_cgi[as.character(modelo)], levels = modelos_cgi)) %>%
  plot_cob("estimate", "Efecto ajustado de la unión dentro del estrato") +
  facet_grid(modelo ~ rep) +
  labs(caption = "Línea: loess robusto sobre log10(n unidos) dentro de cada estrato y réplica.")
save_jpg_pdf(p_estr, "p4_coef_vs_cobertura_estratos", 13, 11)

# p4: percentil de ventana de cada miembro en las cuatro versiones
orden_tf <- cob %>% filter(grupo != "resto") %>% group_by(TF) %>% summarise(m = mean(percentil_ventana)) %>% arrange(m) %>% pull(TF)
p_pct <- cob_cgi %>%
  filter(grupo != "resto") %>%
  mutate(modelo = factor(modelos_cgi[as.character(modelo)], levels = modelos_cgi), TF = factor(TF, levels = orden_tf)) %>%
  ggplot(aes(percentil_ventana, TF, col = modelo)) +
  geom_vline(xintercept = 0.5, col = "grey60") +
  geom_point(size = 2.6, position = position_dodge(0.6)) +
  scale_color_manual(values = c("black", "#D6741F", "#1B8C8E", "#AD343E")) +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(~rep) +
  tema +
  theme(legend.direction = "vertical", axis.text.y = element_text(color = ifelse(orden_tf %in% grupos$metilacion, grupo_clr["metilacion"], grupo_clr["cpg_no_metilado"]))) +
  labs(x = "Percentil en la ventana de cobertura ±25% (0 = mayor reducción del ruido)", y = NULL, col = NULL,
    caption = "Rojo: asociados a metilación. Verde azulado: asociados a CpG no metilado. Sin punto = factor no evaluable en el estrato.")
save_jpg_pdf(p_pct, "p4_percentil_ventana_miembros", 12, 8)

# p6: efecto por categoria de union DNMT3B / ZBTB33
p_cat <- efecto_cat %>%
  mutate(categoria = fct_inorder(categoria)) %>%
  ggplot(aes(estimate, fct_rev(categoria))) +
  geom_vline(xintercept = 0, col = "grey60") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0) +
  geom_point(size = 2.6, col = "#AD343E") +
  geom_text(aes(x = Inf, label = paste0("n = ", n)), hjust = 1.1, size = 3.2, col = "grey30") +
  facet_wrap(~rep) +
  tema +
  labs(x = "Efecto ajustado sobre el rango de varianza (IC 95%)", y = NULL,
    caption = "Categorías: frente a los promotores sin ninguno de los dos. Abajo: efecto de cada factor excluyendo los promotores unidos por el otro.")
save_jpg_pdf(p_cat, "p6_efecto_por_categoria", 12, 5)

message("Exploratorio guardado en ", out_dir)
