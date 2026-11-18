# silentes_ltr (ver docs/mapping_figuras.csv para el numero de figura vigente)
# Promotores silentes en HEK293 (hek_tpm = 0) con actividad en el reportero. El
# analisis de discordancia (discordancia_features, discordancia_remap_ora) los
# excluye porque sin senal endogena el rango endogeno queda empatado.
#
# Grupos, por promotor: estado en HEK293 (activo: hek_tpm > 0; silente: = 0) x
# actividad en el reportero (tercil superior o inferior de la actividad media,
# terciles por replica sobre todos los promotores de alta confianza, en todas las
# replicas en que se midio). Rango en el reportero = rango de la actividad media
# dentro de cada replica / n, promediado entre replicas.
#
# Paneles:
# (A) Esquema de los cuatro grupos con su numero de promotores.
# (B) Efecto de cada caracteristica de secuencia (grupo "seq" de
#     build_tidy_features, mas LTR y transposones de ADN) sobre la actividad alta
#     en el reportero, en silentes (alto vs. bajo) contra activos (alto vs. bajo).
#     Regresion logistica por variable, una observacion por promotor, ajustada por
#     isla CpG (CGI_frag) y G+C. Interaccion = diferencia de log(OR) entre silentes
#     y activos, test z, BH por catalogo. LTR es la que mas se aparta de la diagonal.
# (C) Lo mismo para cada motivo de JASPAR (ajuste ademas por CpG o/e) y la union de
#     cada TF de ReMap (ademas por log(1 + numero de TF unidos)).
# (D) Rango en el reportero de los silentes segun la clase de elemento
#     transponible (TE_superclass de prom_df, superposicion > 10 pb) y, para LTR,
#     la familia del LTR de mayor superposicion (RepeatMasker). Wilcoxon contra los
#     silentes sin TE, BH.
# (E) Silentes con LTR segun si el TSS anotado (posicion 237 del fragmento de 252
#     pb) cae dentro o fuera del LTR. Wilcoxon.
# (F) Metilacion de CpG en HEK293 (WGBS de GEO, 5 muestras control de GSE243785 y
#     GSE233830; R/00_prom_features/build_metilacion_wgbs_hek293.R) en el fragmento de la library, para
#     silentes y activos (todos los de alta confianza, sin separar por reportero),
#     segun isla CpG y LTR. Lecturas metiladas / totales sumadas entre muestras;
#     promotor medido si alguna muestra tiene >= 3 CpG con >= 5 lecturas.
#
# Motivos: JASPAR 2024 CORE vertebrados no redundante, ambas hebras, puntaje
# relativo >= 0.85 (como motivos_tf), cacheados en
# data/processed/silentes_jaspar_presencia_v1.rds (lo genera tambien
# R/64_exploratorio_silentes_reportero/silentes_reportero.R, que tiene el resto
# del analisis exploratorio).
#
# Tablas: tables/silentes_ltr_interaccion.tsv (caracteristicas de secuencia,
# motivos y ReMap) y tables/silentes_ltr_clases.tsv.
# Requiere: data/processed/{activity_stats_highconf,prom_df,remap_tf_hits,
# metilacion_hek293_wgbs_v1}.tsv,
# data/raw/library.bed, RepeatMasker (path_repeatmasker) y, si no esta el cache
# de motivos, el fasta de la library y JASPAR. Run from the TesisDoc repo root.

library(tidyverse)
library(patchwork)
library(ggrepel)
library(GenomicRanges)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")
source("R/00_prom_features/heavy_data_paths.R")

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
slice <- dplyr::slice

slug <- "silentes_ltr"
out_dir <- fig_dir(slug)
cache_motivos <- "data/processed/silentes_jaspar_presencia_v1.rds"
path_fa <- "../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa"
path_jaspar <- "data/external/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
min_con <- 50
min_celda <- 5
clr_silente <- "#C9743A"

save_panel <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot, width = width, height = height, units = "in", device = if (ext == "pdf") cairo_pdf else NULL)
  }
}
th <- theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey90", color = NA), strip.text = element_text(face = "bold"))
fmt_q <- function(q) if_else(q < 0.001, "FDR < 0.001", paste0("FDR = ", formatC(q, format = "fg", digits = 2)))

# --- Datos y grupos -------------------------------------------------------------

stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  group_by(rep) %>%
  mutate(terc = ntile(mean, 3), rr = dense_rank(mean) / n()) %>%
  ungroup()
prom <- data %>%
  group_by(seq_id) %>%
  summarise(n_rep = n(), siempre_alto = all(terc == 3), siempre_bajo = all(terc == 1), rr_medio = mean(rr), .groups = "drop") %>%
  inner_join(distinct(select(data, -rep, -mean, -terc, -rr), seq_id, .keep_all = TRUE), by = "seq_id") %>%
  mutate(
    estado = if_else(hek_tpm > 0, "Activo", "Silente"),
    reportero = case_when(siempre_alto ~ "Alto", siempre_bajo ~ "Bajo"),
    cgi = as.numeric(CGI_frag)
  )
n_grupos <- prom %>% filter(!is.na(reportero)) %>% count(estado, reportero)
message("Promotores: ", nrow(prom), "; silentes: ", sum(prom$estado == "Silente"))
print(n_grupos)

# --- (A) Esquema ---------------------------------------------------------------

p_esquema <- n_grupos %>%
  mutate(
    estado = factor(paste(estado, "en HEK293"), levels = c("Activo en HEK293", "Silente en HEK293")),
    reportero = factor(reportero, levels = c("Bajo", "Alto")),
    destacado = estado == "Silente en HEK293" & reportero == "Alto"
  ) %>%
  ggplot(aes(reportero, estado)) +
  geom_tile(aes(fill = destacado), color = "white", linewidth = 2) +
  geom_text(aes(label = n, color = destacado), size = 7, fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = clr_silente, `FALSE` = "grey88"), guide = "none") +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = "Actividad en el reportero (tercil)", y = NULL) +
  th +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.ticks = element_blank())

# --- (B) Silentes vs. activos ---------------------------------------------------

ajustar <- function(df, x, covs) {
  ok <- !is.na(x)
  y <- df$y[ok]
  x <- x[ok]
  if (sum(x) < min_con || sum(!x) < min_con || sum(x & y == 1) < min_celda || sum(x & y == 0) < min_celda) {
    return(NULL)
  }
  mm <- cbind(1, x = as.numeric(x), as.matrix(df[ok, covs, drop = FALSE]))
  fit <- suppressWarnings(glm.fit(mm, y, family = binomial()))
  if (!fit$converged || fit$rank < ncol(mm)) {
    return(NULL)
  }
  se <- sqrt(diag(chol2inv(fit$qr$qr[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE])))
  s <- se[match("x", colnames(mm)[fit$qr$pivot])]
  b <- unname(fit$coefficients["x"])
  tibble(n_con = sum(x), pct_alto = 100 * mean(x[y == 1]), pct_bajo = 100 * mean(x[y == 0]), b = b, se = s, p_value = 2 * pnorm(-abs(b / s)))
}
datos_estado <- function(e) prom %>% filter(estado == e, !is.na(reportero), !is.na(CGI_frag)) %>% mutate(y = as.numeric(reportero == "Alto"))
barrer <- function(variables, covs, catalogo) {
  map(c("Silente", "Activo"), function(e) {
    df <- datos_estado(e)
    imap(variables, function(ids, nombre) {
      r <- ajustar(df, df$seq_id %in% ids, covs)
      if (is.null(r)) NULL else mutate(r, variable = nombre)
    }) %>%
      list_rbind() %>%
      mutate(estado = e, p_adj = p.adjust(p_value, "BH"))
  }) %>%
    list_rbind() %>%
    mutate(catalogo = catalogo)
}

# motivos
if (file.exists(cache_motivos)) {
  presencia <- readRDS(cache_motivos)
} else {
  library(Biostrings)
  library(universalmotif)
  fa <- readDNAStringSet(path_fa)[prom$seq_id]
  jaspar <- read_jaspar(path_jaspar)
  nombres <- make.unique(map_chr(jaspar, ~ .x@altname))
  presencia <- map2(jaspar, nombres, function(pcm, nombre) {
    pwm <- convert_type(pcm, "PWM")@motif
    s_min <- sum(apply(pwm, 2, min))
    s_max <- sum(apply(pwm, 2, max))
    r <- scan_sequences(pcm, fa, threshold = s_min + 0.85 * (s_max - s_min), threshold.type = "logodds.abs", RC = TRUE, calc.pvals = FALSE, warn.NA = FALSE, verbose = 0, nthreads = 12)
    tibble(motivo = nombre, seq_id = unique(as.character(as.data.frame(r)$sequence)))
  }) %>% list_rbind()
  saveRDS(presencia, cache_motivos)
}

# ReMap, todas las lineas, sin marcas de histona
remap <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE) %>%
  filter(seq_id %in% prom$seq_id, !str_detect(TF, "^H[0-9]")) %>%
  distinct(seq_id, TF)
prom <- prom %>%
  left_join(count(remap, seq_id, name = "n_tf"), by = "seq_id") %>%
  mutate(log_tf = log1p(replace_na(n_tf, 0)))

# Caracteristicas de secuencia: las del grupo "seq" de build_tidy_features (sin
# "Elementos transponibles", que agrupa a las clases), mas LTR y transposones de
# ADN. Islas CpG se ajusta solo por G+C y alto G+C solo por isla CpG. LTR tiene
# pocos promotores entre los activos, por eso el minimo baja a 15 en este catalogo.
feats <- prom %>%
  mutate(rep = "todos") %>%
  build_tidy_features(keep = "seq_id") %>%
  mutate(across(-c(seq_id, rep), ~ .x == "TRUE")) %>%
  select(seq_id, all_of(setdiff(feature_groups$feature[feature_groups$group == "seq"], "Elementos transponibles"))) %>%
  left_join(prom %>% transmute(seq_id, LTR = TE_superclass %in% "LTR", `Transposón de ADN` = TE_superclass %in% "DNA"), by = "seq_id")
barrer_caracteristicas <- function() {
  min_con <<- 15
  on.exit(min_con <<- 50)
  map(c("Silente", "Activo"), function(e) {
    df <- datos_estado(e) %>% left_join(feats, by = "seq_id")
    map(setdiff(names(feats), "seq_id"), function(f) {
      covs <- switch(f, "Islas CpG" = "g_c", "Alto contenido G+C" = "cgi", c("cgi", "g_c"))
      r <- ajustar(df, df[[f]], covs)
      if (is.null(r)) NULL else mutate(r, variable = f)
    }) %>%
      list_rbind() %>%
      mutate(estado = e, p_adj = p.adjust(p_value, "BH"))
  }) %>%
    list_rbind() %>%
    mutate(catalogo = "Características de secuencia")
}
res <- bind_rows(
  barrer_caracteristicas(),
  barrer(split(presencia$seq_id, presencia$motivo), c("cgi", "g_c", "cpg_oe"), "Motivos de JASPAR"),
  barrer(split(remap$seq_id, remap$TF), c("cgi", "g_c", "log_tf"), "Unión de TF (ReMap)")
)
inter <- res %>%
  select(catalogo, variable, estado, n_con, pct_alto, pct_bajo, b, se, p_adj) %>%
  pivot_wider(names_from = estado, values_from = c(n_con, pct_alto, pct_bajo, b, se, p_adj)) %>%
  mutate(
    OR_silentes = exp(b_Silente), OR_activos = exp(b_Activo),
    p_inter = 2 * pnorm(-abs((b_Silente - b_Activo) / sqrt(se_Silente^2 + se_Activo^2)))
  ) %>%
  group_by(catalogo) %>%
  mutate(p_adj_inter = p.adjust(p_inter, "BH")) %>%
  ungroup()
write_tsv(
  inter %>% select(catalogo, variable, starts_with("n_con"), starts_with("pct_"), OR_silentes, p_adj_silentes = p_adj_Silente, OR_activos, p_adj_activos = p_adj_Activo, p_inter, p_adj_inter),
  "tables/silentes_ltr_interaccion.tsv"
)

d_b <- inter %>%
  filter(!is.na(b_Silente), !is.na(b_Activo)) %>%
  mutate(
    sig = case_when(
      p_adj_inter < 0.05 ~ "Interacción (FDR < 0.05)",
      p_adj_Silente < 0.05 & p_adj_Activo < 0.05 ~ "FDR < 0.05 en ambos",
      p_adj_Silente < 0.05 | p_adj_Activo < 0.05 ~ "FDR < 0.05 en uno",
      TRUE ~ "No significativo"
    ),
    sig = factor(sig, levels = c("FDR < 0.05 en ambos", "FDR < 0.05 en uno", "No significativo", "Interacción (FDR < 0.05)")),
    etiqueta_var = str_remove(variable, "\\.[0-9]+$") # sufijo de make.unique en JASPAR
  ) %>%
  arrange(desc(sig))
resumen_b <- d_b %>%
  group_by(catalogo) %>%
  summarise(n = n(), rho = cor(b_Silente, b_Activo, method = "spearman"), n_inter = sum(p_adj_inter < 0.05), .groups = "drop") %>%
  mutate(etiqueta = paste0("n = ", n, "\nSpearman = ", formatC(rho, format = "f", digits = 2), "\nInteracción FDR < 0.05: ", n_inter))
print(resumen_b)
print(inter %>% filter(p_adj_inter < 0.05) %>% select(catalogo, variable, OR_silentes, OR_activos, p_adj_inter))
print(inter %>% filter(catalogo == "Características de secuencia") %>% select(variable, n_con_Silente, n_con_Activo, OR_silentes, p_adj_Silente, OR_activos, p_adj_Activo, p_adj_inter) %>% arrange(p_adj_inter), n = 30)

colores_sig <- c("FDR < 0.05 en ambos" = thesis_clr, "FDR < 0.05 en uno" = "#9CC6D6", "No significativo" = "grey80", "Interacción (FDR < 0.05)" = clr_silente)
# el resumen va arriba a la izquierda o, si ahi hay puntos (LTR en B), abajo a la derecha
dispersion <- function(d, etiquetas, tam_punto, lim, resumen_abajo = FALSE) {
  pos <- if (resumen_abajo) list(x = lim[2], y = lim[1], h = 1, v = 0) else list(x = lim[1], y = lim[2], h = 0, v = 1)
  ggplot(d, aes(b_Activo / log(2), b_Silente / log(2))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, color = "grey80") +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_point(aes(color = sig), size = tam_punto, alpha = 0.85) +
    geom_text_repel(data = etiquetas, aes(label = etiqueta_var), size = 3.6, min.segment.length = 0, segment.size = 0.25, max.overlaps = Inf, seed = 1) +
    geom_text(data = semi_join(resumen_b, d, by = "catalogo"), aes(label = etiqueta), x = pos$x, y = pos$y, inherit.aes = FALSE, hjust = pos$h, vjust = pos$v, size = 3.6, lineheight = 0.9) +
    facet_wrap(~catalogo) +
    coord_equal(xlim = lim, ylim = lim) +
    scale_color_manual(values = colores_sig, name = NULL, drop = FALSE) +
    labs(x = "log2(OR) en activos:\nreportero alto vs. bajo", y = "log2(OR) en silentes:\nreportero alto vs. bajo") +
    th +
    theme(legend.position = "none")
}

# (B) caracteristicas de secuencia
d_feat <- d_b %>% filter(catalogo == "Características de secuencia")
etiq_feat <- d_feat %>% filter(p_adj_inter < 0.05 | variable %in% c("CCAAT", "TATA-box", "GC-box"))
lim_feat <- range(c(d_feat$b_Silente, d_feat$b_Activo, 0)) / log(2) + c(-0.3, 0.3)
p_caracteristicas <- dispersion(d_feat, etiq_feat, 2.6, lim_feat, resumen_abajo = TRUE)

# (C) motivos y union de TF
d_mt <- d_b %>% filter(catalogo != "Características de secuencia")
etiq_mt <- d_mt %>%
  filter(sig %in% c("FDR < 0.05 en ambos", "Interacción (FDR < 0.05)")) %>%
  group_by(catalogo) %>%
  slice_max(pmin(-log10(p_adj_Silente), -log10(p_adj_Activo)), n = 6, with_ties = FALSE) %>%
  ungroup()
p_silentes_activos <- dispersion(d_mt, etiq_mt, 1.3, range(c(d_mt$b_Silente, d_mt$b_Activo)) / log(2))
leyenda_sig <- guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)))

# --- LTR: familia y posicion del TSS ---------------------------------------------

silentes <- prom %>% filter(estado == "Silente")
sin_te <- silentes %>% filter(is.na(TE_superclass))
lib <- read_tsv("data/raw/library.bed", col_names = c("chr", "ini0", "fin", "seq_id", "score", "hebra"), show_col_types = FALSE) %>%
  filter(seq_id %in% silentes$seq_id[silentes$TE_superclass %in% "LTR"]) %>%
  mutate(tss = if_else(hebra == "+", ini0 + 237L, fin - 236L))
rm_ltr <- data.table::fread(path_repeatmasker, sep = " ", fill = TRUE, select = c("seqname", "start", "end", "repeat_class")) %>%
  as_tibble() %>%
  filter(str_starts(repeat_class, "LTR"))
ov <- findOverlaps(GRanges(lib$chr, IRanges(lib$ini0 + 1L, lib$fin)), GRanges(rm_ltr$seqname, IRanges(rm_ltr$start, rm_ltr$end)), minoverlap = 11)
ltr_hits <- bind_cols(lib[queryHits(ov), c("seq_id", "tss", "ini0", "fin")], rm_ltr[subjectHits(ov), ]) %>%
  mutate(solape = pmin(fin, end) - pmax(ini0 + 1L, start) + 1L, tss_dentro = tss >= start & tss <= end, familia = str_remove(repeat_class, "^LTR/?")) %>%
  group_by(seq_id) %>%
  slice_max(solape, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(seq_id, familia, tss_dentro)
stopifnot(nrow(ltr_hits) == sum(silentes$TE_superclass %in% "LTR"))

# (C) clases
niveles_clase <- c("Sin TE", "SINE", "LINE", "DNA", "LTR: ERVL", "LTR: ERV1", "LTR: ERVL-MaLR")
d_c <- silentes %>%
  left_join(ltr_hits, by = "seq_id") %>%
  mutate(clase = case_when(
    is.na(TE_superclass) ~ "Sin TE",
    TE_superclass == "LTR" ~ paste0("LTR: ", familia),
    TRUE ~ TE_superclass
  )) %>%
  filter(clase %in% niveles_clase)
tabla_c <- d_c %>%
  group_by(clase) %>%
  summarise(n = n(), rango_medio = mean(rr_medio), pct_alto = 100 * mean(siempre_alto), .groups = "drop") %>%
  mutate(p_value = map_dbl(clase, ~ if (.x == "Sin TE") NA_real_ else wilcox.test(d_c$rr_medio[d_c$clase == .x], sin_te$rr_medio)$p.value)) %>%
  mutate(p_adj = p.adjust(p_value, "BH"))
write_tsv(tabla_c, "tables/silentes_ltr_clases.tsv")
print(tabla_c)
d_c <- d_c %>%
  left_join(select(tabla_c, clase, n, p_adj), by = "clase") %>%
  mutate(
    etiqueta = factor(paste0(clase, " (", n, ")"), levels = paste0(niveles_clase, " (", tabla_c$n[match(niveles_clase, tabla_c$clase)], ")")),
    es_ltr = str_starts(clase, "LTR")
  )
ref_rango <- median(sin_te$rr_medio)
p_clases <- ggplot(d_c, aes(rr_medio, etiqueta, fill = es_ltr)) +
  geom_vline(xintercept = ref_rango, linetype = "dashed", color = "grey50") +
  geom_violin(color = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_text(data = distinct(d_c, etiqueta, p_adj) %>% filter(!is.na(p_adj)), aes(x = 1.03, y = etiqueta, label = fmt_q(p_adj)), inherit.aes = FALSE, hjust = 0, size = 3.8) +
  scale_fill_manual(values = c(`TRUE` = clr_silente, `FALSE` = "grey80"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.38)) +
  labs(x = "Rango en el reportero", y = NULL) +
  th

# (D) TSS dentro o fuera del LTR
d_d <- silentes %>%
  inner_join(ltr_hits, by = "seq_id") %>%
  mutate(grupo = if_else(tss_dentro, "TSS dentro\ndel LTR", "TSS fuera\ndel LTR")) %>%
  add_count(grupo) %>%
  mutate(grupo = factor(paste0(grupo, "\n(", n, ")")))
p_tss <- wilcox.test(rr_medio ~ tss_dentro, data = d_d)$p.value
print(d_d %>% group_by(grupo) %>% summarise(n = n(), rango_medio = mean(rr_medio), pct_alto = 100 * mean(siempre_alto)))
message("TSS dentro vs. fuera, Wilcoxon p = ", signif(p_tss, 2))
p_tss_ltr <- ggplot(d_d, aes(grupo, rr_medio, fill = tss_dentro)) +
  geom_hline(yintercept = ref_rango, linetype = "dashed", color = "grey50") +
  geom_violin(color = NA, scale = "width", width = 0.8, alpha = 0.85) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.6, color = "grey25") +
  annotate("text", x = 1.5, y = 1.06, label = paste0("p = ", formatC(p_tss, format = "e", digits = 1)), size = 3.8) +
  scale_fill_manual(values = c(`TRUE` = clr_silente, `FALSE` = "grey80"), guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.1)) +
  labs(x = NULL, y = "Rango en el reportero") +
  th

# --- (F) Metilacion en HEK293, silentes vs. activos -------------------------------

clr_activo <- thesis_clr
metilacion <- read_tsv("data/processed/metilacion_hek293_wgbs_v1.tsv", show_col_types = FALSE) %>%
  filter(linea == "HEK293", ventana == "Fragmento") %>%
  mutate(met = metilacion * lecturas / 100) %>%
  group_by(seq_id) %>%
  summarise(metilacion = 100 * sum(met) / sum(lecturas), n_cpg5 = max(n_cpg5), .groups = "drop") %>%
  filter(n_cpg5 >= 3)
niveles_sub <- c("Con isla CpG", "Sin isla CpG,\nsin LTR", "LTR")
d_f <- prom %>%
  filter(!is.na(CGI_frag)) %>%
  mutate(sub = factor(case_when(TE_superclass %in% "LTR" ~ "LTR", CGI_frag ~ "Con isla CpG", TRUE ~ "Sin isla CpG,\nsin LTR"), levels = niveles_sub)) %>%
  left_join(metilacion, by = "seq_id")
tabla_f <- d_f %>%
  group_by(sub, estado) %>%
  summarise(n = n(), n_medido = sum(!is.na(metilacion)), mediana = median(metilacion, na.rm = TRUE), pct_metilado = 100 * mean(metilacion >= 50, na.rm = TRUE), .groups = "drop")
print(tabla_f)
write_tsv(tabla_f, "tables/silentes_ltr_metilacion.tsv")
d_f <- filter(d_f, !is.na(metilacion))
p_metilacion <- ggplot(d_f, aes(sub, metilacion, fill = estado)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(data = tabla_f, aes(sub, -9, label = n_medido, group = estado), position = position_dodge(width = 0.8), size = 3.4, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Silente = clr_silente, Activo = clr_activo), breaks = c("Silente", "Activo"), labels = c("Silentes en HEK293", "Activos en HEK293"), name = NULL) +
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Metilación de CpG en HEK293 (%)") +
  th +
  theme(legend.position = "top")

# --- Figura completa y paneles sueltos -------------------------------------------

save_panel(p_esquema, "silentes_ltr_grupos", 5, 3.5)
p_caracteristicas <- p_caracteristicas + theme(legend.position = "top") + leyenda_sig
p_silentes_activos <- p_silentes_activos + theme(legend.position = "top") + leyenda_sig
save_panel(p_esquema, "silentes_ltr_grupos", 5, 3.5)
save_panel(p_caracteristicas, "silentes_ltr_caracteristicas", 6.5, 7)
save_panel(p_silentes_activos, "silentes_ltr_silentes_vs_activos", 10, 6)
save_panel(p_clases, "silentes_ltr_clases", 7.5, 5)
save_panel(p_tss_ltr, "silentes_ltr_tss", 4.5, 5)
save_panel(p_metilacion, "silentes_ltr_metilacion", 6, 5)

# B y C comparten la escala de color: en la figura completa la leyenda va solo en B
fila_1 <- (p_esquema | p_caracteristicas | (p_silentes_activos + theme(legend.position = "none"))) + plot_layout(widths = c(0.75, 1, 2))
fila_2 <- (p_clases | p_tss_ltr | p_metilacion) + plot_layout(widths = c(1.4, 0.7, 1.1))
figura <- (fila_1 / fila_2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))
save_panel(figura, slug, 18, 11)

message("Figuras guardadas en ", out_dir)
