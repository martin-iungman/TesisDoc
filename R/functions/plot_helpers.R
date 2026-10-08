# Shared ggplot helpers/theme for thesis figures (Spanish labels, ggpubr style).

library(ggplot2)
library(ggpubr)

thesis_clr <- "#358AAA"
thesis_clr_bg <- "#dfedf6"
thesis_clr_label <- "#AADAD4"

# Barplot of counts against a total: a full-height background bar (the
# total), a foreground bar (the count) and a value label on top. Used for
# "N of total promoters have this feature" figures (EPD motifs, CGI, TSS
# motifs, conservation classes...).
count_bar_labeled <- function(data, x, n, total, base_size = 16,
                               label_nudge_frac = 0.07, coord_flip = FALSE) {
  x <- rlang::enquo(x)
  n <- rlang::enquo(n)
  total <- rlang::enquo(total)
  # `total` may be a scalar from the caller's env (e.g. n_promoters) or a
  # column in `data` (e.g. total_n); eval_tidy resolves either.
  total_val <- max(rlang::eval_tidy(total, data), na.rm = TRUE)
  p <- ggplot(data, aes(x = !!x, label = !!n)) +
    geom_col(aes(y = !!total), fill = thesis_clr_bg) +
    geom_col(aes(y = !!n), fill = thesis_clr) +
    geom_label(aes(y = !!n), nudge_y = total_val * label_nudge_frac, fill = thesis_clr_label, size = 4) +
    theme_pubclean(base_size = base_size)
  if (coord_flip) p <- p + coord_flip()
  p
}

# Bins de 100 promotores segun rango de actividad media, por replica (se
# descarta el resto de la division para que los bins queden parejos). Usado
# para los scatters de union de TF/motivo vs. actividad (R/16_tf_contexto_
# endogeno, R/17_remap_activity_gsea).
add_mean_sw_bins <- function(df) {
  excess <- df %>% group_by(rep) %>% summarise(excess = n() %% 100, .groups = "drop")
  df %>%
    group_split(rep) %>%
    purrr::map2(excess$excess, ~ .x %>%
      mutate(mean_rank_sw = row_number(mean)) %>%
      filter(mean_rank_sw > .y) %>%
      mutate(mean_rank_sw = row_number(mean), mean_sw = ceiling(mean_rank_sw / 100)) %>%
      ungroup()) %>%
    purrr::list_rbind()
}

# Proporcion de promotores con `feature`==TRUE por bin de actividad
# (mean_sw, ver add_mean_sw_bins), con smooth loess por replica.
tf_scatter <- function(df, feature, titulo) {
  df %>%
    group_by(rep, mean_sw) %>%
    summarise(prop = sum(.data[[feature]]) / n(), .groups = "drop") %>%
    ggplot(aes(mean_sw, prop)) +
    geom_point(col = thesis_clr) +
    geom_smooth(col = "#216869") +
    facet_wrap(~rep) +
    ylim(0, 1) +
    labs(x = "Actividad media (bins de 100, ordenados por rango)", y = "Proporción de promotores", title = titulo) +
    theme_bw(base_size = 16)
}

# Actividad media de un promotor endogeno (CAGE TPM, FANTOM5) por bin de
# actividad del reportero (mean_sw, ver add_mean_sw_bins), con smooth
# lineal por replica y escala log10 en Y. Usado por R4 (HEK293/HeLa/
# musculo esqueletico - ver R/00_prom_features/build_endo_cage_activity.R
# y build_prom_features.R para hek_tpm).
endo_activity_scatter <- function(data, tpm_col, titulo) {
  data %>%
    group_by(rep, mean_sw) %>%
    summarise(tpm = mean(.data[[tpm_col]], na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(mean_sw, tpm)) +
    geom_point(col = "#AD343E") +
    geom_smooth(col = "#216869", method = "lm") +
    facet_wrap(~rep) +
    scale_y_log10() +
    labs(x = "Bins de actividad del promotor reportero", y = "Actividad media del promotor\nendógeno (CAGE TPM)", title = titulo) +
    theme_bw(base_size = 16)
}

# Rango de varianza (ruido) dentro de cada bin de actividad media
# (mean_sw, ver add_mean_sw_bins), por replica. add_noise_rank() hace
# ambos pasos juntos - usado por todo R3 (ruido) para construir
# var_rank_sw desde cero, incluida la re-binning de subconjuntos (ej.
# solo promotores no-CGI en R3.4).
add_var_rank_sw <- function(df) {
  df %>%
    group_by(rep, mean_sw) %>%
    mutate(var_rank_sw = row_number(var)) %>%
    ungroup()
}
add_noise_rank <- function(df) {
  add_mean_sw_bins(df) %>% add_var_rank_sw()
}

# Descompone el par (rank_reporter, rank_endo) en "nivel" (avg_rank,
# promedio de ambos rangos normalizados 0-1 por replica - mejor proxy de
# actividad real que cualquiera de los dos solo, promedia el ruido de
# cada assay) y "discordancia con signo" (dif_signed = rr - re; positivo
# = el reportero sobreestima respecto a la actividad endogena, negativo
# = el reportero subestima - convencion elegida por el autor 2026-07-31,
# antes era re - rr con el signo invertido).
# avg_rank y dif_signed son ~ortogonales por construccion (cor~-0.02 en
# estos datos) - analisis tipo Bland-Altman (mean-difference), evita el
# confounding de controlar una feature/TF por uno solo de los dos rangos
# (rank_reporter unicamente) al modelar la discordancia. Requiere `mean`
# (actividad del reportero) y `max_tpm` (actividad endogena maxima,
# FANTOM5 - ver R/00_prom_features/build_fantom_endo_activity_summary.R)
# ya presentes en df. Usado por R4.6/R4.7 (rank_dif_lmm y derivados) -
# ver R/27_rank_dif_lmm/rank_dif_lmm.R para la discusion completa de por
# que reemplaza al diseño anterior (controlar solo por rank_reporter).
add_rank_discordance <- function(df) {
  df %>%
    group_by(rep) %>%
    mutate(
      max_tpm = replace_na(max_tpm, 0),
      rr = dense_rank(mean) / n(),
      re = row_number(max_tpm) / n(),
      avg_rank = (rr + re) / 2,
      dif_signed = rr - re
    ) %>%
    ungroup()
}

# ROC curve (por umbral de rango) de una feature booleana contra
# var_rank_sw (rango de varianza dentro de cada bin de actividad media -
# ver add_mean_sw_bins), para testear si una feature de secuencia predice
# ruido inusualmente alto/bajo. Ported from transcriptional_library/
# Analysis/scripts/noise_analysis.qmd.
roc_curve <- function(data, feature) {
  purrr::map(1:99, ~ data %>%
    count(gr = var_rank_sw > .x, across(all_of(feature))) %>%
    mutate(threshold = .x)) %>%
    purrr::list_rbind() %>%
    pivot_wider(names_from = c(gr, all_of(feature)), values_from = n) %>%
    mutate(TPR = TRUE_TRUE / (TRUE_TRUE + FALSE_TRUE), FPR = TRUE_FALSE / (FALSE_FALSE + TRUE_FALSE)) %>%
    arrange(FPR)
}

# Area bajo la curva (regla del trapecio) de un roc_curve().
noise_auc <- function(roc) {
  total <- 0
  for (i in seq_len(nrow(roc) - 1)) {
    deltax <- roc$FPR[i + 1] - roc$FPR[i]
    deltay <- roc$TPR[i + 1] - roc$TPR[i]
    total <- sum(total, deltax * roc$TPR[i] + deltax * deltay / 2, na.rm = TRUE)
  }
  total
}

# Scatter de ruido (var_rank_sw) vs. bin de actividad media, coloreado por
# una feature booleana, con smooth por nivel (paleta separada para puntos
# vs. curvas via ggnewscale, igual que el original).
noise_scatter <- function(data, feature, titulo, leyenda) {
  data %>%
    ggplot(aes(mean_sw, var_rank_sw, col = .data[[feature]])) +
    geom_point(size = 0.2, alpha = 0.5) +
    facet_wrap(~rep) +
    labs(x = "Bin de actividad media", y = "Rango de varianza (dentro del bin)", col = leyenda, title = titulo) +
    scale_color_manual(values = c("grey", thesis_clr), guide = guide_legend(override.aes = list(size = 3))) +
    ggnewscale::new_scale_color() +
    geom_smooth(aes(col = .data[[feature]]), se = FALSE, linewidth = 1.2, show.legend = FALSE) +
    scale_color_manual(values = c("#525252", "#136869")) +
    theme_pubr(base_size = 16)
}

# Curva ROC (una linea por replica, con AUC en el titulo) de una feature
# booleana prediciendo var_rank_sw alto.
noise_roc <- function(data, feature, titulo) {
  data %>%
    group_split(rep) %>%
    purrr::map(~ roc_curve(.x, feature) %>% mutate(rep = unique(.x$rep), AUC = noise_auc(.))) %>%
    purrr::list_rbind() %>%
    ggplot(aes(FPR, TPR, col = rep)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", col = "grey50") +
    geom_text(
      data = . %>% distinct(AUC, rep) %>% mutate(AUC = paste("AUC:", round(AUC, 3))) %>%
        bind_cols(tibble(TPR = c(0.9, 0.8), FPR = 0.02)),
      aes(label = AUC), hjust = "left", show.legend = FALSE
    ) +
    ggtitle(titulo) +
    scale_color_manual(values = c("#0D2C54", "#AD343E")) +
    theme_pubr(base_size = 16) +
    labs(col = "")
}

# Como noise_roc(), pero separando ademas por nivel de actividad (bajo el
# percentil 75 de mean_sw vs. por encima) via linetype - una curva por
# combinacion replica x nivel.
noise_roc_by_expr <- function(data, feature, titulo) {
  data %>%
    mutate(expr_level = ifelse(mean_sw < 0.75 * max(mean_sw), "Baja actividad", "Alta actividad")) %>%
    group_split(rep, expr_level) %>%
    purrr::map(~ roc_curve(.x, feature) %>% mutate(rep = unique(.x$rep), expr_level = unique(.x$expr_level), AUC = noise_auc(.))) %>%
    purrr::list_rbind() %>%
    ggplot(aes(FPR, TPR, col = rep, linetype = expr_level)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", col = "grey50") +
    geom_text(
      data = . %>% distinct(AUC, rep, expr_level) %>% arrange(rep, expr_level) %>%
        mutate(AUC = paste("AUC:", round(AUC, 3))) %>%
        bind_cols(tibble(TPR = c(0.9, 0.8, 0.7, 0.6), FPR = 0.02)),
      aes(label = AUC), hjust = "left", show.legend = FALSE
    ) +
    ggtitle(titulo) +
    scale_color_manual(values = c("#0D2C54", "#AD343E")) +
    theme_pubr(base_size = 16) +
    labs(col = "", linetype = "")
}

# Tabla de ~32 features booleanas del promotor (nombres en espanol),
# compartida entre R7 (summary_features.R, efecto sobre actividad) y R3.2
# (ruido_summary.R, efecto sobre ruido). `keep` son las columnas no-feature
# a conservar ademas de rep (ej. "mean" para R7, "var_rank_sw" para R3.2).
# `cgi_col` es la columna de prom_df que define "Islas CpG": "CGI" (anotacion
# genomica UCSC, default) o "CGI_frag" (composicion del fragmento de 252 pb,
# usada en summary_features_secuencia desde 2026-10-08).
build_tidy_features <- function(data, keep, cgi_col = "CGI") {
  data %>%
    fastDummies::dummy_cols("TE_superclass", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = FALSE) %>%
    fastDummies::dummy_cols("sample_specificity_class", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = TRUE) %>%
    group_by(rep) %>%
    mutate(
      `Alto contenido G+C` = (cut_number(g_c, n = 3) %>% as.numeric()) == 3,
      across(c(DNA, SINE, LINE, LTR), ~ replace_na(.x, 0)),
      `Sin actividad en ratón` = turnover %in% c("expression-turnover", "mouse-diminished"),
      `Insertado en humanos` = turnover == "human-inserted",
      `Elementos transponibles` = !is.na(TE_superclass),
      `Repeticiones de baja complejidad` = LCR_overlap > 0,
      `Alta conservación (16 a -50pb)` = (cut_number(phylop100_close, n = 3) %>% as.numeric()) == 3,
      `Alta conservación (-50 a -150)` = (cut_number(phylop100_intermediate, n = 3) %>% as.numeric()) == 3,
      `Alta conservación (-150 a -235)` = (cut_number(phylop100_far, n = 3) %>% as.numeric()) == 3,
      `Alta especificidad tisular` = (cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) == 3,
      `Baja especificidad tisular` = (cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) == 1,
      shape_class_n = cut_number(interquantile_width, 3) %>% as_factor() %>% as.numeric(),
      `Promotores angostos` = shape_class_n == 1,
      `Promotores anchos` = shape_class_n == 3,
      `Alta actividad en HEK293` = hek_tpm > median(data$hek_tpm[data$hek_tpm > 0], na.rm = TRUE),
      `Sin módulo cis-regulatorio` = N_TF_CRM == 0,
      `Alta accesibilidad de cromatina (DNase-seq)` = (cut_number(mean_dnase, n = 3) %>% as.numeric()) == 3,
      across(starts_with("enh"), ~ .x > 0),
      `Promotor unidireccional` = case_when(
        orientacion == "unidireccional" ~ TRUE,
        orientacion == "bidireccional" ~ FALSE,
        TRUE ~ NA
      )
    ) %>%
    rename(
      `TATA-box` = TATA_EPD,
      CCAAT = CCAAT_EPD,
      `GC-box` = GCbox_EPD,
      `Islas CpG` = all_of(cgi_col),
      Retrotransposón = DNA,
      TCT = TCT_TSS,
      `CG en TSS` = CG_TSS,
      `TA en TSS` = TA_TSS,
      `TG en TSS` = TG_TSS,
      `CA en TSS` = CA_TSS,
      `TSS no canónico` = other_TSS,
      `TSS fuerte` = INR_strong_TSS,
      `No detectado (FANTOM5)` = non_detected,
      `Enhancers a 10kb` = enh10kb,
      `Enhancers a 50kb` = enh50kb,
      `Enhancers a 100kb` = enh100kb
    ) %>%
    select(
      all_of(keep), rep,
      LINE, SINE, LTR,
      `Alto contenido G+C`, `Sin actividad en ratón`, `Insertado en humanos`,
      `Elementos transponibles`, `Repeticiones de baja complejidad`,
      `Alta conservación (16 a -50pb)`, `Alta conservación (-50 a -150)`, `Alta conservación (-150 a -235)`,
      `TATA-box`, CCAAT, `GC-box`, `Islas CpG`, Retrotransposón, TCT,
      `CG en TSS`, `TA en TSS`, `TG en TSS`, `CA en TSS`,
      `TSS no canónico`, `TSS fuerte`,
      `No detectado (FANTOM5)`,
      `Alta especificidad tisular`, `Baja especificidad tisular`,
      `Promotores angostos`, `Promotores anchos`,
      `Alta actividad en HEK293`, `Alta accesibilidad de cromatina (DNase-seq)`,
      `Enhancers a 50kb`, `Sin módulo cis-regulatorio`,
      `Promotor unidireccional`
    ) %>%
    mutate(
      across(c(where(is.numeric), -all_of(keep)), as.logical),
      rep = as.factor(rep),
      across(where(is.logical), ~ .x %>% factor(levels = c("TRUE", "FALSE")))
    ) %>%
    ungroup()
}

# Barras de AUC-0.5 por feature/TF, coloreadas por sentido del efecto
# sobre el ruido (ruido alto/bajo), filtradas a las que son significativas
# con el mismo sentido en ambas replicas. Usado por R3.2 (features
# curadas, show_labels=TRUE), R3.3 (TFs de ReMap, show_labels=FALSE,
# show_errorbars=FALSE - demasiados para etiquetar o mostrar IC, igual
# que remap_act.jpg en R2.6) y R3.4.
# Dos criterios de significancia soportados segun las columnas de
# auc_df: si trae `pval_corr` (p-valor de Wald via varianza de DeLong,
# corregido por BH - ver R3.2/ruido_summary.R, resync 2026-09-30) se usa
# pval_corr<0.05; si no (R3.3/R3.4, que todavia no migraron) se usa el
# criterio viejo de IC (delong) que no cruza 0.5 - equivalentes en
# esencia (mismo SE de DeLong), pero el de p-valor permite corregir por
# tests multiples.
plot_noise_auc_summary <- function(auc_df, show_labels = TRUE, show_errorbars = TRUE, base_size = 20) {
  auc_df_signif <- if ("pval_corr" %in% names(auc_df)) {
    auc_df %>% filter(pval_corr < 0.05)
  } else {
    auc_df %>% filter((ci2.5 > 0.5 & ci97.5 > 0.5) | (ci2.5 < 0.5 & ci97.5 < 0.5))
  }
  p <- auc_df_signif %>%
    mutate(noise = ifelse(AUC > 0.5, "Ruido alto", "Ruido bajo")) %>%
    group_by(feature) %>%
    mutate(n_dir = length(unique(noise))) %>%
    filter(n_dir == 1, n() == 2) %>%
    ungroup() %>%
    arrange(desc(rep), AUC) %>%
    mutate(feature = fct_inorder(feature)) %>%
    ggplot(aes(x = AUC - 0.5, y = feature, group = fct_inorder(rep))) +
    geom_col(orientation = "y", position = "dodge", aes(fill = noise, alpha = rep)) +
    scale_x_continuous(labels = function(x) x + 0.5) +
    scale_alpha_manual(values = c("Rep 1" = 1, "Rep 2" = 0.7)) +
    scale_fill_manual(values = c("Ruido alto" = "#D6741F", "Ruido bajo" = "#7FB800")) +
    theme_pubr(base_size = base_size) +
    theme(legend.position = "top") +
    labs(fill = "Efecto", x = "AUC (efecto sobre el ruido)", y = "Features", alpha = "")
  if (show_errorbars) p <- p + geom_errorbarh(aes(xmax = ci2.5 - 0.5, xmin = ci97.5 - 0.5), position = position_dodge(1), height = 0.05, col = "#777777", linewidth = 1.5)
  if (!show_labels) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}

# Agrupa cada feature de build_tidy_features() en "seq" (secuencia del
# promotor) o "endo" (contexto endogeno/genomico) - migrado de
# transcriptional_library/Analysis/Tables/tidy_names.tsv. Usado por R7.
feature_groups <- tribble(
  ~feature, ~group,
  "LINE", "seq", "SINE", "seq", "LTR", "seq",
  "Alto contenido G+C", "seq",
  "Elementos transponibles", "seq",
  "Repeticiones de baja complejidad", "seq",
  "TATA-box", "seq", "CCAAT", "seq", "GC-box", "seq", "Islas CpG", "seq",
  "Retrotransposón", "seq", "TCT", "seq",
  "CG en TSS", "seq", "TA en TSS", "seq", "TG en TSS", "seq", "CA en TSS", "seq",
  "TSS no canónico", "seq", "TSS fuerte", "seq",
  "Sin actividad en ratón", "endo", "Insertado en humanos", "endo",
  "Alta conservación (16 a -50pb)", "endo", "Alta conservación (-50 a -150)", "endo", "Alta conservación (-150 a -235)", "endo",
  "No detectado (FANTOM5)", "endo",
  "Alta especificidad tisular", "endo", "Baja especificidad tisular", "endo",
  "Promotores angostos", "endo", "Promotores anchos", "endo",
  "Alta actividad en HEK293", "endo",
  "Alta accesibilidad de cromatina (DNase-seq)", "endo",
  "Enhancers a 50kb", "endo",
  "Sin módulo cis-regulatorio", "endo",
  "Promotor unidireccional", "endo"
)

# Carga y parsea los .cells.csv (export FlowJo) de la validacion por
# citometria de 8 promotores individuales + 2 controles (Control="US",
# Strong) - compartido por R1.5 (densidad_promotores_individuales) y R1.6
# (distribuciones_expresion). El orden del vector de nombres tiene que
# coincidir con el orden en que unique(df$prom_name) devuelve los
# short_name - fragil pero heredado tal cual de
# transcriptional_library/Tesis/tesis.R.
load_citometry_stable_validation <- function(path_citometry_stable_validation) {
  files <- list.files(path_citometry_stable_validation, pattern = ".cells", full.names = TRUE, recursive = TRUE)
  prom_name <- files %>%
    str_remove("^.+Tables/(lvs-)?") %>%
    str_remove("_Data .+$")
  df <- map2(files, prom_name, ~ read_csv(.x, show_col_types = FALSE) %>% mutate(sample_name = .y)) %>% list_rbind()
  names(df) <- str_replace_all(names(df), "-", "_")
  df$prom_name <- str_remove(df$sample_name, " -.+$")

  name_df <- tibble(
    short_name = unique(df$prom_name),
    name = c("BTG1_1", "ETS1_1", "KIAA0753_1", "LSM1_1", "METAP2_1", "PPP1R14B_3", "TMEM87A_1", "ZKSCAN2_1", "Control", "Strong")
  )
  inner_join(name_df, df, by = c("short_name" = "prom_name"))
}

# Efecto de cada feature booleana (columnas factor TRUE/FALSE) de
# tidy_data sobre `mean` (actividad media): p-valor combinado (test de
# Wilcoxon pareado por replica, coin::wilcox_test(mean~.x|rep), con
# correccion BH) y tamano de efecto + IC95% por separado en cada
# replica (coin::wilcox_test(mean~.x, conf.int=TRUE)). Usado por R7
# (summary_features_secuencia) y R5.3 (promalt_summary, con "X main"/
# "X" como las features booleanas en vez de features de secuencia).
wilcox_effect_summary <- function(tidy_data) {
  wilcox <- map(tidy_data %>% select(-rep, -contains("mean")), ~ coin::wilcox_test(formula = mean ~ .x | rep, data = tidy_data) %>% pvalue())
  wilcox <- tibble(feature = names(wilcox), pval = list_c(wilcox), pval_corr = p.adjust(pval, "BH", length(wilcox)))

  wilcox_by_rep <- function(rep_id) {
    map(tidy_data %>% filter(rep == rep_id) %>% select(-mean, -rep), ~ coin::wilcox_test(formula = mean ~ .x, data = tidy_data %>% filter(rep == rep_id), conf.int = TRUE))
  }
  wilcox_rep1 <- wilcox_by_rep("Rep 1")
  wilcox_rep2 <- wilcox_by_rep("Rep 2")

  map2(
    list(wilcox_rep1, wilcox_rep2), c("Rep 1", "Rep 2"),
    ~ tibble(
      feature = names(.x),
      estimate = map(.x, ~ confint(.x)$estimate) %>% list_c(),
      P2.5 = map(.x, ~ confint(.x)$conf.int[1]) %>% list_c(),
      P97.5 = map(.x, ~ confint(.x)$conf.int[2]) %>% list_c(),
      rep = .y
    ) %>%
      pivot_longer(c(starts_with("estimate"), starts_with("P2.5"), starts_with("P97.5")), names_to = "val", values_to = "estimate") %>%
      arrange(desc(estimate)) %>%
      mutate(feature = fct_inorder(feature))
  ) %>%
    list_rbind() %>%
    left_join(wilcox, by = "feature")
}

# Tabla de features booleanas del promotor (mismos nombres en espanol
# que build_tidy_features()/R7, pero por seq_id en vez de por rep -
# usada para co-ocurrencia, no para asociacion con actividad/ruido).
# Requiere prom_df.tsv completo (no filtrado a type=="promoter"). Usado
# por M14 (coocurrencia_motivos) y R5.7 (promalt_coocurrencia).
build_cooccurrence_features <- function(prom_df) {
  prom_df %>%
    fastDummies::dummy_cols("TE_superclass", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = FALSE) %>%
    fastDummies::dummy_cols("sample_specificity_class", ignore_na = TRUE, omit_colname_prefix = TRUE, remove_selected_columns = TRUE) %>%
    mutate(
      `Alto contenido G+C` = (cut_number(g_c, n = 3) %>% as.numeric()) == 3,
      across(c(DNA, SINE, LINE, LTR), ~ replace_na(.x, 0)),
      `Sin actividad en ratón` = turnover %in% c("expression-turnover", "mouse-diminished"),
      `Insertado en humanos` = turnover == "human-inserted",
      `Elementos transponibles` = !is.na(TE_superclass),
      `Repeticiones de baja complejidad` = LCR_overlap > 0,
      `Alta conservación (16 a -50pb)` = (cut_number(phylop100_close, n = 3) %>% as.numeric()) == 3,
      `Alta conservación (-50 a -150)` = (cut_number(phylop100_intermediate, n = 3) %>% as.numeric()) == 3,
      `Alta conservación (-150 a -235)` = (cut_number(phylop100_far, n = 3) %>% as.numeric()) == 3,
      `Alta especificidad tisular` = (cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) == 3,
      `Baja especificidad tisular` = (cut_number(sample_specificity_gini, n = 3) %>% as.numeric()) == 1,
      shape_class_n = cut_number(interquantile_width, 3) %>% as_factor() %>% as.numeric(),
      `Promotores angostos` = shape_class_n == 1,
      `Promotores anchos` = shape_class_n == 3,
      `Alta actividad en HEK293` = hek_tpm > median(prom_df$hek_tpm[prom_df$hek_tpm > 0], na.rm = TRUE),
      `Sin módulo cis-regulatorio` = N_TF_CRM == 0,
      `Alta accesibilidad de cromatina (DNase-seq)` = (cut_number(mean_dnase, n = 3) %>% as.numeric()) == 3,
      across(starts_with("enh"), ~ .x > 0),
      `Promotor unidireccional` = case_when(
        orientacion == "unidireccional" ~ TRUE,
        orientacion == "bidireccional" ~ FALSE,
        TRUE ~ NA
      )
    ) %>%
    rename(
      `TATA-box` = TATA_EPD,
      CCAAT = CCAAT_EPD,
      `GC-box` = GCbox_EPD,
      `Islas CpG` = CGI,
      Retrotransposón = DNA,
      TCT = TCT_TSS,
      `CG en TSS` = CG_TSS,
      `TA en TSS` = TA_TSS,
      `TG en TSS` = TG_TSS,
      `CA en TSS` = CA_TSS,
      `TSS no canónico` = other_TSS,
      `TSS fuerte` = INR_strong_TSS,
      `No detectado (FANTOM5)` = non_detected,
      `Enhancers a 10kb` = enh10kb,
      `Enhancers a 50kb` = enh50kb,
      `Enhancers a 100kb` = enh100kb
    ) %>%
    select(
      seq_id,
      LINE, SINE, LTR,
      `Alto contenido G+C`, `Sin actividad en ratón`, `Insertado en humanos`,
      `Elementos transponibles`, `Repeticiones de baja complejidad`,
      `Alta conservación (16 a -50pb)`, `Alta conservación (-50 a -150)`, `Alta conservación (-150 a -235)`,
      `TATA-box`, CCAAT, `GC-box`, `Islas CpG`, Retrotransposón, TCT,
      `CG en TSS`, `TA en TSS`, `TG en TSS`, `CA en TSS`,
      `TSS no canónico`, `TSS fuerte`,
      `No detectado (FANTOM5)`,
      `Alta especificidad tisular`, `Baja especificidad tisular`,
      `Promotores angostos`, `Promotores anchos`,
      `Alta actividad en HEK293`, `Alta accesibilidad de cromatina (DNase-seq)`,
      `Enhancers a 50kb`, `Sin módulo cis-regulatorio`,
      `Promotor unidireccional`
    )
}

# Matriz de coeficiente phi (correlacion de Pearson entre columnas
# booleanas 0/1) y de similitud de Jaccard, mas su heatmap - usado por
# M14 y R5.7 para co-ocurrencia entre features booleanas del promotor
# (incluida la categoria de promotor alternativo en R5.7).
cooccurrence_phi_matrix <- function(feat_mat) cor(feat_mat, use = "pairwise.complete.obs")

# Variante de cooccurrence_phi_matrix() SOLO para calcular el orden del
# clustering jerarquico (ver plot_cooccurrence_phi() / mismo patron que
# cooccurrence_lift_matrix(..., mask_zero=FALSE)): a veces cor() da NA en
# una celda puntual porque, dentro del subconjunto de filas donde ese par
# especifico tiene ambas columnas no-NA, una de las dos queda con
# varianza cero (ej. sesion 2026-09-30: "No detectado (FANTOM5)", "Alta
# especificidad tisular", "Baja especificidad tisular"). Como
# cooccurrence_cluster_levels() descarta del clustering CUALQUIER feature
# con al menos una celda NA, una sola celda asi tira afuera a toda la
# feature. Se imputan esas celdas puntuales a 0 (sin evidencia de
# asociacion) solo para definir el orden - el heatmap real se sigue
# pintando con cooccurrence_phi_matrix() sin imputar (NA queda gris).
cooccurrence_phi_matrix_for_clustering <- function(feat_mat) {
  m <- cooccurrence_phi_matrix(feat_mat)
  m[is.na(m)] <- 0
  m
}

cooccurrence_jaccard_matrix <- function(feat_mat) {
  m <- matrix(NA, nrow = ncol(feat_mat), ncol = ncol(feat_mat), dimnames = list(colnames(feat_mat), colnames(feat_mat)))
  for (i in seq_len(ncol(feat_mat))) {
    for (j in seq_len(ncol(feat_mat))) {
      a <- feat_mat[, i]
      b <- feat_mat[, j]
      intersection <- sum(a == 1 & b == 1, na.rm = TRUE)
      union <- sum(a == 1 | b == 1, na.rm = TRUE)
      m[i, j] <- if (union == 0) NA else intersection / union
    }
  }
  m
}

# Coeficiente de overlap (Szymkiewicz-Simpson): |A∩B| / min(|A|,|B|), en
# vez de la union (Jaccard). A diferencia de Jaccard, no penaliza cuando
# un conjunto chico esta casi totalmente contenido en uno mucho mas
# grande (ej. TATA-box practicamente un subconjunto de "Promotores
# angostos", pero Jaccard da bajo por la asimetria de tamanios - ver
# sesion 2026-09-30, M15/coocurrencia_motivos).
cooccurrence_overlap_matrix <- function(feat_mat) {
  m <- matrix(NA, nrow = ncol(feat_mat), ncol = ncol(feat_mat), dimnames = list(colnames(feat_mat), colnames(feat_mat)))
  for (i in seq_len(ncol(feat_mat))) {
    for (j in seq_len(ncol(feat_mat))) {
      a <- feat_mat[, i]
      b <- feat_mat[, j]
      intersection <- sum(a == 1 & b == 1, na.rm = TRUE)
      min_size <- min(sum(a == 1, na.rm = TRUE), sum(b == 1, na.rm = TRUE))
      m[i, j] <- if (min_size == 0) NA else intersection / min_size
    }
  }
  m
}

# Lift/enriquecimiento (observado/esperado bajo independencia):
# P(A∩B) / (P(A)*P(B)). Simetrico (a diferencia del overlap), pero NO
# tiene el "techo" de phi cuando las prevalencias estan lejos de 50/50 -
# aporta informacion que phi no muestra para relaciones fuertes entre
# features raras (ver sesion 2026-09-30, M15/coocurrencia_motivos: ej.
# LTR+Insertado en humanos da lift=15x pero phi=0.16, mientras que
# Baja especificidad tisular+Alta accesibilidad da phi=0.51 pero
# lift=1.9x - el orden por cada metrica es bien distinto, Spearman~0.74
# entre las dos). Correccion de Haldane-Anscombe (pseudocount EPS=0.5 en
# las 4 celdas de la tabla 2x2) para pares con interseccion real >0 pero
# chica (evita inestabilidad numerica); para pares con interseccion real
# EXACTAMENTE 0 (mutuamente excluyentes, ej. subtipos de TE) NO se usa
# pseudocount - el lift "gris" (NA) en vez de un valor cerca de 0 pero
# artificialmente inflado/opacando la escala de color del resto (ver
# sesion 2026-09-30, feedback del autor).
# mask_zero=TRUE (default, para pintar el heatmap): pares con interseccion
# real 0 quedan en NA/gris, sin pseudocount. mask_zero=FALSE (para calcular
# el orden del clustering unicamente - ver plot_cooccurrence_lift()): usa
# el pseudocount tambien en esos pares, matriz completa sin NA, para que
# esas features no queden afuera del clustering jerarquico.
cooccurrence_lift_matrix <- function(feat_mat, eps = 0.5, mask_zero = TRUE) {
  m <- matrix(NA, nrow = ncol(feat_mat), ncol = ncol(feat_mat), dimnames = list(colnames(feat_mat), colnames(feat_mat)))
  for (i in seq_len(ncol(feat_mat))) {
    for (j in seq_len(ncol(feat_mat))) {
      a <- feat_mat[, i]
      b <- feat_mat[, j]
      ok <- !is.na(a) & !is.na(b)
      n11_raw <- sum(a[ok] == 1 & b[ok] == 1)
      if (n11_raw == 0 && mask_zero) next # deja NA (gris) - sin pseudocount para no distorsionar la escala
      n10 <- sum(a[ok] == 1 & b[ok] == 0) + eps
      n01 <- sum(a[ok] == 0 & b[ok] == 1) + eps
      n00 <- sum(a[ok] == 0 & b[ok] == 0) + eps
      n11 <- n11_raw + eps
      n_tot <- n11 + n10 + n01 + n00
      pa <- (n11 + n10) / n_tot
      pb <- (n11 + n01) / n_tot
      pab <- n11 / n_tot
      m[i, j] <- if (pa == 0 || pb == 0) NA else pab / (pa * pb)
    }
  }
  m
}

cooccurrence_cluster_levels <- function(mat) {
  has_na <- colSums(is.na(mat)) > 0
  complete <- mat[!has_na, !has_na]
  clust_order <- hclust(as.dist(1 - complete))$order
  c(colnames(complete)[clust_order], colnames(mat)[has_na])
}

# levels=NULL (default): orden de clustering calculado sobre la misma
# matriz que se va a pintar (comportamiento original). Si se pasa un
# vector de niveles ya calculado (ej. desde una matriz sin enmascarar,
# ver plot_cooccurrence_lift()), se usa ese orden en vez de recalcularlo -
# permite desacoplar "con que datos clusterizo" de "que pinto".
cooccurrence_tidy_matrix <- function(mat, value_name, levels = NULL) {
  if (is.null(levels)) levels <- cooccurrence_cluster_levels(mat)
  mat %>%
    as_tibble(rownames = "feature1") %>%
    pivot_longer(-feature1, names_to = "feature2", values_to = value_name) %>%
    mutate(feature1 = factor(feature1, levels = levels), feature2 = factor(feature2, levels = levels))
}

plot_cooccurrence_phi <- function(feat_mat, titulo = "Co-ocurrencia de features del promotor (coeficiente phi)") {
  levels <- cooccurrence_cluster_levels(cooccurrence_phi_matrix_for_clustering(feat_mat))
  cooccurrence_tidy_matrix(cooccurrence_phi_matrix(feat_mat), "phi", levels = levels) %>%
    ggplot(aes(x = feature2, y = feature1, fill = phi)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0, limits = c(-1, 1), name = "Coeficiente\nphi") +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
      axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
    ) +
    labs(title = titulo)
}

plot_cooccurrence_jaccard <- function(feat_mat, titulo = "Co-ocurrencia de features del promotor (similitud Jaccard)") {
  cooccurrence_tidy_matrix(cooccurrence_jaccard_matrix(feat_mat), "jaccard") %>%
    ggplot(aes(x = feature2, y = feature1, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#d6604d", limits = c(0, 1), na.value = "grey80", name = "Similitud\nJaccard") +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
      axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
    ) +
    labs(title = titulo)
}

plot_cooccurrence_overlap <- function(feat_mat, titulo = "Co-ocurrencia de features del promotor (coeficiente de overlap)") {
  cooccurrence_tidy_matrix(cooccurrence_overlap_matrix(feat_mat), "overlap") %>%
    ggplot(aes(x = feature2, y = feature1, fill = overlap)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#d6604d", limits = c(0, 1), na.value = "grey80", name = "Coeficiente\nde overlap") +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
      axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
    ) +
    labs(title = titulo)
}

plot_cooccurrence_lift <- function(feat_mat, titulo = "Co-ocurrencia de features del promotor (log2 lift)") {
  # el orden de clustering se calcula sobre la matriz COMPLETA (con
  # pseudocount tambien en pares de interseccion 0), para que esos pares
  # (ej. subtipos de TE mutuamente excluyentes) no queden afuera del
  # clustering jerarquico solo por estar enmascarados en el heatmap
  levels <- cooccurrence_cluster_levels(log2(cooccurrence_lift_matrix(feat_mat, mask_zero = FALSE)))
  lift_mat <- cooccurrence_lift_matrix(feat_mat)
  log2_lift_mat <- log2(lift_mat)
  lim <- max(abs(log2_lift_mat), na.rm = TRUE)
  cooccurrence_tidy_matrix(log2_lift_mat, "log2_lift", levels = levels) %>%
    ggplot(aes(x = feature2, y = feature1, fill = log2_lift)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d", midpoint = 0, limits = c(-lim, lim), na.value = "grey80", name = "log2\n(lift)") +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
      axis.title = element_blank(), panel.grid = element_blank(), legend.position = "right"
    ) +
    labs(title = titulo)
}

# Frecuencia de cada feature booleana de feat_mat sobre SU PROPIO universo
# (promotores no-NA para esa columna especifica, que varia entre
# features - ej. "Promotor unidireccional" solo esta definido para los
# promotores con senal suficiente, MIN_SENSE, ver R/28_orientacion_
# promotor). Reusa count_bar_labeled(): barra clara = universo (n no-NA),
# barra oscura = n con feature=TRUE, ordenadas por frecuencia.
plot_feature_frequency <- function(feat_mat, titulo = "Frecuencia de cada feature (barra clara = universo no-NA)") {
  freq_df <- tibble(
    feature = colnames(feat_mat),
    n_universo = colSums(!is.na(feat_mat)),
    n_true = colSums(feat_mat == 1, na.rm = TRUE)
  ) %>%
    mutate(feature = fct_reorder(feature, n_true / n_universo))
  count_bar_labeled(freq_df, x = feature, n = n_true, total = n_universo, base_size = 12, coord_flip = TRUE) +
    labs(title = titulo, x = NULL, y = "Promotores (n)")
}
