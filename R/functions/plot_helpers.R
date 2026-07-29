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
build_tidy_features <- function(data, keep) {
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
      across(starts_with("enh"), ~ .x > 0)
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
      `Enhancers a 50kb`, `Sin módulo cis-regulatorio`
    ) %>%
    mutate(
      across(c(where(is.numeric), -all_of(keep)), as.logical),
      rep = as.factor(rep),
      across(where(is.logical), ~ .x %>% factor(levels = c("TRUE", "FALSE")))
    ) %>%
    ungroup()
}

# Barras de AUC-0.5 por feature/TF, coloreadas por sentido del efecto
# sobre el ruido (ruido alto/bajo), filtradas a las que tienen IC que no
# cruza 0.5 y el mismo sentido en ambas replicas. Usado por R3.2
# (features curadas, show_labels=TRUE), R3.3 (TFs de ReMap,
# show_labels=FALSE, show_errorbars=FALSE - demasiados para etiquetar o
# mostrar IC, igual que remap_act.jpg en R2.6) y R3.4.
plot_noise_auc_summary <- function(auc_df, show_labels = TRUE, show_errorbars = TRUE, base_size = 20) {
  p <- auc_df %>%
    mutate(noise = ifelse(AUC > 0.5, "Ruido alto", "Ruido bajo")) %>%
    filter((ci2.5 > 0.5 & ci97.5 > 0.5) | (ci2.5 < 0.5 & ci97.5 < 0.5)) %>%
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
  "Sin módulo cis-regulatorio", "endo"
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
