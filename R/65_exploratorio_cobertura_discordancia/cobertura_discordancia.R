# Exploratorio de discordancia_remap_ora (R/60_discordancia_remap_ora):
# diagnostico de cobertura de ReMap sobre el efecto de cada factor en la
# discordancia reportero - endogeno, igual que el que se hizo para el ruido en
# R/51_exploratorio_remap_cobertura.
#
# Modelo (el mismo que discordancia_remap_ora): dif_signed ~ avg_rank + unido,
# en promotores activos en HEK293 (hek_tpm > 0), union con picos de cualquier
# linea celular, 949 factores (> 100 promotores unidos en ambas replicas, sin
# marcas de histona).
#   - Por replica: MCO dentro de cada replica (una fila por promotor), errores
#     HC1. Se calcula en forma matricial (Frisch-Waugh-Lovell) y se verifica
#     contra lm + sandwich.
#   - Conjunto: el modelo de la figura (dif_signed ~ avg_rank + rep + unido,
#     SE agrupados por promotor), leido de tables/discordancia_remap_lmm.tsv
#     (salida de R/60) y verificado contra un recalculo para los marcadores.
#   - Crudo (sin avg_rank): media de dif_signed de los unidos - la de los no
#     unidos, para ver si el control por avg_rank atenua el problema.
#
# p1 Spearman coeficiente - n de promotores unidos, por replica.
# p2 Descomposicion: media de dif_signed (y del residuo de dif_signed ~
#    avg_rank) de los unidos y de los no unidos vs. cobertura.
# p3 Crudo vs. ajustado por avg_rank: Spearman con la cobertura de cada uno.
# p4 Correccion por cobertura (como R/51): loess robusto (family =
#    "symmetric", span 0.5) del coeficiente sobre log10(n unidos). Se usan:
#      - z_cobertura = (coeficiente - esperado por cobertura) / SE del factor,
#        p normal bilateral y BH (significancia ajustada, principal);
#      - residuo / escala local (loess de |residuo|) para el ranking;
#      - percentil del residuo en la ventana de cobertura +-25%, p por
#        permutacion exacta bilateral y BH (no parametrico).
#    Posicion de SETDB1, NCOR1 y L3MBTL2 en cobertura y en el ranking antes y
#    despues del ajuste.
# p5 Conjunto significativo de cada lado (FDR < 0.05) antes y despues, y ORA
#    (enrichGO, ont = "all", universo = factores testeados, como R/60) con los
#    conjuntos ajustados.
# p6 Cobertura de los 263 factores que pasan el filtro solo-HEK293 del
#    criterio de robustez frente a los 949.
#
# Sin remuestreo aleatorio (los p por permutacion son exactos sobre la
# ventana), asi que no hace falta set.seed.
#
# Salida: figures/exploratorio_cobertura_discordancia/ (jpg y pdf, tsv) y
# handoff.md escrito a mano. Requiere data/processed/{remap_tf_hits,
# activity_stats_highconf,prom_df,fantom_endo_activity_summary}.tsv y
# tables/discordancia_remap_{lmm,robustez_hek}.tsv. Run from the TesisDoc repo
# root.

library(tidyverse)
library(Matrix)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
source("R/functions/plot_helpers.R")

select <- dplyr::select
filter <- dplyr::filter

out_dir <- "figures/exploratorio_cobertura_discordancia"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

marcadores <- c("SETDB1", "NCOR1", "L3MBTL2")
cov_window <- 0.25
loess_span <- 0.5
tramos <- c(100, 300, 1000, 2000, 3000, 4000, 5000, Inf)

save_jpg_pdf <- function(plot, name, width, height) {
  for (ext in c("jpg", "pdf")) {
    ggsave(file.path(out_dir, paste0(name, ".", ext)), plot,
      width = width, height = height, units = "in",
      device = if (ext == "pdf") cairo_pdf else "jpeg"
    )
  }
}
tema <- theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey90", color = NA), strip.text = element_text(face = "bold"))
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = "")
rho <- function(x, y) if (all(is.na(x))) NA_real_ else cor(x, y, method = "spearman", use = "complete.obs")
p_rho <- function(x, y) cor.test(x, y, method = "spearman", exact = FALSE)$p.value

# --- Datos (igual que R/60) -------------------------------------------------

remap_todos <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE)
hits_full <- remap_todos %>% distinct(seq_id, TF)
hits_hek <- remap_todos %>% filter(str_detect(sample, "HEK293")) %>% distinct(seq_id, TF)
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fantom_summary <- read_tsv("data/processed/fantom_endo_activity_summary.tsv", show_col_types = FALSE)

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  left_join(fantom_summary, by = "seq_id") %>%
  filter(!is.na(hek_tpm), hek_tpm > 0) %>%
  add_rank_discordance() %>%
  select(seq_id, rep, avg_rank, dif_signed, rr, re)
repname <- sort(unique(data$rep))
etiqueta_rep <- set_names(str_replace(repname, "Rep", "Réplica"), repname)
message("n promotores: ", paste(repname, table(data$rep)[repname], collapse = " / "))

factores_con <- function(hits) {
  hits %>%
    semi_join(data, by = "seq_id") %>%
    inner_join(distinct(data, seq_id, rep), by = "seq_id", relationship = "many-to-many") %>%
    count(rep, TF) %>%
    filter(n > 100) %>%
    count(TF, name = "n_rep") %>%
    filter(n_rep == 2, !str_detect(TF, "^H\\d")) %>%
    pull(TF) %>%
    sort()
}
tfs <- factores_con(hits_full)
tfs_hek <- factores_con(hits_hek)

lmm <- read_tsv("tables/discordancia_remap_lmm.tsv", show_col_types = FALSE)
robustez <- read_tsv("tables/discordancia_remap_robustez_hek.tsv", show_col_types = FALSE)
stopifnot(setequal(tfs, lmm$TF), length(tfs) == 949, setequal(tfs_hek, robustez$TF), length(tfs_hek) == 263)
sets <- hits_full %>% filter(TF %in% tfs) %>% with(split(seq_id, TF))
sets <- sets[tfs]

# --- Coeficientes por replica (FWL, HC1) -------------------------------------

coef_rep <- function(d) {
  N <- nrow(d)
  X <- cbind(1, d$avg_rank)
  XtXinv <- solve(crossprod(X))
  y <- d$dif_signed
  y_res <- as.numeric(y - X %*% (XtXinv %*% crossprod(X, y)))
  idx <- map(sets, ~ which(d$seq_id %in% .x))
  map2(names(idx), idx, function(tf, b) {
    x <- numeric(N)
    x[b] <- 1
    cc <- XtXinv %*% c(length(b), sum(d$avg_rank[b]))
    x_res <- x - as.numeric(X %*% cc)
    xx <- sum(x_res^2)
    beta <- sum(y_res[b]) / xx
    e <- y_res - beta * x_res
    se <- sqrt(sum(x_res^2 * e^2) / xx^2 * N / (N - 3))
    tibble(
      TF = tf, n_unidos = length(b), frac_unidos = length(b) / N,
      estimate = beta, std.error = se, statistic = beta / se,
      p_value = 2 * pt(-abs(beta / se), N - 3),
      crudo = mean(y[b]) - mean(y[-b]),
      media_unidos = mean(y[b]), media_no_unidos = mean(y[-b]),
      media_res_unidos = mean(y_res[b]), media_res_no_unidos = mean(y_res[-b]),
      avg_rank_unidos = mean(d$avg_rank[b]), avg_rank_no_unidos = mean(d$avg_rank[-b])
    )
  }) %>% list_rbind()
}

por_rep <- map(repname, ~ coef_rep(filter(data, rep == .x)) %>% mutate(rep = .x)) %>%
  list_rbind() %>%
  group_by(rep) %>%
  mutate(p_adj = p.adjust(p_value, "BH")) %>%
  ungroup()

# verificacion contra lm + sandwich
for (r in repname) {
  d <- filter(data, rep == r) %>% mutate(presente = seq_id %in% sets[["SETDB1"]])
  fit <- lm(dif_signed ~ avg_rank + presente, data = d)
  se <- sqrt(sandwich::vcovHC(fit, type = "HC1")["presenteTRUE", "presenteTRUE"])
  m <- filter(por_rep, rep == r, TF == "SETDB1")
  stopifnot(abs(coef(fit)["presenteTRUE"] - m$estimate) < 1e-10, abs(se - m$std.error) < 1e-10)
}

# conjunto: modelo de la figura (R/60) + cobertura en promotores distintos
seqs <- unique(data$seq_id)
conjunto <- lmm %>%
  select(TF, estimate, std.error, statistic, p_value, p_adj) %>%
  mutate(n_unidos = map_int(sets[TF], ~ sum(seqs %in% .x)), frac_unidos = n_unidos / length(seqs), rep = "Conjunto")
chk <- data %>%
  mutate(presente = factor(seq_id %in% sets[["SETDB1"]], levels = c(FALSE, TRUE))) %>%
  ajustar_discordancia_robusto()
stopifnot(abs(chk$estimate - lmm$estimate[lmm$TF == "SETDB1"]) < 1e-10)

# --- Correccion por cobertura -------------------------------------------------

pos_ref <- function(x, ref) {
  n_le <- sum(ref <= x)
  n_ge <- sum(ref >= x)
  c(percentil = n_le / length(ref), p_perm = min(1, 2 * (min(n_le, n_ge) + 1) / (length(ref) + 1)))
}

corregir <- function(df) {
  df$log_n <- log10(df$n_unidos)
  fit <- loess(estimate ~ log_n, data = df, family = "symmetric", span = loess_span)
  df$esperado <- predict(fit, df)
  df$residuo <- df$estimate - df$esperado
  df$residuo_z <- df$residuo / predict(loess(abs(residuo) ~ log_n, data = df, span = 0.75), df)
  df$z_cobertura <- df$residuo / df$std.error
  df$p_cobertura <- 2 * pnorm(-abs(df$z_cobertura))
  win <- map(seq_len(nrow(df)), function(i) {
    w <- setdiff(which(abs(df$n_unidos - df$n_unidos[i]) <= cov_window * df$n_unidos[i]), i)
    pr <- pos_ref(df$residuo[i], df$residuo[w])
    tibble(n_ventana = length(w), percentil_ventana_residuo = pr[["percentil"]], p_perm_residuo = pr[["p_perm"]])
  }) %>% list_rbind()
  bind_cols(df, win) %>%
    mutate(
      p_adj_cobertura = p.adjust(p_cobertura, "BH"),
      p_adj_perm = p.adjust(p_perm_residuo, "BH"),
      percentil_cobertura = percent_rank(n_unidos),
      puesto_original = rank(-estimate, ties.method = "first"),
      puesto_ajustado = rank(-residuo_z, ties.method = "first"),
      lado = case_when(p_adj < 0.05 & estimate > 0 ~ "sobreestima", p_adj < 0.05 ~ "subestima", TRUE ~ "ns"),
      lado_ajustado = case_when(p_adj_cobertura < 0.05 & residuo > 0 ~ "sobreestima", p_adj_cobertura < 0.05 ~ "subestima", TRUE ~ "ns")
    )
}

todos <- bind_rows(por_rep, conjunto) %>%
  mutate(rep = factor(rep, levels = c(repname, "Conjunto"))) %>%
  group_split(rep) %>%
  map(corregir) %>%
  list_rbind()
niveles_rep <- c(etiqueta_rep, Conjunto = "Conjunto (modelo de la figura)")
todos <- todos %>% mutate(rep_lab = factor(niveles_rep[as.character(rep)], levels = niveles_rep))
write_tsv(todos %>% select(-rep_lab), file.path(out_dir, "estadisticos_por_factor.tsv"))

# --- p1 + p3: Spearman con la cobertura ---------------------------------------

resumen_rho <- todos %>%
  group_by(rep) %>%
  summarise(
    n_factores = n(),
    rho_ajustado = rho(estimate, n_unidos), p_rho_ajustado = p_rho(estimate, n_unidos),
    rho_crudo = rho(crudo, n_unidos),
    rho_t = rho(statistic, n_unidos),
    rho_media_unidos = rho(media_unidos, n_unidos), rho_media_no_unidos = rho(media_no_unidos, n_unidos),
    rho_res_unidos = rho(media_res_unidos, n_unidos), rho_res_no_unidos = rho(media_res_no_unidos, n_unidos),
    rho_avg_rank_unidos = rho(avg_rank_unidos, n_unidos), rho_avg_rank_no_unidos = rho(avg_rank_no_unidos, n_unidos),
    rho_residuo_loess = rho(residuo, n_unidos), rho_residuo_z = rho(residuo_z, n_unidos),
    rho_percentil_ventana = rho(percentil_ventana_residuo, n_unidos),
    .groups = "drop"
  )
write_tsv(resumen_rho, file.path(out_dir, "p1_p3_spearman_cobertura.tsv"))
print(resumen_rho %>% mutate(across(where(is.numeric) & !starts_with("p_") & !n_factores, ~ round(.x, 3))), width = Inf)

por_tramo <- todos %>%
  mutate(tramo = cut(n_unidos, tramos, right = FALSE, dig.lab = 6)) %>%
  group_by(rep, tramo) %>%
  summarise(
    n = n(), mediana_estimate = median(estimate), mediana_crudo = median(crudo),
    media_unidos = mean(media_unidos), media_no_unidos = mean(media_no_unidos),
    media_res_unidos = mean(media_res_unidos), media_res_no_unidos = mean(media_res_no_unidos),
    avg_rank_unidos = mean(avg_rank_unidos), avg_rank_no_unidos = mean(avg_rank_no_unidos),
    frac_sobreestima = mean(lado == "sobreestima"), frac_subestima = mean(lado == "subestima"),
    .groups = "drop"
  )
write_tsv(por_tramo, file.path(out_dir, "p2_por_tramo_cobertura.tsv"))
print(por_tramo %>% mutate(across(where(is.double), ~ round(.x, 3))), n = Inf, width = Inf)

lab_marc <- todos %>% filter(TF %in% marcadores)

p1 <- ggplot(todos, aes(n_unidos, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(color = lado), size = 1, alpha = 0.6) +
  geom_line(aes(y = esperado), color = "black", linewidth = 0.8) +
  geom_point(data = lab_marc, color = "black", size = 2.2, shape = 21, fill = "gold") +
  geom_text_repel(data = lab_marc, aes(label = TF), size = 3.5, min.segment.length = 0) +
  geom_text(
    data = resumen_rho %>% mutate(rep_lab = factor(niveles_rep[as.character(rep)], levels = niveles_rep)),
    aes(x = 120, y = Inf, label = paste0("ρ = ", fmt(rho_ajustado))), hjust = 0, vjust = 1.5, size = 4.5, inherit.aes = FALSE
  ) +
  scale_x_log10() +
  scale_color_manual(values = c(sobreestima = "#c0392b", subestima = "#2980b9", ns = "grey70"), labels = c(sobreestima = "Sobreestima (FDR < 0.05)", subestima = "Subestima (FDR < 0.05)", ns = "ns")) +
  facet_wrap(~rep_lab, nrow = 1) +
  labs(x = "Promotores unidos (log10)", y = "Efecto (β) sobre la discordancia", color = NULL) +
  tema +
  theme(legend.position = "bottom")
save_jpg_pdf(p1, "p1_coef_vs_cobertura", 14, 5.5)

# --- p2: descomposicion -------------------------------------------------------

decomp <- todos %>%
  filter(rep != "Conjunto") %>%
  select(TF, rep_lab, n_unidos, media_unidos, media_no_unidos, media_res_unidos, media_res_no_unidos) %>%
  pivot_longer(starts_with("media"), names_to = "var", values_to = "valor") %>%
  mutate(
    grupo = if_else(str_detect(var, "no_unidos"), "No unidos", "Unidos"),
    escala = if_else(str_detect(var, "res"), "Residuo de discordancia ~ avg_rank", "Discordancia (dif_signed)") %>%
      factor(levels = c("Discordancia (dif_signed)", "Residuo de discordancia ~ avg_rank"))
  )
p2 <- ggplot(decomp, aes(n_unidos, valor, color = grupo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1, formula = y ~ x) +
  scale_x_log10() +
  scale_color_manual(values = c("Unidos" = "#c0392b", "No unidos" = "grey40")) +
  facet_grid(escala ~ rep_lab, labeller = label_wrap_gen(width = 22)) +
  labs(x = "Promotores unidos (log10)", y = "Media del grupo", color = NULL) +
  tema +
  theme(legend.position = "bottom")
save_jpg_pdf(p2, "p2_media_unidos_no_unidos", 11, 8)

p2b <- todos %>%
  filter(rep != "Conjunto") %>%
  select(TF, rep_lab, n_unidos, avg_rank_unidos, avg_rank_no_unidos) %>%
  pivot_longer(starts_with("avg"), names_to = "grupo", values_to = "valor") %>%
  mutate(grupo = if_else(grupo == "avg_rank_unidos", "Unidos", "No unidos")) %>%
  ggplot(aes(n_unidos, valor, color = grupo)) +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1, formula = y ~ x) +
  scale_x_log10() +
  scale_color_manual(values = c("Unidos" = "#c0392b", "No unidos" = "grey40")) +
  facet_wrap(~rep_lab) +
  labs(x = "Promotores unidos (log10)", y = "avg_rank medio del grupo", color = NULL) +
  tema +
  theme(legend.position = "bottom")
save_jpg_pdf(p2b, "p2_avg_rank_unidos_no_unidos", 10, 5)

# --- p4: marcadores ------------------------------------------------------------

marc <- todos %>%
  filter(TF %in% marcadores) %>%
  select(rep, TF, n_unidos, frac_unidos, percentil_cobertura, estimate, p_adj, puesto_original,
         esperado, residuo, z_cobertura, p_adj_cobertura, residuo_z, puesto_ajustado,
         percentil_ventana_residuo, n_ventana, p_adj_perm) %>%
  arrange(TF, rep)
write_tsv(marc, file.path(out_dir, "p4_marcadores.tsv"))
print(marc %>% mutate(across(where(is.double), ~ signif(.x, 3))), width = Inf)

# correpresores / remodeladores en general: terminos de R/60 del lado sobreestima
p4 <- todos %>%
  ggplot(aes(puesto_original, puesto_ajustado)) +
  geom_abline(linetype = "dashed", color = "grey60") +
  geom_point(size = 0.6, alpha = 0.4, color = "grey40") +
  geom_point(data = lab_marc, color = "black", size = 2.2, shape = 21, fill = "gold") +
  geom_text_repel(data = lab_marc, aes(label = TF), size = 3.5, min.segment.length = 0) +
  facet_wrap(~rep_lab, nrow = 1) +
  labs(x = "Puesto por β (1 = mayor sobreestimación)", y = "Puesto ajustado por cobertura") +
  tema
save_jpg_pdf(p4, "p4_puesto_original_vs_ajustado", 14, 5)

# --- p5: conjuntos significativos y ORA -----------------------------------------

conteo <- todos %>%
  group_by(rep) %>%
  summarise(
    sobre_original = sum(lado == "sobreestima"), sobre_ajustado = sum(lado_ajustado == "sobreestima"),
    sobre_ambos = sum(lado == "sobreestima" & lado_ajustado == "sobreestima"),
    sobre_perm = sum(p_adj_perm < 0.05 & residuo > 0),
    sub_original = sum(lado == "subestima"), sub_ajustado = sum(lado_ajustado == "subestima"),
    sub_ambos = sum(lado == "subestima" & lado_ajustado == "subestima"),
    sub_perm = sum(p_adj_perm < 0.05 & residuo < 0),
    invierten = sum(lado != "ns" & lado_ajustado != "ns" & lado != lado_ajustado),
    .groups = "drop"
  )
write_tsv(conteo, file.path(out_dir, "p5_conteo_significativos.tsv"))
print(conteo, width = Inf)

cob_lado <- todos %>%
  group_by(rep, lado) %>%
  summarise(n = n(), mediana_n_unidos = median(n_unidos), .groups = "drop")
write_tsv(cob_lado, file.path(out_dir, "p5_cobertura_por_lado.tsv"))
print(cob_lado)

ids <- AnnotationDbi::select(org.Hs.eg.db, keys = tfs, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL") %>%
  distinct() %>%
  filter(!is.na(ENTREZID))
correr_ora <- function(genes_sym) {
  genes <- ids$ENTREZID[ids$SYMBOL %in% genes_sym]
  ora <- enrichGO(gene = genes, universe = ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "all", pAdjustMethod = "BH")
  if (is.null(ora)) return(tibble())
  as_tibble(as.data.frame(ora))
}
conj <- todos %>% filter(rep == "Conjunto")
ora <- expand_grid(version = c("original", "ajustado"), direccion = c("sobreestima", "subestima")) %>%
  mutate(res = map2(version, direccion, function(v, dir) {
    lado_col <- if (v == "original") conj$lado else conj$lado_ajustado
    correr_ora(conj$TF[lado_col == dir]) %>% mutate(version = v, direccion = dir)
  }))
ora_tab <- list_rbind(ora$res)
write_tsv(ora_tab, file.path(out_dir, "p5_ora_original_vs_ajustado.tsv"))

ora_comp <- ora_tab %>%
  filter(p.adjust < 0.05) %>%
  select(direccion, ONTOLOGY, ID, Description, version, p.adjust, Count) %>%
  pivot_wider(names_from = version, values_from = c(p.adjust, Count)) %>%
  arrange(direccion, p.adjust_original)
# p ajustado del ajustado tambien para los terminos que no llegan a 0.05
ora_comp <- ora_comp %>%
  select(-p.adjust_ajustado, -Count_ajustado) %>%
  left_join(ora_tab %>% filter(version == "ajustado") %>% select(direccion, ID, p.adjust_ajustado = p.adjust, Count_ajustado = Count), by = c("direccion", "ID")) %>%
  left_join(ora_tab %>% filter(version == "original") %>% select(direccion, ID, p_orig_todos = p.adjust), by = c("direccion", "ID")) %>%
  mutate(p.adjust_original = coalesce(p.adjust_original, p_orig_todos)) %>%
  select(-p_orig_todos) %>%
  mutate(estado = case_when(
    p.adjust_original < 0.05 & coalesce(p.adjust_ajustado, 1) < 0.05 ~ "se mantiene",
    p.adjust_original < 0.05 ~ "se pierde",
    TRUE ~ "aparece"
  ))
write_tsv(ora_comp, file.path(out_dir, "p5_ora_comparacion.tsv"))
print(ora_comp %>% count(direccion, estado))
print(ora_comp %>% mutate(across(starts_with("p.adj"), ~ signif(.x, 2))) %>% select(direccion, Description, p.adjust_original, p.adjust_ajustado, Count_original, Count_ajustado, estado), n = 80, width = Inf)

# --- p6: sesgo de cobertura del filtro solo-HEK293 --------------------------------

hek_n <- hits_hek %>%
  filter(TF %in% tfs) %>%
  filter(seq_id %in% seqs) %>%
  count(TF, name = "n_unidos_hek")
cob6 <- conj %>%
  select(TF, n_unidos, estimate, lado, lado_ajustado) %>%
  left_join(hek_n, by = "TF") %>%
  mutate(n_unidos_hek = replace_na(n_unidos_hek, 0L), pasa_hek = TF %in% tfs_hek, quintil = ntile(n_unidos, 5)) %>%
  left_join(robustez %>% select(TF, estimate_hek, concordante), by = "TF")
write_tsv(cob6, file.path(out_dir, "p6_cobertura_filtro_hek.tsv"))
w6 <- wilcox.test(n_unidos ~ pasa_hek, data = cob6)
p6_res <- tibble(
  n_todos = nrow(cob6), n_pasan = sum(cob6$pasa_hek),
  mediana_todos = median(cob6$n_unidos), mediana_pasan = median(cob6$n_unidos[cob6$pasa_hek]),
  mediana_no_pasan = median(cob6$n_unidos[!cob6$pasa_hek]),
  p_wilcoxon = w6$p.value,
  rho_cobertura_hek_vs_total = rho(cob6$n_unidos_hek[cob6$pasa_hek], cob6$n_unidos[cob6$pasa_hek]),
  rho_estimate_hek_vs_cobertura_hek = rho(cob6$estimate_hek[cob6$pasa_hek], cob6$n_unidos_hek[cob6$pasa_hek]),
  frac_sobre_original_que_pasan = mean(cob6$pasa_hek[cob6$lado == "sobreestima"]),
  frac_sub_original_que_pasan = mean(cob6$pasa_hek[cob6$lado == "subestima"]),
  concordantes_sobre = sum(cob6$concordante & cob6$estimate > 0, na.rm = TRUE),
  concordantes_sub = sum(cob6$concordante & cob6$estimate < 0, na.rm = TRUE),
  concordantes_sobre_que_sobreviven_ajuste = sum(cob6$concordante & cob6$estimate > 0 & cob6$lado_ajustado == "sobreestima", na.rm = TRUE),
  concordantes_sub_que_sobreviven_ajuste = sum(cob6$concordante & cob6$estimate < 0 & cob6$lado_ajustado == "subestima", na.rm = TRUE)
)
write_tsv(p6_res, file.path(out_dir, "p6_resumen.tsv"))
print(p6_res %>% mutate(across(where(is.double), ~ signif(.x, 3))), width = Inf)
p6_quintil <- cob6 %>%
  group_by(quintil) %>%
  summarise(rango = paste0(min(n_unidos), "-", max(n_unidos)), n = n(), frac_pasa_hek = mean(pasa_hek), .groups = "drop")
write_tsv(p6_quintil, file.path(out_dir, "p6_paso_por_quintil.tsv"))
print(p6_quintil)

p6 <- cob6 %>%
  mutate(grupo = if_else(pasa_hek, "Pasa el filtro HEK293 (263)", "No pasa (686)")) %>%
  ggplot(aes(n_unidos, fill = grupo)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  scale_fill_manual(values = c("#16a085", "grey60")) +
  labs(x = "Promotores unidos con todas las líneas (log10)", y = "Densidad", fill = NULL) +
  tema +
  theme(legend.position = "bottom")
save_jpg_pdf(p6, "p6_cobertura_filtro_hek", 7, 5)

message("Salidas en ", out_dir)
