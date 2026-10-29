# Exploratorio complementario de remap_ruido_ajustado.R (mismo pedido,
# puntos 5 y 6), en la misma carpeta de salida.
#
# p5 - Ocupacion cuantitativa y ruido. Numero de factores distintos de
#   ReMap2022 con al menos un pico sobre cada promotor, contado sobre el
#   catalogo no redundante completo (remap_tf_hits.tsv, sin marcas de
#   histona - no hay ninguna) y sobre los 982 factores del analisis
#   principal (remap_tf_noise_auc.tsv). "Sin modulo cis-regulatorio" del
#   panorama (N_TF_CRM == 0) sale del mismo bed de ReMap, asi que deberia
#   ser exactamente n_full == 0 (se verifica). Relacion con var_rank_sw:
#   Spearman por replica, media por bin (0 + deciles de los unidos) con IC
#   95% (t), cruda y ajustada (residuo de var_rank_sw ~ log10(hek_tpm + 1) +
#   mean_dnase + CGI_frag, sumado a la media). Modelos lineales con
#   log2(n + 1), la binaria y ambas (con y sin covariables), F parcial para
#   ver si el conteo agrega informacion sobre la binaria.
#
# p6 - Efecto sobre la actividad vs. efecto sobre el ruido, por factor:
#   HL de la union sobre la actividad media (remap_tf_activity_wilcoxon.tsv,
#   remap_tf_activity_gsea) vs. AUC del ruido (remap_tf_noise_auc.tsv). Ambas
#   como percentil dentro de la replica. Ajuste robusto (MASS::rlm) del
#   percentil de ruido sobre el de actividad; residuo estandarizado por su
#   MAD. Desviados = |residuo estandarizado| >= 2 en ambas replicas con el
#   mismo signo. Spearman y Spearman parcial controlando por el numero de
#   promotores unidos (correlacion de los residuos de los rangos).
#
# Requiere las tablas de remap_tf_noise_gsea y remap_tf_activity_gsea, y
# data/processed/{remap_tf_hits,activity_stats_highconf,prom_df}.tsv. Corre
# en ~1 min. Run from the TesisDoc repo root.

library(tidyverse)
library(ggrepel)
source("R/functions/plot_helpers.R")
source("R/functions/fig_paths.R")

out_dir <- "figures/exploratorio_remap_ruido_ajustado"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

auc_orig <- read_tsv(file.path(fig_dir("remap_tf_noise_gsea"), "remap_tf_noise_auc.tsv"), show_col_types = FALSE)
tfs_982 <- unique(auc_orig$feature)

remap_hits <- read_tsv("data/processed/remap_tf_hits.tsv", show_col_types = FALSE) %>% distinct(seq_id, TF)
stopifnot(!any(str_detect(remap_hits$TF, "^H[0-9]")))
stats_highconf <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE)
prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")

# --- p5: ocupacion cuantitativa ---------------------------------------------

ocupacion <- remap_hits %>%
  group_by(seq_id) %>%
  summarise(n_full = n_distinct(TF), n_982 = n_distinct(TF[TF %in% tfs_982]))

data <- inner_join(stats_highconf, prom_df, by = c("seq_id", "name")) %>%
  add_noise_rank() %>%
  left_join(ocupacion, by = "seq_id") %>%
  mutate(
    across(c(n_full, n_982), ~ replace_na(.x, 0L)),
    log_hek = log10(hek_tpm + 1),
    cgi = as.numeric(CGI_frag),
    sin_crm = as.numeric(N_TF_CRM == 0)
  ) %>%
  select(seq_id, rep, var_rank_sw, n_full, n_982, sin_crm, log_hek, mean_dnase, cgi)
stopifnot(!anyNA(data))

# comparacion de las dos versiones del conteo y con la binaria
versiones <- data %>%
  group_by(rep) %>%
  summarise(
    n_prom = n(),
    rho_full_982 = cor(n_full, n_982, method = "spearman"),
    mediana_full = median(n_full), mediana_982 = median(n_982),
    max_full = max(n_full), max_982 = max(n_982),
    prop_cero_full = mean(n_full == 0), prop_cero_982 = mean(n_982 == 0),
    sin_crm_igual_a_cero_full = all((n_full == 0) == (sin_crm == 1)),
    rho_ruido_full = cor(n_full, var_rank_sw, method = "spearman"),
    rho_ruido_982 = cor(n_982, var_rank_sw, method = "spearman"),
    rho_ruido_full_solo_unidos = cor(n_full[n_full > 0], var_rank_sw[n_full > 0], method = "spearman"),
    rho_ocupacion_hek = cor(n_full, log_hek, method = "spearman"),
    rho_ocupacion_dnase = cor(n_full, mean_dnase, method = "spearman"),
    rho_ocupacion_cgi = cor(n_full, cgi, method = "spearman")
  )
write_tsv(versiones, file.path(out_dir, "p5_ocupacion_versiones.tsv"))
print(as.data.frame(versiones))

# var_rank_sw ajustado por covariables: residuo + media, dentro de cada replica
data <- data %>%
  group_by(rep) %>%
  mutate(var_rank_aj = residuals(lm(var_rank_sw ~ log_hek + mean_dnase + cgi)) + mean(var_rank_sw)) %>%
  ungroup()

# rho parcial: Spearman de los residuos de los rangos sobre las covariables
rho_parcial <- data %>%
  group_by(rep) %>%
  summarise(across(c(n_full, n_982), function(n) {
    cov <- cbind(rank(log_hek), rank(mean_dnase), cgi)
    cor(residuals(lm(rank(n) ~ cov)), residuals(lm(rank(var_rank_sw) ~ cov)))
  }, .names = "rho_parcial_{.col}"))
write_tsv(rho_parcial, file.path(out_dir, "p5_rho_parcial_covariables.tsv"))
print(rho_parcial)

# bins: 0 factores + deciles entre los promotores unidos (cortes comunes a
# las dos replicas)
bins_por_version <- function(var) {
  unidos <- data[[var]][data[[var]] > 0]
  cortes <- unique(quantile(unidos, seq(0, 1, 0.1)))
  data %>%
    mutate(
      n = .data[[var]],
      bin = ifelse(n == 0, "0", as.character(cut(n, cortes, include.lowest = TRUE, dig.lab = 4)))
    ) %>%
    pivot_longer(c(var_rank_sw, var_rank_aj), names_to = "medida", values_to = "y") %>%
    group_by(rep, medida, bin) %>%
    summarise(
      n_prom = n(), n_mediano = median(n), media = mean(y),
      ic_low = media - qt(0.975, n_prom - 1) * sd(y) / sqrt(n_prom),
      ic_high = media + qt(0.975, n_prom - 1) * sd(y) / sqrt(n_prom),
      .groups = "drop"
    ) %>%
    mutate(version = var)
}
bins <- bind_rows(bins_por_version("n_full"), bins_por_version("n_982"))
write_tsv(bins, file.path(out_dir, "p5_ruido_por_bin_ocupacion.tsv"))

p_bins <- bins %>%
  mutate(
    medida = factor(medida, levels = c("var_rank_sw", "var_rank_aj"), labels = c("Crudo", "Ajustado por actividad endógena,\naccesibilidad e isla CpG")),
    version = factor(version, levels = c("n_full", "n_982"), labels = c("Catálogo completo (1205 factores)", "982 factores del análisis principal"))
  ) %>%
  ggplot(aes(n_mediano + 1, media, col = rep, shape = medida, linetype = medida)) +
  geom_hline(yintercept = 50.5, col = "grey70") +
  geom_errorbar(aes(ymin = ic_low, ymax = ic_high), width = 0.04) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_x_log10() +
  scale_color_manual(values = c("Rep 1" = "black", "Rep 2" = "grey55")) +
  facet_wrap(~version) +
  ggpubr::theme_pubr(base_size = 14) +
  theme(legend.position = "top", legend.direction = "vertical", strip.background = element_blank()) +
  labs(
    x = "Factores de ReMap unidos al promotor + 1 (mediana del bin, escala log)", y = "Rango de varianza medio (IC 95%)",
    col = NULL, shape = NULL, linetype = NULL, caption = "Bins: 0 factores y deciles de los promotores unidos"
  )
ggsave(file.path(out_dir, "p5_ruido_por_bin_ocupacion.jpg"), p_bins, width = 13, height = 6.5, units = "in")

# modelos: binaria vs. continua (log2(n + 1), catalogo completo), con y sin
# covariables; F parcial de agregar cada una sobre la otra
modelos_p5 <- map(unique(data$rep), function(r) {
  d <- data %>% filter(rep == r) %>% mutate(log_n = log2(n_full + 1))
  map(c(sin_cov = "", con_cov = " + log_hek + mean_dnase + cgi"), function(cv) {
    f <- function(rhs) lm(as.formula(paste("var_rank_sw ~", rhs, cv)), data = d)
    m_bin <- f("sin_crm")
    m_cont <- f("log_n")
    m_ambos <- f("sin_crm + log_n")
    co <- summary(m_ambos)$coefficients
    ci <- confint(m_ambos)
    tibble(
      rep = r, covariables = ifelse(cv == "", "no", "si"),
      coef_sin_crm_solo = coef(m_bin)["sin_crm"], coef_log_n_solo = coef(m_cont)["log_n"],
      coef_sin_crm_ambos = co["sin_crm", 1], ic_sin_crm_ambos = sprintf("[%.2f, %.2f]", ci["sin_crm", 1], ci["sin_crm", 2]), p_sin_crm_ambos = co["sin_crm", 4],
      coef_log_n_ambos = co["log_n", 1], ic_log_n_ambos = sprintf("[%.2f, %.2f]", ci["log_n", 1], ci["log_n", 2]), p_log_n_ambos = co["log_n", 4],
      r2_bin = summary(m_bin)$r.squared, r2_cont = summary(m_cont)$r.squared, r2_ambos = summary(m_ambos)$r.squared,
      F_p_conteo_sobre_binaria = anova(m_bin, m_ambos)$`Pr(>F)`[2],
      F_p_binaria_sobre_conteo = anova(m_cont, m_ambos)$`Pr(>F)`[2]
    )
  }) %>% list_rbind()
}) %>% list_rbind()
write_tsv(modelos_p5, file.path(out_dir, "p5_modelos_binaria_vs_conteo.tsv"))
print(as.data.frame(modelos_p5), digits = 3)

# --- p6: efecto sobre la actividad vs. efecto sobre el ruido ----------------

act <- read_tsv(file.path(fig_dir("remap_tf_activity_gsea"), "remap_tf_activity_wilcoxon.tsv"), show_col_types = FALSE) %>%
  filter(val == "estimate") %>%
  select(TF = feature, rep, HL_actividad = estimate)
n_unidos <- data %>%
  select(seq_id, rep) %>%
  inner_join(remap_hits %>% filter(TF %in% tfs_982), by = "seq_id", relationship = "many-to-many") %>%
  count(TF, rep, name = "n_unidos")
comp <- auc_orig %>%
  select(TF = feature, rep, AUC_ruido = AUC) %>%
  inner_join(act, by = c("TF", "rep")) %>%
  left_join(n_unidos, by = c("TF", "rep")) %>%
  group_by(rep) %>%
  mutate(
    pct_actividad = percent_rank(HL_actividad),
    pct_ruido = percent_rank(AUC_ruido),
    pct_n = percent_rank(n_unidos)
  ) %>%
  ungroup()
message(n_distinct(comp$TF), " factores con efecto sobre actividad y ruido")

# grupos destacados
mll <- read_tsv(file.path(out_dir, "p4_mll1_2_miembros.tsv"), show_col_types = FALSE) %>% pull(TF) %>% unique()
# Polycomb: los que remap_tf_activity_gsea marca en actividad baja (RING1,
# JARID2, CBX8, PHF19) + el resto de los componentes canonicos de PRC1 y
# PRC2 presentes. No se usa el GO PcG protein complex entero porque incluye
# YY1, HDAC2, BAP1, ASXL y otros que no son Polycomb en sentido estricto.
polycomb <- intersect(c(
  "RING1", "JARID2", "CBX8", "PHF19",
  "RNF2", "PCGF1", "PCGF2", "BMI1", "CBX2", "CBX4", "CBX7", "RYBP",
  "EZH1", "EZH2", "SUZ12", "EED", "MTF2"
), unique(comp$TF))
# RNF2 esta en los dos GO; se lo cuenta como Polycomb (es RING1B, PRC1)
grupo_de <- function(tf) case_when(tf == "NFYA" ~ "NF-YA", tf %in% polycomb ~ "Polycomb", tf %in% mll ~ "MLL1/2", TRUE ~ NA_character_)

ajustes <- comp %>%
  group_by(rep) %>%
  group_modify(function(d, k) {
    fit <- MASS::rlm(pct_ruido ~ pct_actividad, data = d, maxit = 100)
    res <- residuals(fit)
    # segunda version: residuo controlando tambien por la cobertura
    fit_n <- MASS::rlm(pct_ruido ~ pct_actividad + pct_n, data = d, maxit = 100)
    res_n <- residuals(fit_n)
    d %>% mutate(
      resid = res, resid_z = res / mad(res), intercepto = coef(fit)[1], pendiente = coef(fit)[2],
      resid_z_cob = res_n / mad(res_n)
    )
  }) %>%
  ungroup()

cors <- comp %>%
  group_by(rep) %>%
  summarise(
    n = n(),
    rho = cor(HL_actividad, AUC_ruido, method = "spearman"),
    rho_parcial_n = cor(residuals(lm(rank(HL_actividad) ~ rank(n_unidos))), residuals(lm(rank(AUC_ruido) ~ rank(n_unidos)))),
    rho_actividad_n = cor(HL_actividad, n_unidos, method = "spearman"),
    rho_ruido_n = cor(AUC_ruido, n_unidos, method = "spearman")
  ) %>%
  left_join(ajustes %>% distinct(rep, intercepto, pendiente), by = "rep")
write_tsv(cors, file.path(out_dir, "p6_correlaciones.tsv"))
print(as.data.frame(cors), digits = 3)

tabla <- ajustes %>%
  select(TF, rep, n_unidos, pct_n, HL_actividad, AUC_ruido, pct_actividad, pct_ruido, resid_z, resid_z_cob) %>%
  pivot_wider(names_from = rep, values_from = -c(TF, rep), names_glue = "{.value}_{str_remove(rep, ' ')}") %>%
  mutate(
    resid_z_medio = (resid_z_Rep1 + resid_z_Rep2) / 2,
    resid_z_cob_medio = (resid_z_cob_Rep1 + resid_z_cob_Rep2) / 2,
    desvio_cob = case_when(
      resid_z_cob_Rep1 <= -2 & resid_z_cob_Rep2 <= -2 ~ "menos ruido de lo esperado",
      resid_z_cob_Rep1 >= 2 & resid_z_cob_Rep2 >= 2 ~ "más ruido de lo esperado",
      TRUE ~ "dentro de la tendencia"
    ),
    desvio = case_when(
      resid_z_Rep1 <= -2 & resid_z_Rep2 <= -2 ~ "menos ruido de lo esperado",
      resid_z_Rep1 >= 2 & resid_z_Rep2 >= 2 ~ "más ruido de lo esperado",
      TRUE ~ "dentro de la tendencia"
    ),
    grupo = grupo_de(TF)
  ) %>%
  arrange(resid_z_medio)
write_tsv(tabla, file.path(out_dir, "p6_residuos_por_factor.tsv"))
# Con |z| >= 2 en ambas replicas puede no haber ninguno: ademas, los 15 de
# cada extremo por residuo medio entre los que tienen el mismo signo en ambas
top_desv <- tabla %>%
  filter(sign(resid_z_Rep1) == sign(resid_z_Rep2)) %>%
  mutate(lado = ifelse(resid_z_medio < 0, "menos ruido de lo esperado", "más ruido de lo esperado")) %>%
  group_by(lado) %>%
  slice_max(abs(resid_z_medio), n = 15) %>%
  ungroup() %>%
  select(lado, TF, grupo, n_unidos_Rep1, n_unidos_Rep2, pct_n_Rep1, resid_z_Rep1, resid_z_Rep2, pct_actividad_Rep1, pct_ruido_Rep1, pct_actividad_Rep2, pct_ruido_Rep2)
write_tsv(top_desv, file.path(out_dir, "p6_top_desviados_consistentes.tsv"))
print(as.data.frame(top_desv), digits = 2)
top_desv_cob <- tabla %>%
  filter(sign(resid_z_cob_Rep1) == sign(resid_z_cob_Rep2)) %>%
  mutate(lado = ifelse(resid_z_cob_medio < 0, "menos ruido de lo esperado", "más ruido de lo esperado")) %>%
  group_by(lado) %>%
  slice_max(abs(resid_z_cob_medio), n = 15) %>%
  ungroup() %>%
  select(lado, TF, grupo, n_unidos_Rep1, n_unidos_Rep2, pct_n_Rep1, resid_z_cob_Rep1, resid_z_cob_Rep2, resid_z_Rep1, resid_z_Rep2)
write_tsv(top_desv_cob, file.path(out_dir, "p6_top_desviados_consistentes_ajustado_cobertura.tsv"))
print(as.data.frame(top_desv_cob), digits = 2)
print(tabla %>% count(desvio_cob))
print(top_desv_cob %>% group_by(lado) %>% summarise(pct_n_mediano = median(pct_n_Rep1)))
print(tabla %>% filter(!is.na(grupo)) %>% group_by(grupo) %>% summarise(n = n(), z_medio = median(resid_z_medio), z_cob_medio = median(resid_z_cob_medio)))
print(top_desv %>% group_by(lado) %>% summarise(pct_n_mediano = median(pct_n_Rep1), n_min = min(n_unidos_Rep1), n_max = max(n_unidos_Rep1)))
tabla %>%
  filter(desvio != "dentro de la tendencia" | !is.na(grupo)) %>%
  write_tsv(file.path(out_dir, "p6_desviados_y_destacados.tsv"))
print(tabla %>% count(desvio))
print(tabla %>% filter(desvio != "dentro de la tendencia") %>% select(TF, desvio, n_unidos_Rep1, pct_n_Rep1, resid_z_Rep1, resid_z_Rep2, grupo) %>% as.data.frame(), digits = 2)
print(tabla %>% filter(!is.na(grupo)) %>% select(TF, grupo, n_unidos_Rep1, resid_z_Rep1, resid_z_Rep2, pct_actividad_Rep1, pct_ruido_Rep1) %>% arrange(grupo) %>% as.data.frame(), digits = 2)
# cobertura de los desviados frente al resto
print(tabla %>% group_by(desvio) %>% summarise(n = n(), pct_n_mediano = median((pct_n_Rep1 + pct_n_Rep2) / 2), n_unidos_mediano = median(n_unidos_Rep1)))

plot_df <- ajustes %>%
  left_join(tabla %>% select(TF, desvio, grupo), by = "TF")
# percentil en el que cae el efecto nulo (HL = 0, AUC = 0.5) en cada replica
ceros <- plot_df %>%
  group_by(rep) %>%
  summarise(x0 = mean(HL_actividad < 0), y0 = mean(AUC_ruido < 0.5))
p_comp <- plot_df %>%
  ggplot(aes(pct_actividad, pct_ruido)) +
  geom_vline(data = ceros, aes(xintercept = x0), col = "grey60", linetype = "dotted", linewidth = 0.7) +
  geom_hline(data = ceros, aes(yintercept = y0), col = "grey60", linetype = "dotted", linewidth = 0.7) +
  geom_abline(data = cors, aes(slope = pendiente, intercept = intercepto), col = "grey40", linetype = "dashed") +
  geom_point(aes(col = log10(n_unidos)), size = 1.4, alpha = 0.8) +
  geom_point(data = ~ filter(.x, desvio != "dentro de la tendencia"), shape = 1, size = 3.2, col = "black") +
  geom_point(data = ~ filter(.x, !is.na(grupo)), aes(fill = grupo), shape = 23, size = 3, col = "black") +
  geom_text_repel(data = ~ filter(.x, !is.na(grupo) | desvio != "dentro de la tendencia"), aes(label = TF), size = 2.8, max.overlaps = 40, min.segment.length = 0) +
  geom_text(data = cors, aes(x = 0, y = 1.04, label = sprintf("rho = %.2f; rho parcial (n unidos) = %.2f", rho, rho_parcial_n)), hjust = 0, size = 4) +
  scale_color_viridis_c(name = "log10(promotores unidos)") +
  scale_fill_manual(values = c("MLL1/2" = "#AD343E", "Polycomb" = "#1B8C8E", "NF-YA" = "#FFB400"), name = NULL) +
  facet_wrap(~rep) +
  coord_equal() +
  ggpubr::theme_pubr(base_size = 13) +
  theme(legend.position = "top", strip.background = element_blank()) +
  labs(
    x = "Efecto sobre la actividad (percentil del HL; alto = más actividad)",
    y = "Efecto sobre el ruido (percentil de la AUC; alto = más ruido)",
    caption = "Círculo negro: residuo robusto |z| >= 2 con el mismo signo en ambas réplicas. Línea discontinua: ajuste robusto (rlm). Líneas punteadas: efecto nulo (HL = 0; AUC = 0.5)."
  )
ggsave(file.path(out_dir, "p6_actividad_vs_ruido.jpg"), p_comp, width = 15, height = 8.5, units = "in")

# residuo vs. cobertura: si los desviados son todos de cobertura extrema
p_cob <- tabla %>%
  pivot_longer(c(resid_z_Rep1, resid_z_Rep2), names_to = "rep", values_to = "resid_z") %>%
  mutate(rep = str_replace(str_remove(rep, "resid_z_"), "Rep", "Rep "), n_unidos = ifelse(rep == "Rep 1", n_unidos_Rep1, n_unidos_Rep2)) %>%
  ggplot(aes(n_unidos, resid_z)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", col = "grey50") +
  geom_point(aes(col = desvio), alpha = 0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("dentro de la tendencia" = "grey70", "menos ruido de lo esperado" = "#7FB800", "más ruido de lo esperado" = "#D6741F"), name = NULL) +
  facet_wrap(~rep) +
  ggpubr::theme_pubr(base_size = 13) +
  theme(strip.background = element_blank()) +
  labs(x = "Promotores unidos (escala log)", y = "Residuo robusto estandarizado")
ggsave(file.path(out_dir, "p6_residuo_vs_cobertura.jpg"), p_cob, width = 12, height = 5.5, units = "in")

message("Exploratorio guardado en ", out_dir)
