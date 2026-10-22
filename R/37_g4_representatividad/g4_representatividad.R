# EXPLORATORIO - ¿en que paso de la medicion se pierden los promotores con
# G4 en la hebra codificante, y depende del contenido de G+C?
#
# En R/36_polyAT_g4 los G4 en la hebra codificante (la del promotor) son mas
# frecuentes que en la molde en la library disenada, pero la relacion se
# invierte entre los promotores de actividad de alta confianza. Aca se
# sigue la representacion de cada clase de G4 a lo largo de las etapas:
#   0. Library disenada (20.851 promotores).
#   1. Presort 2022: lecturas > 0 en el pool pre-sorting del reportero de
#      splicing (mismo pool de oligos, otro vector; data/external/
#      2022_presort_counts.tsv, el mismo insumo que sesgo_representatividad).
#      Es la unica medida previa al sorting; no es un paso anidado de esta
#      medicion sino un control del pool independiente de la actividad.
#   2. Detectado en el sorting: lecturas > 0 en alguna compuerta (por replica).
#   3. Cuentas minimas: sum_counts >= 100 y sum_norm_counts > 1000
#      (activity_stats.tsv; por replica).
#   4. Mediana consistente entre replicas (activity_stats_consistent.tsv).
#   5. Alta confianza: ademas no bimodal (activity_stats_highconf.tsv).
# Las etapas 2-5 son anidadas dentro de cada replica (R/01_activity_stats).
#
# G4: regex canonico (clases de figures/exploratorio_polyAT_g4/
# polyAT_g4_por_promotor.tsv) y G4-seq K+ / K+ + PDS por hebra
# (data/processed/g4seq_lib_v1.tsv, solo liftover_ok). Clases exclusivas:
# Sin G4, Solo codificante, Solo molde, Ambas.
#
# Control de G+C: la representacion depende fuertemente del G+C del
# fragmento (sesgo_representatividad: casi no se recuperan los fragmentos
# con G+C < 0,40) y, dentro de un tercil de G+C, los Sin G4 son mas AT que
# los que tienen G4 - los terciles no alcanzan (la referencia queda
# castigada). Por eso todo se controla con bins de G+C de ancho fijo
# (0,05). Todos los promotores juntos (no se estratifica por isla CpG):
#   - las comparaciones por bin se hacen dentro de cada bin;
#   - las comparaciones globales contra Sin G4 usan una referencia apareada
#     por bin (por cada promotor de la clase, 3 Sin G4 del mismo bin de G+C,
#     sorteados con semilla fija; con reposicion solo si el bin no tiene
#     suficientes);
#   - la retencion global se testea con Cochran-Mantel-Haenszel
#     estratificado por bin (OR comun) y se informa la diferencia contra la
#     retencion de Sin G4 estandarizada a la distribucion de G+C de la clase;
#   - Solo codificante vs. Solo molde se compara directo (tienen la misma
#     distribucion de G+C; se informan sus medianas).
#
# Analisis:
#   A. Retencion desde la library por clase x bin de G+C x etapa (grafico
#      principal) y por paso (CMH), mas la prevalencia de G4 por hebra en
#      cada etapa, cruda y estandarizada por G+C a la library.
#   B. Abundancia de lecturas: log10(lecturas + 1) en el presort 2022 y en
#      el sorting (suma de compuertas, por replica, toda la library); HL
#      (IC 95%) y p de Wilcoxon por bin y global.
#   C. ¿Puede la perdida sesgar el efecto de los G4 sobre la actividad? HL
#      (IC 95%) de la actividad media en: alta confianza (referencia);
#      cuentas minimas sin filtros de consistencia; solo promotores con
#      >= 100 lecturas en el presort 2022 (umbral fijo); y un peor caso en
#      el que los promotores que cada clase pierde "de mas" respecto de Sin
#      G4 del mismo bin de G+C se reincorporan con la actividad minima (1).
# Cada replica por separado; BH dentro de metodo x replica (o conjunto) x
# analisis (y escenario en C).
#
# Salidas: la figura g4_representatividad (retencion por G+C, jpg + pdf)
# va a fig_dir("g4_representatividad") (ver docs/mapping_figuras.csv); el
# resto (tablas y graficos de soporte, exploratorios) a
# figures/exploratorio_g4_representatividad/.
# Requiere: figures/exploratorio_polyAT_g4/polyAT_g4_por_promotor.tsv
# (R/36_polyAT_g4/polyAT_g4.R), data/processed/g4seq_lib_v1.tsv,
# data/processed/data_long.tsv y activity_stats{,_consistent,_highconf}.tsv
# (R/01_activity_stats) y data/external/2022_presort_counts.tsv.
# Run from the TesisDoc repo root.

library(tidyverse)
library(patchwork)
source("R/functions/tata_scan.R") # tata_hl_vs()
source("R/functions/fig_paths.R") # fig_dir()

out_dir <- "figures/exploratorio_g4_representatividad"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261019)

MIN_N <- 10 # minimo por grupo para HL / Fisher
MIN_PLOT <- 20 # minimo de promotores de la clase en el bin para graficarlo
MIN_PRESORT <- 100
RATIO_REF <- 3
REPS <- c("Rep 1", "Rep 2")
METODOS <- c("Regex canónico", "G4-seq K+", "G4-seq K+ + PDS")
METODOS_PLOT <- c("G4-seq K+", "Regex canónico")
G4_CLASES <- c("Sin G4", "Solo codificante", "Solo molde", "Ambas")
CLR_CLASE <- c("Sin G4" = "grey45", "Solo codificante" = "#C0392B", "Solo molde" = "#2471A3")
ETAPAS <- c(
  diseno = "Library diseñada", presort = "Presort 2022", detectado = "Detectado en el sorting",
  conteos = "Cuentas mínimas", consistente = "Consistente entre réplicas", highconf = "Alta confianza"
)
GC_BREAKS <- seq(0.20, 1.00, by = 0.05)
fmt2 <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = ",")
GC_LABS <- paste0(fmt2(head(GC_BREAKS, -1)), "–", fmt2(tail(GC_BREAKS, -1)))
GC_MID <- head(GC_BREAKS, -1) + 0.025

COMPS <- tribble(
  ~comparacion, ~clase, ~ref,
  "Solo codificante vs. Sin G4", "Solo codificante", "Sin G4",
  "Solo molde vs. Sin G4", "Solo molde", "Sin G4",
  "Solo codificante vs. Solo molde", "Solo codificante", "Solo molde"
)
CLR_COMP <- setNames(c("#C0392B", "#2471A3", "#7D3C98"), COMPS$comparacion)

save_tsv <- function(df, file) write_tsv(df, file.path(out_dir, file))
efecto_lab <- function(est, p_adj, pos = "Mayor", neg = "Menor") {
  case_when(p_adj < 0.05 & est > 0 ~ pos, p_adj < 0.05 & est < 0 ~ neg, TRUE ~ "n.s.")
}

# --- Promotores, bins de G+C y clases de G4 ------------------------------------

por_prom <- read_tsv("figures/exploratorio_polyAT_g4/polyAT_g4_por_promotor.tsv", show_col_types = FALSE)
prom <- por_prom %>%
  transmute(
    seq_id, gc_frag,
    gc_i = findInterval(gc_frag, GC_BREAKS, rightmost.closed = TRUE),
    gc_bin = factor(GC_LABS[gc_i], GC_LABS), gc_mid = GC_MID[gc_i],
    mbin = as.character(gc_bin) # estrato de apareamiento
  ) %>%
  select(-gc_i)

g4 <- bind_rows(
  por_prom %>% transmute(seq_id, metodo = "Regex canónico", cod = `n_g4_cod__Regex canónico` > 0, molde = `n_g4_molde__Regex canónico` > 0),
  read_tsv("data/processed/g4seq_lib_v1.tsv", show_col_types = FALSE) %>%
    filter(liftover_ok) %>%
    transmute(seq_id, metodo = paste("G4-seq", condicion), cod = n_oq_cod > 0, molde = n_oq_molde > 0)
) %>%
  mutate(
    clase = case_when(cod & molde ~ "Ambas", cod ~ "Solo codificante", molde ~ "Solo molde", TRUE ~ "Sin G4") %>% factor(G4_CLASES),
    metodo = factor(metodo, METODOS)
  )
stopifnot(all(g4$clase[g4$metodo == "Regex canónico"] == por_prom$`g4_clase__Regex canónico`))

# referencia apareada por bin de G+C: RATIO_REF promotores de `pool`
# por cada uno de `casos` del mismo mbin (con reposicion solo si faltan);
# los casos de un mbin sin referencia se descartan
aparear <- function(casos, pool) {
  pool_split <- split(pool, pool$mbin)
  casos <- casos %>% filter(mbin %in% names(pool_split))
  ref <- casos %>%
    count(mbin) %>%
    pmap_dfr(function(mbin, n) {
      p <- pool_split[[mbin]]
      k <- RATIO_REF * n
      p[sample(nrow(p), k, replace = nrow(p) < k), ]
    })
  list(casos = casos, ref = ref)
}

# HL de `y` para las 3 comparaciones: vs. Sin G4 con referencia apareada por
# G+C; codificante vs. molde directo
comparar_hl <- function(d, y) {
  pmap_dfr(COMPS, function(comparacion, clase, ref) {
    x <- d[d$clase == clase, ]
    r <- d[d$clase == ref, ]
    if (ref == "Sin G4") {
      m <- aparear(x, r)
      x <- m$casos
      r <- m$ref
    }
    tata_hl_vs(x[[y]], r[[y]], MIN_N) %>%
      mutate(comparacion = comparacion, n = nrow(x), n_ref = nrow(r), gc_mediana = median(x$gc_frag), gc_mediana_ref = median(r$gc_frag))
  })
}

# --- Etapas -------------------------------------------------------------------

presort <- read_tsv("data/external/2022_presort_counts.tsv", show_col_types = FALSE) %>%
  transmute(seq_id = seqname, presort_count = original_count)
sort_counts <- read_tsv("data/processed/data_long.tsv", show_col_types = FALSE) %>%
  group_by(seq_id, rep) %>%
  summarise(sort_count = sum(counts), .groups = "drop")
ids_de <- function(file) {
  read_tsv(file.path("data/processed", file), show_col_types = FALSE) %>%
    distinct(seq_id, rep) %>%
    mutate(key = paste(seq_id, rep)) %>%
    pull(key)
}
k_conteos <- ids_de("activity_stats.tsv")
k_consist <- ids_de("activity_stats_consistent.tsv")
k_highconf <- ids_de("activity_stats_highconf.tsv")

etapas <- expand_grid(seq_id = prom$seq_id, rep = REPS) %>%
  left_join(presort, by = "seq_id") %>%
  left_join(sort_counts, by = c("seq_id", "rep")) %>%
  mutate(
    across(c(presort_count, sort_count), ~ replace_na(.x, 0)),
    key = paste(seq_id, rep),
    diseno = TRUE,
    presort = presort_count > 0,
    detectado = sort_count > 0,
    conteos = key %in% k_conteos,
    consistente = key %in% k_consist,
    highconf = key %in% k_highconf
  ) %>%
  select(-key)
stopifnot(with(etapas, all(!conteos | detectado), all(!consistente | conteos), all(!highconf | consistente)))
print(etapas %>% group_by(rep) %>% summarise(across(all_of(names(ETAPAS)), sum)))

base <- etapas %>%
  inner_join(prom, by = "seq_id") %>%
  inner_join(g4, by = "seq_id", relationship = "many-to-many")

# --- A1. Retencion desde la library por bin de G+C -----------------------------

retencion_bin <- map_dfr(names(ETAPAS)[-1], function(e) {
  base %>%
    group_by(metodo, rep, gc_bin, gc_mid, clase) %>%
    summarise(n = n(), n_pasa = sum(.data[[e]]), .groups = "drop") %>%
    mutate(etapa = ETAPAS[[e]])
}) %>%
  mutate(pct_pasa = 100 * n_pasa / n, etapa = factor(etapa, ETAPAS)) %>%
  relocate(metodo, rep, etapa, gc_bin)
save_tsv(retencion_bin, "g4_retencion_por_bin_gc.tsv")

# --- A2. Retencion por paso, controlada por G+C (CMH) ------------------------------

PASOS <- tribble(
  ~desde, ~hasta,
  "diseno", "presort",
  "diseno", "detectado",
  "detectado", "conteos",
  "conteos", "consistente",
  "consistente", "highconf",
  "diseno", "highconf"
) %>% mutate(paso = paste(ETAPAS[desde], "→", ETAPAS[hasta]))

# OR comun de pasar (clase vs. ref) estratificado por mbin; estratos con
# ambos grupos presentes. Ademas la retencion de la ref estandarizada a la
# distribucion de bins de la clase.
cmh_retencion <- function(x, clase, ref) {
  d <- x %>% filter(clase %in% c(.env$clase, .env$ref))
  ok <- d %>%
    group_by(mbin) %>%
    summarise(ok = n_distinct(clase) == 2, .groups = "drop") %>%
    filter(ok) %>%
    pull(mbin)
  d <- d %>% filter(mbin %in% ok)
  a <- d %>% filter(clase == .env$clase)
  r <- d %>% filter(clase == .env$ref)
  out <- tibble(
    n = nrow(a), n_pasa = sum(a$pasa), n_ref = nrow(r), n_ref_pasa = sum(r$pasa),
    pct_ref_estandarizado = 100 * sum(table(factor(a$mbin, ok)) * tapply(r$pasa, factor(r$mbin, ok), mean)) / nrow(a)
  )
  if (out$n < MIN_N || out$n_ref < MIN_N) {
    return(out %>% mutate(or_comun = NA_real_, or_low = NA_real_, or_high = NA_real_, p = NA_real_, n_bins = length(ok)))
  }
  tab <- table(factor(d$clase == clase, c(TRUE, FALSE)), factor(d$pasa, c(TRUE, FALSE)), d$mbin)
  t <- if (dim(tab)[3] >= 2) mantelhaen.test(tab, exact = FALSE) else fisher.test(tab[, , 1])
  out %>% mutate(or_comun = unname(t$estimate), or_low = t$conf.int[1], or_high = t$conf.int[2], p = t$p.value, n_bins = length(ok))
}

retencion_paso <- pmap_dfr(PASOS, function(desde, hasta, paso) {
  base %>%
    filter(.data[[desde]], hasta != "presort" | rep == REPS[1]) %>% # el presort no depende de la replica
    mutate(pasa = .data[[hasta]]) %>%
    group_by(metodo, rep) %>%
    group_modify(function(x, k) pmap_dfr(COMPS, function(comparacion, clase, ref) cmh_retencion(x, clase, ref) %>% mutate(comparacion = comparacion))) %>%
    ungroup() %>%
    mutate(paso = paso, rep = if (hasta == "presort") "Presort 2022" else rep)
}) %>%
  mutate(
    pct_pasa = 100 * n_pasa / n, pct_ref_crudo = 100 * n_ref_pasa / n_ref,
    dif_pct_ajustada = pct_pasa - pct_ref_estandarizado
  ) %>%
  group_by(metodo, rep) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup() %>%
  mutate(retencion = efecto_lab(log(or_comun), p_adj, "Mayor que la ref.", "Menor que la ref.")) %>%
  relocate(metodo, rep, paso, comparacion)
save_tsv(retencion_paso, "g4_retencion_por_paso_cmh.tsv")

# --- A3. Prevalencia por etapa, cruda y estandarizada por G+C ----------------------
# estandarizada = sum_bin w_bin * prevalencia_bin, con w = distribucion de
# bins de G+C de la library disenada; los bins sin
# promotores en la etapa se excluyen y los pesos se renormalizan.

prev_bin <- map_dfr(names(ETAPAS), function(e) {
  base %>%
    filter(.data[[e]]) %>%
    group_by(metodo, rep, mbin) %>%
    summarise(n = n(), n_cod = sum(cod), n_molde = sum(molde), .groups = "drop") %>%
    mutate(etapa = ETAPAS[[e]])
})
pesos <- prev_bin %>%
  filter(etapa == ETAPAS[["diseno"]]) %>%
  group_by(metodo, rep) %>%
  transmute(mbin, w = n / sum(n)) %>%
  ungroup()
prevalencia <- prev_bin %>%
  left_join(pesos, by = c("metodo", "rep", "mbin")) %>%
  group_by(metodo, rep, etapa) %>%
  summarise(
    pct_cod_std = 100 * sum(w * n_cod / n) / sum(w),
    pct_molde_std = 100 * sum(w * n_molde / n) / sum(w),
    n_tot = sum(n), pct_cod = 100 * sum(n_cod) / n_tot, pct_molde = 100 * sum(n_molde) / n_tot,
    .groups = "drop"
  ) %>%
  rename(n = n_tot) %>%
  mutate(etapa = factor(etapa, ETAPAS)) %>%
  select(metodo, rep, etapa, n, pct_cod, pct_molde, pct_cod_std, pct_molde_std) %>%
  arrange(metodo, rep, etapa)
save_tsv(prevalencia, "g4_prevalencia_por_etapa.tsv")

# --- B. Abundancia de lecturas -----------------------------------------------------

abund <- bind_rows(
  base %>% filter(rep == REPS[1]) %>% mutate(conjunto = "Presort 2022", y = log10(presort_count + 1)),
  base %>% mutate(conjunto = rep, y = log10(sort_count + 1))
) %>%
  mutate(conjunto = factor(conjunto, c("Presort 2022", REPS)))

abundancia_bin <- abund %>%
  group_by(metodo, conjunto, gc_bin, gc_mid) %>%
  group_modify(function(d, k) {
    pmap_dfr(COMPS, function(comparacion, clase, ref) {
      x <- d$y[d$clase == clase]
      r <- d$y[d$clase == ref]
      tata_hl_vs(x, r, MIN_N) %>% mutate(comparacion = comparacion, n = length(x), n_ref = length(r))
    })
  }) %>%
  group_by(metodo, conjunto) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup() %>%
  mutate(efecto = efecto_lab(estimate, p_adj))
save_tsv(abundancia_bin, "g4_abundancia_por_bin_gc_hl.tsv")

abundancia_global <- abund %>%
  group_by(metodo, conjunto) %>%
  group_modify(~ comparar_hl(.x, "y")) %>%
  group_by(metodo, conjunto) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup() %>%
  mutate(efecto = efecto_lab(estimate, p_adj), veces = 10^estimate)
save_tsv(abundancia_global, "g4_abundancia_global_hl.tsv")

# --- C. Sesgo potencial sobre los efectos de G4 en la actividad -----------------------

act_base <- function(file) {
  read_tsv(file.path("data/processed", file), show_col_types = FALSE) %>%
    select(seq_id, rep, mean, sum_counts) %>%
    inner_join(prom, by = "seq_id") %>%
    left_join(presort, by = "seq_id") %>%
    mutate(presort_count = replace_na(presort_count, 0)) %>%
    inner_join(g4, by = "seq_id", relationship = "many-to-many")
}
act_conteos <- act_base("activity_stats.tsv")
hc <- act_base("activity_stats_highconf.tsv")

# ¿las lecturas del sorting dependen de la actividad?
cor_lecturas_act <- act_conteos %>%
  filter(metodo == METODOS[1]) %>%
  group_by(rep) %>%
  summarise(
    rho_spearman = cor(log10(sum_counts), mean, method = "spearman"),
    p = cor.test(log10(sum_counts), mean, method = "spearman", exact = FALSE)$p.value,
    n = n(), .groups = "drop"
  )
print(cor_lecturas_act)

# peor caso: por bin de G+C y clase, los promotores perdidos de mas
# respecto de Sin G4 del mismo bin (library -> alta confianza) se agregan
# con actividad = 1
exceso <- retencion_bin %>%
  filter(etapa == ETAPAS[["highconf"]], rep %in% REPS) %>%
  group_by(metodo, rep, gc_bin, gc_mid) %>%
  mutate(ret_ref = sum(n_pasa[clase == "Sin G4"]) / sum(n[clase == "Sin G4"])) %>%
  ungroup() %>%
  filter(clase != "Sin G4", is.finite(ret_ref)) %>%
  transmute(metodo, rep, gc_bin, gc_mid, clase, n_extra = pmax(0, round(n * ret_ref - n_pasa)))
peor_caso <- hc %>%
  bind_rows(
    exceso %>%
      filter(n_extra > 0) %>%
      uncount(n_extra) %>%
      mutate(mean = 1, mbin = as.character(gc_bin), gc_frag = gc_mid)
  )

ESCENARIOS <- c(
  "Alta confianza (referencia)", "Cuentas mínimas (sin filtros de consistencia)",
  "Alta confianza, presort ≥ 100 lecturas", "Peor caso: pérdidas en exceso con actividad mínima"
)
sesgo <- bind_rows(
  hc %>% mutate(escenario = ESCENARIOS[1]),
  act_conteos %>% mutate(escenario = ESCENARIOS[2]),
  hc %>% filter(presort_count >= MIN_PRESORT) %>% mutate(escenario = ESCENARIOS[3]),
  peor_caso %>% mutate(escenario = ESCENARIOS[4])
) %>%
  group_by(metodo, rep, escenario) %>%
  group_modify(~ comparar_hl(.x, "mean")) %>%
  group_by(metodo, rep, escenario) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%
  ungroup() %>%
  mutate(efecto = efecto_lab(estimate, p_adj, "Mayor actividad", "Menor actividad"), escenario = factor(escenario, ESCENARIOS)) %>%
  relocate(metodo, rep, escenario, comparacion)
save_tsv(sesgo, "g4_sesgo_efectos_actividad_hl.tsv")
save_tsv(cor_lecturas_act, "g4_spearman_lecturas_vs_actividad.tsv")

# --- Graficos -------------------------------------------------------------------

# principal: % de la library que llega a cada etapa, por clase y bin de G+C
p_ret <- retencion_bin %>%
  filter(metodo %in% METODOS_PLOT, clase != "Ambas", !etapa %in% ETAPAS[c("consistente", "conteos")]) %>%
  filter(etapa != ETAPAS[["presort"]] | rep == REPS[1]) %>%
  group_by(metodo, rep, etapa, gc_bin, gc_mid, clase) %>%
  summarise(n = sum(n), n_pasa = sum(n_pasa), .groups = "drop") %>%
  filter(n >= MIN_PLOT) %>%
  mutate(pct_pasa = 100 * n_pasa / n, metodo = factor(metodo, METODOS_PLOT)) %>%
  ggplot(aes(gc_mid, pct_pasa, color = clase, linetype = rep, group = interaction(clase, rep))) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(size = n), alpha = 0.8) +
  facet_grid(metodo ~ etapa, labeller = labeller(
    metodo = c("G4-seq K+" = "G4 observados (G4-seq)", "Regex canónico" = "G4 predichos (secuencia)"),
    etapa = ~ if_else(.x == ETAPAS[["presort"]], "Pool previo al sorting", .x)
  )) +
  scale_color_manual(values = CLR_CLASE) +
  scale_size_area(max_size = 5, breaks = c(50, 500, 2000)) +
  scale_x_continuous(labels = ~ fmt2(.x)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = "Contenido de G+C del fragmento (bins de 0,05)", y = "% de la library que llega a la etapa",
    color = NULL, linetype = "Réplica", size = "n en el bin"
  ) +
  guides(color = guide_legend(order = 1, override.aes = list(linewidth = 1.5, size = 4)), linetype = guide_legend(order = 2, override.aes = list(linewidth = 1.1)), size = guide_legend(order = 3)) +
  theme_bw(base_size = 22) +
  theme(legend.position = "top", legend.box = "vertical", legend.key.width = unit(1.6, "cm"), panel.spacing = unit(1, "lines"))
for (ext in c("jpg", "pdf")) {
  ggsave(file.path(fig_dir("g4_representatividad"), paste0("g4_retencion_por_gc.", ext)), p_ret, width = 15, height = 11, units = "in", dpi = 200,
         device = if (ext == "pdf") cairo_pdf else NULL)
}

# abundancia de lecturas por bin (HL con IC)
p_abund <- abundancia_bin %>%
  filter(metodo == "G4-seq K+", !is.na(estimate)) %>%
  ggplot(aes(gc_mid, estimate, color = comparacion, group = comparacion)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high), position = position_dodge(width = 0.02)) +
  geom_point(aes(shape = p_adj < 0.05), position = position_dodge(width = 0.02), size = 2) +
  geom_line(position = position_dodge(width = 0.02), alpha = 0.5) +
  facet_grid(. ~ conjunto) +
  scale_color_manual(values = CLR_COMP) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p ajustado < 0,05", `FALSE` = "n.s.")) +
  scale_x_continuous(labels = ~ fmt2(.x)) +
  labs(x = "Contenido de G+C del fragmento (bins de 0,05)", y = "Hodges-Lehmann de log10(lecturas + 1) (IC 95%)", color = NULL, shape = NULL, title = "G4-seq K+") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "g4_abundancia_por_gc.jpg"), p_abund, width = 13, height = 7.5, units = "in", dpi = 200)

# sesgo: HL de actividad por escenario
p_sesgo <- sesgo %>%
  filter(metodo %in% METODOS_PLOT) %>%
  mutate(comparacion = factor(comparacion, COMPS$comparacion)) %>%
  ggplot(aes(estimate, escenario, color = comparacion, group = comparacion)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), position = position_dodge(width = 0.7), linewidth = 0.8) +
  geom_point(aes(shape = p_adj < 0.05), position = position_dodge(width = 0.7), size = 2.3) +
  facet_grid(metodo ~ rep) +
  scale_color_manual(values = CLR_COMP) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "p ajustado < 0,05", `FALSE` = "n.s.")) +
  scale_y_discrete(limits = rev) +
  labs(
    x = "Hodges-Lehmann de la actividad media (IC 95%); vs. Sin G4 con referencia apareada por G+C",
    y = NULL, color = NULL, shape = NULL
  ) +
  theme_bw() +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "g4_sesgo_efectos_actividad.jpg"), p_sesgo, width = 14, height = 11, units = "in", dpi = 200)

message("Salidas en ", out_dir)
