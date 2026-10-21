# EXPLORATORIO - tractos poli(dA:dT) y G-cuadruplex (G4) en la library y
# actividad.
#
# Parte A - poli(dA:dT): por promotor, el tracto homopolimerico de A o de T
# mas largo (A_n o T_n en la hebra del promotor; es el mismo tracto visto
# desde cualquiera de las dos hebras), la posicion de su centro respecto del
# TSS y el numero de tractos >= 6 y >= 8. Categorias por el tracto mas
# largo: < 6, 6-7, 8-9, >= 10. Efecto de cada categoria vs. "< 6"; y por
# posicion: promotores con algun tracto >= 8 centrado en la region proximal
# (-100 a -1, donde se espera la zona libre de nucleosomas) vs. con tractos
# >= 8 solo mas rio arriba (-236 a -101), cada uno vs. "< 6".
# Si hay empate en el tracto mas largo se toma el mas cercano al TSS. Un
# centro semientero en -100,5 cuenta como rio arriba; tractos centrados en
# -0,5 o despues no entran en ninguna de las dos regiones.
#
# Parte B - G4: patron canonico G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,} en
# la hebra del promotor (G4 en la hebra codificante, sin molde) y su
# equivalente C{3,}N{1,7}C{3,}N{1,7}C{3,}N{1,7}C{3,} (G4 en la hebra molde);
# coincidencias no solapadas (gregexpr). Si esta instalado pqsfinder se usa
# tambien, con el puntaje por defecto (hebra + = codificante, - = molde).
# Ademas, como otros dos "metodos", los G4 observados in vitro por G4-seq
# (Marsico et al. 2019, K+ y K+ + PDS; R/00_prom_features/build_g4seq.R),
# por hebra respecto del promotor; se reporta su concordancia con el regex.
# Clases: sin G4, solo codificante, solo molde, ambas. Prevalencia con /
# sin isla CpG (OR y p de Fisher) y efecto de cada clase (y de "cualquier
# hebra") vs. sin G4.
#
# Parte C - G4 plegado en HEK293T: senal de G4-CUT&Tag (BG4 / IgG,
# GSE178668) por fragmento y ATAC-seq de la misma serie; ¿la senal se
# asocia con el G4 por secuencia o por G4-seq? y efecto de los terciles de
# senal y de G4 x senal alta sobre la actividad, tambien dentro de terciles
# de ATAC (ver el encabezado de la seccion).
#
# Parte D - poli(dA:dT) y elementos transponibles: fraccion de promotores
# con tracto >= 8 (y >= 8 proximal) que solapan un TE (TE_superclass de
# prom_df; SINE / LINE por separado), si el tracto cae dentro del TE
# (coordenadas de RepeatMasker) y el efecto de "≥ 8 proximal" / "≥ 8 solo
# rio arriba" vs. "< 6" excluyendo promotores con TE.
#
# Efecto = Hodges-Lehmann (IC 95%) y p de Wilcoxon de la actividad media,
# cada replica por separado; BH dentro de replica (familia = todas las
# comparaciones x estratos de una parte; en G4 ademas por metodo).
# Estratos (la referencia es la del mismo estrato): todos, con / sin isla
# CpG (CGI_frag) y terciles del contenido de G+C del fragmento completo
# (252 pb), calculados entre los promotores con actividad de cada replica.
#
# Geometria: fragmento de 252 pb de -236 a +15, TSS (= 0) en la posicion
# 237 de la secuencia orientada (TATA_TSS_POS en R/functions/tata_scan.R).
# Salidas (exploratorio, sin fig_dir ni mapping_figuras.csv) en
# figures/exploratorio_polyAT_g4/.
# Requiere: data/processed/activity_stats_highconf.tsv,
# data/processed/prom_df.tsv, data/processed/g4_cuttag_hek293t_v1.tsv
# (R/00_prom_features/build_g4_cuttag.R), data/processed/g4seq_lib_v1.tsv y
# g4seq_hits_lib_v1.tsv (R/00_prom_features/build_g4seq.R) y la secuencia de la library
# (../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa);
# para la Parte D, data/raw/library.bed y RepeatMasker (path_repeatmasker
# en R/00_prom_features/heavy_data_paths.R).
# Run from the TesisDoc repo root.

library(tidyverse)
library(ggpubr)
library(patchwork)
library(Biostrings)
library(rtracklayer) # import.bed (Parte D)
source("R/functions/tata_scan.R") # TATA_TSS_POS, TATA_FRAG_LEN, tata_hl_vs()

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename # Biostrings/S4Vectors las enmascaran
first <- dplyr::first

out_dir <- "figures/exploratorio_polyAT_g4"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
path_fa <- "../promoter_MPRA/Library_data/data/promoters_wo_dupl_SpikeIns.fa"

TSS_POS <- TATA_TSS_POS
MIN_N <- 10
PROXIMAL <- c(-100, -1)
RIO_ARRIBA <- c(-236, -101)
G4_REGEX <- c(
  Codificante = "G{3,}[ACGT]{1,7}?G{3,}[ACGT]{1,7}?G{3,}[ACGT]{1,7}?G{3,}",
  Molde = "C{3,}[ACGT]{1,7}?C{3,}[ACGT]{1,7}?C{3,}[ACGT]{1,7}?C{3,}"
)

AT_CATS <- c("< 6", "6–7", "8–9", "≥ 10")
AT_COMPS <- c(cat_6_7 = "6–7", cat_8_9 = "8–9", cat_ge10 = "≥ 10", prox8 = "≥ 8 proximal (−100 a −1)", arriba8 = "≥ 8 solo río arriba (−236 a −101)")
G4_CLASES <- c("Sin G4", "Solo codificante", "Solo molde", "Ambas")
G4_COMPS <- c(g4_cualquiera = "Cualquier hebra", g4_solo_cod = "Solo codificante", g4_solo_molde = "Solo molde", g4_ambas = "Ambas")
TERCILES <- c("G+C tercil bajo", "G+C tercil medio", "G+C tercil alto")
ESTRATOS <- c("Todos", "Con isla CpG", "Sin isla CpG", TERCILES)
GRUPOS <- c("Todos", "Con isla CpG", "Sin isla CpG")

col_alto <- "#D6741F"
col_bajo <- "#7FB800"
col_hebra <- c(Codificante = "#D6741F", Molde = "#216869")
col_grupo <- c("Todos" = "grey30", "Con isla CpG" = "#216869", "Sin isla CpG" = "#14AFB2")
col_efecto <- scale_color_manual(values = c("Mayor actividad" = col_alto, "Menor actividad" = col_bajo, "n.s." = "grey60"), drop = FALSE)
th <- theme_pubclean(base_size = 13) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")
minus <- function(x) sub("^-", "−", format(x, trim = TRUE))
x_pos <- scale_x_continuous(breaks = c(-200, -150, -100, -50, 0), labels = minus, expand = expansion(0.01))
save_tsv <- function(x, name) write_tsv(x, file.path(out_dir, name))
efecto_lab <- function(estimate, p_adj) {
  case_when(
    p_adj < 0.05 & estimate > 0 ~ "Mayor actividad",
    p_adj < 0.05 & estimate < 0 ~ "Menor actividad",
    TRUE ~ "n.s."
  ) %>% factor(levels = c("Mayor actividad", "Menor actividad", "n.s."))
}

# --- Datos -------------------------------------------------------------------

prom_df <- read_tsv("data/processed/prom_df.tsv", show_col_types = FALSE) %>% filter(type == "promoter")
fa <- readDNAStringSet(path_fa)
fa <- fa[names(fa) %in% prom_df$seq_id]
stopifnot(all(width(fa) == TATA_FRAG_LEN), all(prom_df$seq_id %in% names(fa)))
seqs <- as.character(fa)

# todas las coincidencias de un patron, con su centro respecto del TSS
buscar <- function(patron, seqs) {
  m <- gregexpr(patron, seqs, perl = TRUE)
  tibble(
    seq_id = rep(names(seqs), lengths(m)),
    start = unlist(m),
    len = unlist(map(m, attr, "match.length"))
  ) %>%
    filter(start > 0) %>%
    mutate(pos = start + (len - 1) / 2 - TSS_POS)
}

gc_frag <- letterFrequency(fa, "GC", as.prob = TRUE)[, 1] %>% set_names(names(fa))

# --- Parte A: tractos poli(dA:dT) ---------------------------------------------

tractos <- buscar("A+|T+", seqs) %>%
  mutate(base = substr(seqs[seq_id], start, start)) # A o T

at_prom <- tractos %>%
  group_by(seq_id) %>%
  arrange(desc(len), abs(pos), .by_group = TRUE) %>%
  summarise(
    at_max_len = first(len), at_max_pos = first(pos), at_max_base = first(base),
    n_tractos_ge6 = sum(len >= 6), n_tractos_ge8 = sum(len >= 8),
    ge8_proximal = any(len >= 8 & pos >= PROXIMAL[1] & pos <= PROXIMAL[2]),
    ge8_rio_arriba = any(len >= 8 & pos < PROXIMAL[1]),
    .groups = "drop"
  ) %>%
  right_join(tibble(seq_id = names(seqs)), by = "seq_id") %>%
  mutate(
    at_max_len = replace_na(at_max_len, 0L),
    across(c(n_tractos_ge6, n_tractos_ge8), ~ replace_na(.x, 0L)),
    across(c(ge8_proximal, ge8_rio_arriba), ~ replace_na(.x, FALSE)),
    at_cat = cut(at_max_len, c(-Inf, 5, 7, 9, Inf), labels = AT_CATS),
    at_pos_cat = case_when(
      ge8_proximal ~ "proximal",
      ge8_rio_arriba ~ "rio_arriba",
      n_tractos_ge8 > 0 ~ "solo_tss_o_despues",
      TRUE ~ NA_character_
    ),
    cat_lt6 = at_cat == "< 6", cat_6_7 = at_cat == "6–7", cat_8_9 = at_cat == "8–9", cat_ge10 = at_cat == "≥ 10",
    prox8 = at_pos_cat %in% "proximal", arriba8 = at_pos_cat %in% "rio_arriba"
  )

# --- Parte B: G4 --------------------------------------------------------------

g4_hits <- imap_dfr(G4_REGEX, ~ buscar(.x, seqs) %>% mutate(hebra = .y)) %>% mutate(metodo = "Regex canónico")

usa_pqsfinder <- requireNamespace("pqsfinder", quietly = TRUE)
if (usa_pqsfinder) {
  pqs <- tryCatch(
    imap_dfr(seqs, function(s, id) {
      v <- pqsfinder::pqsfinder(DNAString(s))
      if (length(v) == 0) return(NULL)
      tibble(seq_id = id, start = start(v), len = width(v), strand = as.character(pqsfinder::strand(v)), score = pqsfinder::score(v))
    }),
    error = function(e) {
      message("pqsfinder fallo: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(pqs) && nrow(pqs) > 0) {
    g4_hits <- bind_rows(
      g4_hits,
      pqs %>% mutate(pos = start + (len - 1) / 2 - TSS_POS, hebra = if_else(strand == "+", "Codificante", "Molde"), metodo = "pqsfinder") %>% select(-strand)
    )
  } else {
    usa_pqsfinder <- FALSE
  }
}
if (!usa_pqsfinder) message("pqsfinder no esta disponible: no se reporta.")

# G4 observados in vitro (G4-seq, Marsico et al. 2019), K+ y K+ + PDS, por
# hebra respecto del promotor (R/00_prom_features/build_g4seq.R). Los
# fragmentos que no pasan a hg19 quedan fuera de estos dos metodos.
g4seq_ok <- read_tsv("data/processed/g4seq_lib_v1.tsv", show_col_types = FALSE) %>%
  filter(liftover_ok) %>%
  distinct(seq_id) %>%
  pull()
g4_hits <- bind_rows(
  g4_hits,
  read_tsv("data/processed/g4seq_hits_lib_v1.tsv", show_col_types = FALSE) %>%
    filter(seq_id %in% names(seqs)) %>%
    transmute(seq_id, len = largo, pos, hebra, metodo = paste("G4-seq", condicion))
)
METODOS <- unique(g4_hits$metodo)
METODO_REGEX <- "Regex canónico"

# un registro por promotor x metodo; posicion = la coincidencia mas cercana al TSS
g4_prom <- expand_grid(seq_id = names(seqs), metodo = METODOS) %>%
  left_join(
    g4_hits %>%
      group_by(seq_id, metodo) %>%
      summarise(
        n_g4_cod = sum(hebra == "Codificante"), n_g4_molde = sum(hebra == "Molde"),
        pos_g4_cod = pos[hebra == "Codificante"][which.min(abs(pos[hebra == "Codificante"]))][1],
        pos_g4_molde = pos[hebra == "Molde"][which.min(abs(pos[hebra == "Molde"]))][1],
        .groups = "drop"
      ),
    by = c("seq_id", "metodo")
  ) %>%
  mutate(
    across(c(n_g4_cod, n_g4_molde), ~ replace_na(.x, 0L)),
    g4_clase = case_when(
      n_g4_cod > 0 & n_g4_molde > 0 ~ "Ambas",
      n_g4_cod > 0 ~ "Solo codificante",
      n_g4_molde > 0 ~ "Solo molde",
      TRUE ~ "Sin G4"
    ) %>% factor(G4_CLASES),
    g4_ninguno = g4_clase == "Sin G4", g4_cualquiera = !g4_ninguno,
    g4_solo_cod = g4_clase == "Solo codificante", g4_solo_molde = g4_clase == "Solo molde", g4_ambas = g4_clase == "Ambas",
    metodo = factor(metodo, METODOS)
  ) %>%
  filter(!(str_starts(metodo, "G4-seq") & !seq_id %in% g4seq_ok))

# concordancia regex vs. G4-seq por hebra (promotores con liftover)
g4_concordancia <- g4_prom %>%
  filter(seq_id %in% g4seq_ok) %>%
  select(seq_id, metodo, n_g4_cod, n_g4_molde) %>%
  pivot_longer(c(n_g4_cod, n_g4_molde), names_to = "hebra", values_to = "n") %>%
  mutate(hebra = if_else(hebra == "n_g4_cod", "Codificante", "Molde"), tiene = n > 0) %>%
  select(-n) %>%
  pivot_wider(names_from = metodo, values_from = tiene) %>%
  pivot_longer(starts_with("G4-seq"), names_to = "g4seq", values_to = "obs") %>%
  group_by(g4seq, hebra) %>%
  summarise(
    n = n(), n_regex = sum(.data[[METODO_REGEX]]), n_g4seq = sum(obs), n_ambos = sum(.data[[METODO_REGEX]] & obs),
    pct_regex_observado = 100 * n_ambos / n_regex, pct_g4seq_con_regex = 100 * n_ambos / n_g4seq,
    jaccard = n_ambos / (n_regex + n_g4seq - n_ambos), .groups = "drop"
  )
print(g4_concordancia)

# tabla por promotor (insumo para otros analisis)
save_tsv(
  prom_df %>%
    select(seq_id, gene_sym, CGI_frag) %>%
    mutate(gc_frag = gc_frag[seq_id]) %>%
    left_join(at_prom %>% select(seq_id, at_max_len, at_max_pos, at_max_base, n_tractos_ge6, n_tractos_ge8, at_cat, at_pos_cat), by = "seq_id") %>%
    left_join(g4_prom %>% select(seq_id, metodo, g4_clase, n_g4_cod, n_g4_molde, pos_g4_cod, pos_g4_molde) %>%
      pivot_wider(names_from = metodo, values_from = -c(seq_id, metodo), names_glue = "{.value}__{metodo}"), by = "seq_id"),
  "polyAT_g4_por_promotor.tsv"
)

# --- Actividad y estratos ----------------------------------------------------

data <- read_tsv("data/processed/activity_stats_highconf.tsv", show_col_types = FALSE) %>%
  inner_join(prom_df %>% select(seq_id, name, CGI_frag), by = c("seq_id", "name")) %>%
  mutate(gc_frag = gc_frag[seq_id]) %>%
  group_by(rep) %>%
  mutate(tercil_gc = factor(TERCILES[ntile(gc_frag, 3)], TERCILES)) %>%
  ungroup()

terciles_rango <- data %>%
  group_by(rep, tercil_gc) %>%
  summarise(n = n(), gc_min = min(gc_frag), gc_max = max(gc_frag), gc_mediana = median(gc_frag),
            pct_con_isla = 100 * mean(CGI_frag, na.rm = TRUE), .groups = "drop")
print(terciles_rango)

con_estratos <- function(d) {
  bind_rows(
    d %>% mutate(estrato = "Todos"),
    d %>% filter(!is.na(CGI_frag)) %>% mutate(estrato = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")),
    d %>% mutate(estrato = as.character(tercil_gc))
  ) %>% mutate(estrato = factor(estrato, ESTRATOS))
}

# HL de cada indicador logico (`comps`, nombre = columna) vs. `ref`, por
# replica x estrato (x `extra`); BH dentro de replica (x `extra`).
efectos <- function(df, comps, ref, extra = character()) {
  df %>%
    group_by(across(all_of(c("rep", extra, "estrato")))) %>%
    group_modify(function(d, k) {
      r <- d$mean[d[[ref]]]
      imap_dfr(comps, ~ tata_hl_vs(d$mean[d[[.y]]], r, MIN_N) %>%
        mutate(comparacion = .x, n = sum(d[[.y]]), n_ref = length(r), n_estrato = nrow(d)))
    }) %>%
    group_by(across(all_of(c("rep", extra)))) %>%
    mutate(p_adj = p.adjust(p, "BH")) %>%
    ungroup() %>%
    mutate(efecto = efecto_lab(estimate, p_adj))
}

# grafico de puntos con IC (formato de las clases de TATA)
forest <- function(res, x_lab, facets = estrato ~ rep) {
  ggplot(res, aes(estimate, comparacion, color = efecto)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_linerange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.9) +
    geom_point(size = 2.8) +
    geom_text(aes(label = paste0("n=", n)), vjust = -1, size = 3, color = "grey30") +
    facet_grid(facets) +
    col_efecto +
    th +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(x = x_lab, y = NULL, color = NULL, caption = "BH dentro de réplica; la referencia es la del mismo estrato.")
}

# universos para las prevalencias: library completa y promotores con actividad de cada replica
universos <- bind_rows(
  prom_df %>% select(seq_id, CGI_frag) %>% mutate(universo = "Library (promotores FP)"),
  data %>% select(seq_id, CGI_frag, rep) %>% mutate(universo = paste0("Con actividad (", rep, ")")) %>% select(-rep)
)
por_grupo <- function(u) {
  bind_rows(
    u %>% mutate(grupo = "Todos"),
    u %>% filter(!is.na(CGI_frag)) %>% mutate(grupo = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG"))
  ) %>% mutate(grupo = factor(grupo, GRUPOS))
}

# =============================================================================
# PARTE A - resultados
# =============================================================================

at_prev <- universos %>%
  left_join(at_prom, by = "seq_id") %>%
  por_grupo() %>%
  group_by(universo, grupo) %>%
  mutate(n_grupo = n(), media_tractos_ge6 = mean(n_tractos_ge6), media_tractos_ge8 = mean(n_tractos_ge8)) %>%
  group_by(universo, grupo, n_grupo, media_tractos_ge6, media_tractos_ge8, categoria = at_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = 100 * n / n_grupo)

at_pos_n <- universos %>%
  left_join(at_prom, by = "seq_id") %>%
  por_grupo() %>%
  filter(n_tractos_ge8 > 0) %>%
  count(universo, grupo, categoria = paste("≥ 8:", at_pos_cat), name = "n") %>%
  group_by(universo, grupo) %>%
  mutate(n_grupo = sum(n), pct = 100 * n / n_grupo) %>%
  ungroup()

data_at <- data %>% left_join(at_prom, by = "seq_id") %>% con_estratos()
at_ef <- efectos(data_at, AT_COMPS, "cat_lt6")

save_tsv(
  bind_rows(
    at_prev %>% mutate(analisis = "prevalencia_categorias"),
    at_pos_n %>% mutate(analisis = "prevalencia_posicion_ge8"),
    at_ef %>% mutate(analisis = "efecto_HL_vs_menor6"),
    terciles_rango %>% mutate(analisis = "terciles_gc")
  ) %>% relocate(analisis),
  "polyAT.tsv"
)

p_at <- at_ef %>%
  mutate(comparacion = factor(comparacion, rev(AT_COMPS))) %>%
  forest("Diferencia de actividad media vs. tracto más largo < 6\n(Hodges-Lehmann, IC 95%)")
ggsave(file.path(out_dir, "polyAT_efecto.jpg"), p_at, width = 12, height = 15, units = "in")

# posicion de todos los tractos >= 8 en la library
at_hist <- tractos %>%
  filter(len >= 8) %>%
  inner_join(prom_df %>% select(seq_id, CGI_frag), by = "seq_id") %>%
  mutate(CGI_frag = factor(if_else(CGI_frag, "Con isla CpG", "Sin isla CpG"), GRUPOS[-1])) %>%
  filter(!is.na(CGI_frag)) %>%
  ggplot(aes(pos, fill = base)) +
  annotate("rect", xmin = PROXIMAL[1] - 0.5, xmax = PROXIMAL[2] + 0.5, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_histogram(binwidth = 5, boundary = 0, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c(A = "#7FB800", T = "#B5345A"), labels = c(A = "A_n", T = "T_n"), name = "Tracto (hebra del promotor)") +
  x_pos +
  coord_cartesian(xlim = c(-237, 16)) +
  facet_grid(CGI_frag ~ ., scales = "free_y") +
  th +
  labs(x = "Centro del tracto respecto del TSS (TSS = 0; sombreado: −100 a −1)", y = "Tractos ≥ 8 pb",
       caption = "Todos los tractos ≥ 8 de los promotores FP de la library (puede haber más de uno por promotor).")
ggsave(file.path(out_dir, "polyAT_posicion_hist.jpg"), at_hist, width = 10, height = 6.5, units = "in")

# =============================================================================
# PARTE B - resultados
# =============================================================================

g4_u <- universos %>%
  left_join(g4_prom, by = "seq_id", relationship = "many-to-many")

g4_prev <- g4_u %>%
  por_grupo() %>%
  group_by(metodo, universo, grupo) %>%
  mutate(n_grupo = n()) %>%
  group_by(metodo, universo, grupo, n_grupo, clase = g4_clase) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(nesting(metodo, universo, grupo, n_grupo), clase, fill = list(n = 0L)) %>%
  mutate(pct = 100 * n / n_grupo)

g4_fisher <- g4_u %>%
  filter(!is.na(CGI_frag)) %>%
  group_by(metodo, universo) %>%
  group_modify(function(d, k) {
    imap_dfr(G4_COMPS, function(lab, col) {
      ft <- fisher.test(table(factor(d[[col]], c(TRUE, FALSE)), factor(d$CGI_frag, c(TRUE, FALSE))))
      tibble(
        clase = lab, n_con_cgi = sum(d$CGI_frag), n_g4_con_cgi = sum(d[[col]] & d$CGI_frag),
        n_sin_cgi = sum(!d$CGI_frag), n_g4_sin_cgi = sum(d[[col]] & !d$CGI_frag),
        odds_ratio = unname(ft$estimate), or_conf.low = ft$conf.int[1], or_conf.high = ft$conf.int[2], p_fisher = ft$p.value
      ) %>% mutate(pct_con_cgi = 100 * n_g4_con_cgi / n_con_cgi, pct_sin_cgi = 100 * n_g4_sin_cgi / n_sin_cgi)
    })
  }) %>%
  ungroup()
print(g4_fisher %>% filter(universo == "Library (promotores FP)") %>% select(metodo, clase, pct_con_cgi, pct_sin_cgi, odds_ratio, p_fisher))

data_g4 <- data %>%
  left_join(g4_prom, by = "seq_id", relationship = "many-to-many") %>%
  con_estratos()
g4_ef <- efectos(data_g4, G4_COMPS, "g4_ninguno", extra = "metodo")

save_tsv(
  bind_rows(
    g4_prev %>% mutate(clase = as.character(clase), metodo = as.character(metodo), analisis = "prevalencia_clases"),
    g4_fisher %>% mutate(metodo = as.character(metodo), analisis = "asociacion_isla_fisher"),
    g4_concordancia %>% mutate(analisis = "concordancia_regex_g4seq"),
    g4_ef %>% mutate(metodo = as.character(metodo), analisis = "efecto_HL_vs_sin_G4"),
    terciles_rango %>% mutate(analisis = "terciles_gc")
  ) %>% relocate(analisis, any_of("metodo")),
  "g4.tsv"
)

p_g4 <- g4_ef %>%
  mutate(comparacion = factor(comparacion, rev(G4_COMPS))) %>%
  forest("Diferencia de actividad media vs. sin G4\n(Hodges-Lehmann, IC 95%)",
         facets = if (length(METODOS) > 1) estrato ~ rep + metodo else estrato ~ rep)
ggsave(file.path(out_dir, "g4_efecto.jpg"), p_g4, width = 6 + 5 * length(METODOS), height = 13, units = "in")

g4_hist <- g4_hits %>%
  inner_join(prom_df %>% select(seq_id, CGI_frag), by = "seq_id") %>%
  filter(!is.na(CGI_frag)) %>%
  mutate(CGI_frag = factor(if_else(CGI_frag, "Con isla CpG", "Sin isla CpG"), GRUPOS[-1])) %>%
  ggplot(aes(pos, fill = hebra)) +
  geom_histogram(binwidth = 5, boundary = 0, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = col_hebra, guide = "none") +
  x_pos +
  coord_cartesian(xlim = c(-237, 16)) +
  facet_grid(CGI_frag ~ metodo + hebra, scales = "free_y") +
  th +
  labs(x = "Centro del G4 respecto del TSS (TSS = 0)", y = "G4",
       caption = "Todas las coincidencias de los promotores FP de la library. Codificante: G4 en la hebra del promotor; molde: G4 en la hebra molde.")
ggsave(file.path(out_dir, "g4_posicion_hist.jpg"), g4_hist, width = 5 + 4 * length(METODOS), height = 6.5, units = "in")

# =============================================================================
# PARTE C - G4 plegado en HEK293T (G4-CUT&Tag, GSE178668)
# =============================================================================
# Senal por fragmento = promedio de log2((BG4 + 0.1) / (IgG + 0.1)) de las
# dos replicas de CUT&Tag (R/00_prom_features/build_g4_cuttag.R). Es senal
# del locus endogeno, sin hebra y sin picos (GEO solo tiene bigWig); la
# distribucion es unimodal, asi que no hay un corte natural "plegado / no
# plegado": se usan terciles de la senal calculados dentro de cada replica
# x estrato (alto / medio vs. bajo), para que no reflejen la isla CpG ni
# el G+C del estrato. Estratos extra: terciles de ATAC-seq de HEK293T
# (GSM5395701) dentro de cada replica, para comparar senal alta vs. baja
# entre regiones igual de accesibles. Combinacion: G4 (cualquier hebra; por
# regex o por G4-seq) x senal en el tercil alto, vs. "ninguno".

ct <- read_tsv("data/processed/g4_cuttag_hek293t_v1.tsv", show_col_types = FALSE) %>%
  transmute(seq_id, ct_rep1 = log2enr_rep1, ct_rep2 = log2enr_rep2, ct_enr = (log2enr_rep1 + log2enr_rep2) / 2,
            atac_log2 = log2(atac + 0.1))
CT_COMPS <- c(
  ct_alto = "Señal CUT&Tag alta (T3) vs. baja (T1)", ct_medio = "Señal CUT&Tag media (T2) vs. baja (T1)",
  comb_seq = "Solo G4 (secuencia)", comb_ct = "Solo señal alta", comb_ambos = "G4 (secuencia) + señal alta"
)
TERCILES_ATAC <- c("ATAC tercil bajo", "ATAC tercil medio", "ATAC tercil alto")
ESTRATOS_CT <- c(ESTRATOS, TERCILES_ATAC)

# validacion: ¿los promotores con G4 (regex u observado por G4-seq) tienen
# mas senal de CUT&Tag? En todos, con / sin isla CpG y por tercil de ATAC
# (library completa).
ct_lib <- prom_df %>%
  select(seq_id, CGI_frag) %>%
  inner_join(ct, by = "seq_id") %>%
  mutate(tercil_atac = factor(TERCILES_ATAC[ntile(atac_log2, 3)], TERCILES_ATAC)) %>%
  inner_join(g4_prom, by = "seq_id", relationship = "many-to-many")
ct_vs_seq <- bind_rows(
  ct_lib %>% por_grupo(),
  ct_lib %>% mutate(grupo = tercil_atac)
) %>%
  mutate(grupo = factor(grupo, c(GRUPOS, TERCILES_ATAC))) %>%
  group_by(metodo, grupo) %>%
  group_modify(function(d, k) {
    imap_dfr(G4_COMPS, ~ tata_hl_vs(d$ct_enr[d[[.y]]], d$ct_enr[d$g4_ninguno], MIN_N) %>%
      mutate(clase = .x, n = sum(d[[.y]]), n_ref = sum(d$g4_ninguno),
             mediana = median(d$ct_enr[d[[.y]]]), mediana_ref = median(d$ct_enr[d$g4_ninguno])))
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, "BH"))
print(ct_vs_seq %>% filter(clase == "Cualquier hebra") %>% select(metodo, grupo, n, mediana, mediana_ref, estimate, p_adj), n = Inf)

ct_cor <- tibble(
  medida = c("Spearman CUT&Tag rep1 vs rep2", "Spearman CUT&Tag vs G+C del fragmento", "Spearman CUT&Tag vs ATAC"),
  rho = c(
    cor(ct$ct_rep1, ct$ct_rep2, method = "spearman"),
    cor(ct$ct_enr, gc_frag[ct$seq_id], method = "spearman", use = "complete.obs"),
    cor(ct$ct_enr, ct$atac_log2, method = "spearman")
  )
)

# efectos sobre la actividad: estratos de siempre + terciles de ATAC
# (dentro de cada replica), con los terciles de senal de CUT&Tag calculados
# dentro de replica x estrato (x metodo, que no cambia la senal)
data_ct <- data %>%
  inner_join(ct, by = "seq_id") %>%
  group_by(rep) %>%
  mutate(tercil_atac = TERCILES_ATAC[ntile(atac_log2, 3)]) %>%
  ungroup() %>%
  (function(d) bind_rows(con_estratos(d) %>% mutate(estrato = as.character(estrato)), d %>% mutate(estrato = tercil_atac))) %>%
  mutate(estrato = factor(estrato, ESTRATOS_CT)) %>%
  inner_join(g4_prom %>% select(seq_id, metodo, g4_cualquiera), by = "seq_id", relationship = "many-to-many") %>%
  group_by(rep, estrato, metodo) %>%
  mutate(t_ct = ntile(ct_enr, 3)) %>%
  ungroup() %>%
  mutate(
    ct_bajo = t_ct == 1, ct_medio = t_ct == 2, ct_alto = t_ct == 3,
    comb_ninguno = !g4_cualquiera & !ct_alto, comb_seq = g4_cualquiera & !ct_alto,
    comb_ct = !g4_cualquiera & ct_alto, comb_ambos = g4_cualquiera & ct_alto
  )

ct_ef <- bind_rows(
  efectos(data_ct %>% filter(metodo == METODO_REGEX), CT_COMPS[c("ct_alto", "ct_medio")], "ct_bajo") %>%
    mutate(metodo = "(solo señal)", analisis = "efecto_HL_tercil_senal"),
  efectos(data_ct, CT_COMPS[c("comb_seq", "comb_ct", "comb_ambos")], "comb_ninguno", extra = "metodo") %>%
    mutate(metodo = as.character(metodo), analisis = "efecto_HL_secuencia_x_senal")
)
# rango de la senal en cada tercil
ct_terciles <- data_ct %>%
  group_by(rep, estrato, metodo, t_ct) %>%
  summarise(n = n(), ct_min = min(ct_enr), ct_max = max(ct_enr), ct_mediana = median(ct_enr),
            pct_g4_secuencia = 100 * mean(g4_cualquiera), .groups = "drop") %>%
  mutate(metodo = as.character(metodo))

save_tsv(
  bind_rows(
    ct_cor %>% mutate(analisis = "correlaciones"),
    ct_vs_seq %>% mutate(metodo = as.character(metodo), grupo = as.character(grupo), analisis = "senal_CUTTag_por_clase_G4"),
    ct_terciles %>% mutate(analisis = "terciles_senal_CUTTag"),
    ct_ef
  ) %>% relocate(analisis, metodo),
  "g4_cuttag.tsv"
)

x_lab_ct <- "Diferencia de actividad media (Hodges-Lehmann, IC 95%)"
p_ct_terc <- ct_ef %>%
  filter(analisis == "efecto_HL_tercil_senal") %>%
  mutate(comparacion = factor(comparacion, rev(CT_COMPS))) %>%
  forest(paste0(x_lab_ct, ": terciles de señal G4-CUT&Tag vs. tercil bajo"))
ggsave(file.path(out_dir, "g4_cuttag_efecto_terciles.jpg"), p_ct_terc, width = 12, height = 16, units = "in")
p_ct_comb <- ct_ef %>%
  filter(analisis == "efecto_HL_secuencia_x_senal") %>%
  mutate(comparacion = factor(comparacion, rev(CT_COMPS))) %>%
  forest(paste0(x_lab_ct, ": G4 (cualquier hebra) × señal G4-CUT&Tag alta vs. ninguno"), facets = estrato ~ rep + metodo)
ggsave(file.path(out_dir, "g4_cuttag_efecto_secuencia.jpg"), p_ct_comb, width = 22, height = 18, units = "in")

p_ct_seq <- bind_rows(
  ct_lib %>% filter(!is.na(CGI_frag)) %>% mutate(panel = if_else(CGI_frag, "Con isla CpG", "Sin isla CpG")),
  ct_lib %>% mutate(panel = as.character(tercil_atac))
) %>%
  mutate(panel = factor(panel, c(GRUPOS[-1], TERCILES_ATAC))) %>%
  ggplot(aes(g4_clase, ct_enr, fill = g4_clase)) +
  geom_violin(color = NA, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("Sin G4" = "grey60", "Solo codificante" = "#D6741F", "Solo molde" = "#216869", "Ambas" = "#7A5195"), guide = "none") +
  facet_grid(metodo ~ panel) +
  th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "G4 por secuencia (regex) u observado in vitro (G4-seq)", y = "Señal G4-CUT&Tag en HEK293T\nlog2(BG4 / IgG), media de 2 réplicas",
       caption = "Promotores FP de la library; señal media en los 252 pb del fragmento (GSE178668). Terciles de ATAC-seq de HEK293T (GSM5395701) sobre la library.")
ggsave(file.path(out_dir, "g4_cuttag_vs_secuencia.jpg"), p_ct_seq, width = 16, height = 10, units = "in")

# =============================================================================
# PARTE D - poli(dA:dT) y elementos transponibles (TE)
# =============================================================================
# Los tractos >= 8 estan enriquecidos en promotores con TE (sobre todo Alu,
# cuya cola poli(A) / linker rico en A es un tracto largo). ¿El efecto de
# "≥ 8 proximal" vs. "< 6" se debe a eso? TE por promotor = TE_superclass
# de prom_df (superposicion > 10 pb, sin Low_complexity / Simple_repeat;
# guarda la primera superclase, no la dominante). Para saber si el tracto
# cae dentro del TE se leen las coordenadas de RepeatMasker (misma regla que
# build_prom_features.R) y se pasan al fragmento orientado (1..252).
# Conjuntos: todos; sin ningun TE; sin SINE ni LINE; excluyendo solo los
# promotores con algun tracto >= 8 proximal entero dentro de un TE. En cada
# conjunto la referencia "< 6" es la del mismo conjunto x estrato. BH
# dentro de replica x conjunto.

source("R/00_prom_features/heavy_data_paths.R")
TE_CONJ <- c("Todos los promotores", "Sin TE", "Sin SINE ni LINE", "Sin tracto proximal dentro de TE")

lib <- import.bed("data/raw/library.bed")
lib <- lib[lib$name %in% names(seqs)]
rm_te <- data.table::fread(path_repeatmasker, sep = " ", fill = TRUE, select = c("seqname", "start", "end", "matching_repeat", "repeat_class"))
rm_te <- GRanges(as.data.frame(rm_te[!repeat_class %in% c("Low_complexity", "Simple_repeat") & seqname %in% seqlevels(lib)]))
ov <- findOverlaps(lib, rm_te)
q <- lib[queryHits(ov)]
s <- pmax(start(q), start(rm_te)[subjectHits(ov)])
e <- pmin(end(q), end(rm_te)[subjectHits(ov)])
mas <- as.character(strand(q)) == "+"
te_frag <- tibble(
  seq_id = q$name,
  te_start = if_else(mas, s - start(q) + 1L, end(q) - e + 1L),
  te_end = if_else(mas, e - start(q) + 1L, end(q) - s + 1L),
  te_nombre = rm_te$matching_repeat[subjectHits(ov)], te_clase = rm_te$repeat_class[subjectHits(ov)]
) %>% filter(te_end - te_start + 1 > 10)
rm(rm_te, ov, q, s, e, mas)

# cada tracto >= 8: entero dentro de un TE, solapamiento parcial o fuera
tr_te <- tractos %>%
  filter(len >= 8) %>%
  mutate(end = start + len - 1, proximal = pos >= PROXIMAL[1] & pos <= PROXIMAL[2]) %>%
  left_join(te_frag, by = "seq_id", relationship = "many-to-many") %>%
  mutate(dentro = !is.na(te_start) & start >= te_start & end <= te_end,
         solapa = !is.na(te_start) & start <= te_end & end >= te_start) %>%
  group_by(seq_id, start, len, pos, proximal) %>%
  summarise(dentro_TE = any(dentro), solapa_parcial_TE = any(solapa) & !dentro_TE,
            te_clase_dentro = paste(sort(unique(te_clase[dentro])), collapse = "+"),
            te_nombre_dentro = paste(sort(unique(te_nombre[dentro])), collapse = "+"), .groups = "drop")

te_info <- prom_df %>%
  transmute(seq_id, TE_superclass = as.character(TE_superclass), con_TE = !is.na(TE_superclass)) %>%
  left_join(tr_te %>% group_by(seq_id) %>% summarise(prox_dentro_TE = any(proximal & dentro_TE), algun_ge8_dentro_TE = any(dentro_TE)), by = "seq_id") %>%
  mutate(across(c(prox_dentro_TE, algun_ge8_dentro_TE), ~ replace_na(.x, FALSE)))

# prevalencia de TE por categoria de tracto (library y con actividad)
te_prev <- universos %>%
  left_join(at_prom %>% select(seq_id, at_cat, prox8, n_tractos_ge8), by = "seq_id") %>%
  left_join(te_info, by = "seq_id") %>%
  por_grupo() %>%
  (function(d) bind_rows(
    d %>% filter(n_tractos_ge8 > 0) %>% mutate(categoria = "≥ 8 (cualquier posición)"),
    d %>% filter(prox8) %>% mutate(categoria = AT_COMPS[["prox8"]]),
    d %>% filter(at_cat == "< 6") %>% mutate(categoria = "< 6")
  )) %>%
  group_by(universo, grupo, categoria) %>%
  summarise(
    n = n(), pct_TE = 100 * mean(con_TE), pct_SINE = 100 * mean(TE_superclass %in% "SINE"),
    pct_LINE = 100 * mean(TE_superclass %in% "LINE"), pct_LTR = 100 * mean(TE_superclass %in% "LTR"),
    pct_DNA = 100 * mean(TE_superclass %in% "DNA"), pct_algun_ge8_dentro_TE = 100 * mean(algun_ge8_dentro_TE),
    pct_tracto_prox_dentro_TE = 100 * mean(prox_dentro_TE),
    pct_prox_dentro_entre_con_TE = 100 * sum(prox_dentro_TE & con_TE) / sum(con_TE), .groups = "drop"
  )
print(te_prev %>% filter(universo == "Library (promotores FP)"), width = Inf)

# tractos individuales (library): dentro / parcial / fuera, y familia del TE
tr_te_res <- bind_rows(
  tr_te %>%
    group_by(region = if_else(proximal, "proximal", "fuera de −100 a −1")) %>%
    summarise(n = n(), pct_dentro_TE = 100 * mean(dentro_TE), pct_solapa_parcial = 100 * mean(solapa_parcial_TE), .groups = "drop") %>%
    mutate(tabla = "tractos_dentro_TE"),
  tr_te %>%
    filter(dentro_TE) %>%
    count(region = if_else(proximal, "proximal", "fuera de −100 a −1"), te_clase = te_clase_dentro, sort = TRUE) %>%
    mutate(tabla = "familia_TE_tractos_dentro")
)
print(tr_te_res %>% filter(tabla == "tractos_dentro_TE"))

d_base <- data %>% left_join(at_prom, by = "seq_id") %>% left_join(te_info, by = "seq_id")
data_te <- bind_rows(
  d_base %>% mutate(conjunto = TE_CONJ[1]),
  d_base %>% filter(!con_TE) %>% mutate(conjunto = TE_CONJ[2]),
  d_base %>% filter(!TE_superclass %in% c("SINE", "LINE")) %>% mutate(conjunto = TE_CONJ[3]),
  d_base %>% filter(!prox_dentro_TE) %>% mutate(conjunto = TE_CONJ[4])
) %>%
  con_estratos() %>%
  filter(estrato %in% GRUPOS) %>%
  mutate(estrato = droplevels(estrato))
te_ef <- efectos(data_te, AT_COMPS[c("prox8", "arriba8")], "cat_lt6", extra = "conjunto")
print(te_ef %>% filter(comparacion == AT_COMPS[["prox8"]]) %>% select(rep, conjunto, estrato, n, n_ref, estimate, conf.low, conf.high, p_adj), n = Inf)

save_tsv(
  bind_rows(
    te_prev %>% mutate(analisis = "prevalencia_TE"),
    tr_te_res %>% rename(analisis = tabla),
    te_ef %>% mutate(analisis = "efecto_HL_vs_menor6_por_conjunto")
  ) %>% relocate(analisis),
  "polyAT_TE.tsv"
)

p_te <- te_ef %>%
  mutate(tracto = factor(comparacion, AT_COMPS[c("prox8", "arriba8")]), comparacion = factor(conjunto, rev(TE_CONJ))) %>%
  forest("Diferencia de actividad media vs. tracto más largo < 6\n(Hodges-Lehmann, IC 95%), excluyendo elementos transponibles",
         facets = estrato ~ rep + tracto) +
  labs(caption = "BH dentro de réplica × conjunto; la referencia es la del mismo conjunto y estrato.")
ggsave(file.path(out_dir, "polyAT_TE_efecto.jpg"), p_te, width = 16, height = 9, units = "in")

message("Figuras y tablas guardadas en ", out_dir)
