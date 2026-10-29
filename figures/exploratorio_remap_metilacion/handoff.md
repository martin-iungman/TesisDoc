# Eje de metilación en la unión de factores (ReMap) y el ruido, controlado por cobertura

Scripts: `R/50_exploratorio_remap_metilacion/remap_metilacion.R` (p1–p6, ~5 min) y
`R/50_exploratorio_remap_metilacion/remap_pausa_tfiia.R` (p7–p8, ~4 min; lee las salidas
p1 y p4 del primero). Salidas: esta carpeta, con prefijo por punto del pedido. Todo en JPG y PDF.
Las secciones p7 (pausa / P-TEFb) y p8 (GTF2A2 y TATA-box, SNAPC1) están al final.
Todos los resultados se dan como Rep 1 / Rep 2.

## Resumen

- **El lado «CpG no metilado» se sostiene; el lado «metilación» no.**
- Entre factores de cobertura comparable, MAZ, SP2, ZBTB14, E2F1 y E2F4 están entre los que más reducen el ruido. Como grupo difieren del resto (residuo loess p = 0.0012 / 0.012; percentil de ventana p = 2.2e-04 / 6.6e-04). Siguen difiriendo con la isla CpG como covariable, y dentro de los promotores con isla con el percentil (0.0010 / 0.023), aunque con el residuo la Rep 2 no pasa (0.13).
- DNMT3B y ZBTB33 quedan arriba de su ventana (percentil 0.97 / 0.96 y 0.94 / 0.95), pero tres controles desinflan el resultado:
  1. **DNMT3B lo sostienen 552 promotores con picos en HUES-8, casi todos sin isla.** Con los picos de Hep-G2, que son el 93% de su conjunto, queda en el medio de su ventana: percentil 0.73 / 0.68.
  2. **Dentro de los promotores con isla, ninguno de los dos se destaca.** DNMT3B queda en 0.76 / 0.61 y ZBTB33 en 0.78 / 0.91.
  3. **El grupo de metilación no difiere del resto ni siquiera con ellos** (residuo loess p = 0.10 / 0.60; sin ellos, 0.50 / 0.50). El extremo de percentil alto lo ocupan sobre todo factores de linaje (OTX2, CDX2, IKZF1, TBX5, POU4F2…), no lectores de metilación.
- **Ningún factor individual es significativo** con el p por permutación corregido por BH (mínimo p crudo 0.0056 / 0.0066).
- **El GSEA residualizado no da términos de metilación en ninguno de los dos extremos, y MLL1/2 desaparece** (p ajustado 0.62 / 0.99).
- **Lo único robusto del GSEA es «chromatin», y con signo opuesto al esperado.** Da NES +1.54 / +1.70 (p ajustado 0.0022 / 1.1e-05): los factores anotados a cromatina, que son en su mayoría factores de unión a secuencia, reducen el ruido algo menos que los demás factores de su misma cobertura.

**Corrección respecto de la primera versión de este handoff.** El residuo por cobertura era el de una recta robusta (`rlm`) en log10(n). La relación no es lineal y esa recta dejaba un sesgo sistemático: ~+1 en la franja de 1000–5000 promotores, donde están DNMT3B y ZBTB33, y de −1 a −2 en la de 6000–9000, donde están MAZ, SP1, SP2 y E2F1. Es decir, inflaba el contraste entre los dos grupos. Ahora el residuo principal es el de un loess robusto. Lo que cambió:

- el p del grupo de CpG no metilado se debilita (de 8e-05 / 2e-04 a 0.0012 / 0.012);
- el grupo de metilación deja de ser significativo también en la Rep 1;
- el GSEA cambia (ver p3).

El percentil de ventana no depende de ningún ajuste y no cambió. El residuo lineal se conserva como `residuo_rlm`, solo como referencia.
- **Recomendación:** descartar como panel de R3.7. A lo sumo una frase sobre los factores de CpG no metilado (ver el final).

## Método

- **Insumo.** Coeficiente de `unido` del modelo ajustado de `R/48` (`var_rank_sw ~ unido + log10(hek_tpm + 1) + mean_dnase`) para los 982 factores, por réplica (`p1_coeficientes_por_factor.tsv`). Promotores: 9300 / 11400.
- **Residuo por cobertura (p1).** Loess robusto (`family = "symmetric"`, span 0.5) del coeficiente sobre log10(n unidos), por réplica; el residuo es coef − ajuste. Es la línea de `p1_coef_vs_cobertura`.
  - **Forma de la relación.** Es plana (~0) hasta ~3000 promotores y cae cada vez más rápido por encima de ~4000. Con factores que unen el 40–80% de la library, el grupo no unido se va reduciendo a los promotores más atípicos (sin ocupación, inactivos, sin isla), que son los más ruidosos. A la izquierda la nube se abre en embudo, porque con pocos promotores unidos la estimación es imprecisa.
  - **Por qué no una recta.** Una recta en log10(n) no sigue esa forma. Mediana del residuo lineal (`rlm`) por tramo, Rep 1 / Rep 2:
    - 100–300 promotores: −1.3 / −1.6
    - 1000–5000: +0.8 a +1.2
    - 7000–8000: −2.0 / −0.9

    Con el loess, todas las medianas por tramo quedan entre −0.1 y +0.1.
  - **Dependencia residual de la cobertura.** Spearman residuo–n: −0.01 / −0.001 con el loess y −0.22 / −0.21 con la recta, contra −0.63 / −0.66 del coeficiente crudo.
  - **Referencia.** El residuo lineal queda guardado como `residuo_rlm`, y los tests de grupo también se reportan con él.
- **Ventana de cobertura (p1).** Mismo procedimiento que `remap_ruido_cobertura.R`:
  - Referencia: los otros factores con n unidos a ±25% del focal.
  - Percentil = fracción de la ventana con coeficiente ≤ al del focal (0 = mayor reducción del ruido).
  - p por permutación exacto: si el focal es intercambiable con los factores de su ventana, su rango es uniforme. p bilateral = 2·min(rango inferior, rango superior)/(n_ventana + 1). Es el límite de un test por permutaciones aleatorias, sin error de Monte Carlo.
  - BH sobre los 982 factores en cada réplica.
  - Con ventanas de 200–350 factores, el p mínimo posible es ~0.006, así que con BH sobre 982 nada puede pasar 0.05 salvo que muchos factores sean extremos a la vez.
  - Reproduce los percentiles de la exploración preliminar: DNMT3B 0.97 / 0.96, ZBTB33 0.94 / 0.95, MAZ 0.004 / 0.02, SP2 0.01 / 0.04, ZBTB14 0.04 / 0.02, E2F1 0.05 / 0.04, E2F4 0.04 / 0.06.
- **Tests de grupo (p2).** Wilcoxon de rangos, grupo contra el resto de los 982, por réplica. Reporto el rank-biserial (> 0 = valores más altos que el resto, es decir, menos reducción del ruido) y el desplazamiento de Hodges-Lehmann.
  - Comparaciones adicionales: metilación sin DNMT3B ni ZBTB33, y metilación contra CpG no metilado.
- **GSEA (p3).** Mismos parámetros que `remap_tf_noise_gsea`: `gseGO(ont = "CC")`, `set.seed(2024)` antes de la tanda de réplicas, `seed = TRUE`, el resto por defecto, y `pvalueCutoff = 1` solo para listar todo. Ranking: residuo loess, de mayor a menor.
- **Isla CpG (p4).** `CGI_frag`, la definición por composición del fragmento que usa todo R3. Es CGI el 71% / 69% de los promotores.
  - a) Modelo ajustado + `CGI_frag`.
  - b) Modelo ajustado dentro de cada estrato, con `var_rank_sw` global restringido al estrato. Exigí ≥ 50 unidos y ≥ 50 no unidos en el estrato en ambas réplicas.
    - Con isla: 6608 / 7900 promotores, 982 factores evaluables.
    - Sin isla: 2692 / 3500 promotores, 853 factores evaluables; DNMT1 queda afuera.
  - En cada versión recalculé residuo, ventana (n unidos dentro del estrato) y tests.
- **Procedencia (p5).** Los picos de `remap_tf_hits.tsv` vienen del bed no redundante de ReMap 2022. Su nombre es `TF:linea1,linea2,...` (un pico fusionado puede venir de varias líneas), así que separé las líneas.
  - El bed no redundante no guarda cuántos experimentos hay por línea, y el sitio de ReMap respondía 502. Como referencia cruzada conté los experimentos hg38 de ChIP-Atlas (`experimentList_chipatlas.tab`); no es exactamente el mismo conjunto que ReMap.
  - Para cada línea con > 100 promotores unidos reestimé el efecto usando solo sus picos y lo ubiqué en la ventana de cobertura.
- **Solapamiento (p6).** Conjuntos unidos y no unidos de DNMT3B y ZBTB33, y modelo ajustado con la categoría de unión (ninguno / solo DNMT3B / solo ZBTB33 / ambos).

## p2. Grupos: quién quedó afuera

| Grupo | Incluidos | Afuera |
|---|---|---|
| Metilación | DNMT1, DNMT3B, ZBTB33, MBD1, MBD2, MBD3, TET2 (7 de 13) | DNMT3A, UHRF1, TET1, TET3: no están en el catálogo ReMap 2022. MBD4 (18 / 25 promotores; solo Hep-G2) y MECP2 (17 / 24; solo interneuronas corticales): debajo del umbral de > 100 en ambas réplicas |
| CpG no metilado | MAZ, SP1, SP2, SP3, ZBTB14, E2F1, E2F4, NRF1, GABPA, KLF16 (10 de 12) | CXXC1 y KDM2B: no están en el catálogo ReMap 2022 |

- **No hay control positivo por dominio CXXC.** CXXC1 y KDM2B no están en ReMap.
- **DNMT1 tiene 102 / 119 promotores unidos** (solo K-562), así que su estimación es muy imprecisa. Su ventana tiene 22 factores.

## Resultados

### p1. Posición controlada por cobertura (modelo ajustado)

| Factor | n unidos | Coef. ajustado | Residuo loess | Percentil ventana | p perm. |
|---|---|---|---|---|---|
| DNMT3B | 4204 / 4877 | 0.77 / 0.92 | 1.46 / 1.56 | 0.97 / 0.96 | 0.077 / 0.091 |
| ZBTB33 | 5231 / 6135 | 0.13 / 0.14 | 1.68 / 1.79 | 0.94 / 0.95 | 0.13 / 0.10 |
| MBD1 | 4019 / 4676 | 0.18 / −0.61 | 0.79 / −0.04 | 0.82 / 0.54 | 0.37 / 0.93 |
| DNMT1 | 102 / 119 | 2.70 / 0.91 | 3.75 / 1.21 | 0.82 / 0.64 | 0.43 / 0.78 |
| TET2 | 3281 / 3741 | −0.43 / −0.66 | −0.04 / −0.29 | 0.52 / 0.40 | 0.97 / 0.81 |
| MBD2 | 5656 / 6594 | −1.95 / −3.08 | 0.06 / −0.98 | 0.52 / 0.27 | 0.96 / 0.54 |
| MBD3 | 1258 / 1570 | −0.83 / −1.41 | −1.03 / −1.56 | 0.22 / 0.14 | 0.46 / 0.30 |
| MAZ | 7665 / 9235 | −7.07 / −5.45 | −2.43 / −1.26 | 0.004 / 0.02 | 0.016 / 0.052 |
| SP2 | 6978 / 8313 | −6.31 / −5.04 | −2.63 / −1.50 | 0.01 / 0.04 | 0.025 / 0.095 |
| ZBTB14 | 4218 / 5098 | −3.25 / −3.68 | −2.55 / −2.91 | 0.04 / 0.02 | 0.095 / 0.046 |
| E2F1 | 7529 / 9074 | −5.67 / −5.35 | −1.22 / −1.26 | 0.05 / 0.04 | 0.11 / 0.081 |
| E2F4 | 6643 / 7867 | −5.65 / −4.90 | −2.42 / −1.70 | 0.04 / 0.06 | 0.091 / 0.13 |
| SP1 | 7582 / 9094 | −6.04 / −4.08 | −1.51 / 0.03 | 0.03 / 0.25 | 0.062 / 0.51 |
| NRF1 | 6853 / 8172 | −4.68 / −4.42 | −1.17 / −0.99 | 0.13 / 0.12 | 0.27 / 0.24 |
| SP3 | 6048 / 7149 | −4.74 / −2.19 | −2.27 / 0.42 | 0.08 / 0.60 | 0.16 / 0.81 |
| GABPA | 6850 / 8202 | −2.67 / −3.58 | 0.84 / −0.12 | 0.52 / 0.28 | 0.96 / 0.57 |
| KLF16 | 5205 / 6152 | −1.18 / −2.11 | 0.35 / −0.44 | 0.66 / 0.42 | 0.69 / 0.83 |

Con el residuo loess, los residuos de MAZ, SP2 y E2F1 son de la mitad o menos de lo que daba la recta. La diferencia entre las dos versiones es justamente el sesgo de la recta en la franja de cobertura alta.

- **Qué significa el extremo alto.** Los coeficientes ajustados de DNMT3B y ZBTB33 son cercanos a cero (+0.8 / +0.9 y +0.1 / +0.1). Lo que los distingue no es que aumenten el ruido: son de los pocos factores de su ventana que no lo bajan, en una franja de cobertura donde la mediana del coeficiente es −0.8 y −1.7 / −1.9.
- **Corrección múltiple.** Ningún factor de los 982 pasa BH < 0.1 en ninguna réplica. El p crudo mínimo es 0.0056 / 0.0066.
- **No son únicos en el extremo alto.** Entre los factores con más de 2000 promotores, DNMT3B y ZBTB33 están en los puestos 4 y 8 por percentil medio entre réplicas (13 y 20 entre los 982). Por residuo loess medio quedan en los puestos 53 y 69 de 982. Los acompañan sobre todo factores de linaje: KMT2B-D, OTX2, IKZF1, JUNB, CDX2, PAX3-FOXO1, TBX5, NFATC3, POU4F2 y TCF4 (`p1_ranking_residuo_ambas_replicas.tsv`). Esos factores no tienen relación con la metilación.
- **El extremo bajo, en cambio, es coherente: maquinaria de promotor activo.** RORC, TRIM25, CBX4, CDK9, MAZ, PHF8, CCNT2, SP2, ZBTB14, KDM4A y TAF1. Además de los de CpG no metilado aparecen P-TEFb (CDK9, CCNT2), TAF1 y lectores o borradores de marcas activas (PHF8, KDM4A).

### p2. Test de grupo

Wilcoxon, rank-biserial (p):

| Comparación | Residuo loess (principal) | Percentil ventana | Residuo rlm (referencia, sesgado) |
|---|---|---|---|
| Metilación vs. resto (n = 7) | 0.36 (0.10) / 0.12 (0.60) | 0.41 (0.063) / 0.14 (0.53) | 0.49 (0.026) / 0.21 (0.34) |
| Metilación sin DNMT3B ni ZBTB33 (n = 5) | 0.18 (0.50) / −0.17 (0.50) | 0.20 (0.44) / −0.18 (0.49) | 0.34 (0.20) / −0.06 (0.81) |
| CpG no metilado vs. resto (n = 10) | −0.60 (0.0012) / −0.46 (0.012) | −0.68 (2.2e-04) / −0.63 (6.6e-04) | −0.72 (8.3e-05) / −0.68 (2.1e-04) |
| Metilación vs. CpG no metilado | 0.80 (0.0073) / 0.51 (0.088) | 0.86 (0.0040) / 0.71 (0.017) | 0.94 (0.0015) / 0.83 (0.0054) |

- **CpG no metilado: se sostiene en ambas réplicas** con el residuo loess y con el percentil. Medianas del residuo loess: −1.90 / −1.10 contra 0.04 / 0.01 del resto. Mediana del percentil: 0.05 / 0.09.
- **Metilación: no se sostiene.** No llega a p < 0.05 con el residuo loess ni con el percentil en ninguna réplica. El único p < 0.05 (Rep 1, 0.026) era con el residuo lineal sesgado. Sin DNMT3B y ZBTB33 no hay ninguna señal. Los lectores menos específicos quedan en el medio, como en la exploración preliminar.
- **Metilación contra CpG no metilado.** Da p 0.0073 / 0.088 con el residuo loess y 0.0040 / 0.017 con el percentil, pero la arrastra el grupo de CpG no metilado.

### p3. GSEA sobre el ranking residualizado

Ranking por residuo loess (`p3_gsea_comparacion.tsv`, `p3_gsea_residuo_todos_los_terminos.tsv`, `p3_gsea_nes_comparacion.jpg`). NES (p ajustado):

| Término | AUC original | Coeficiente ajustado | Residuo por cobertura (loess) |
|---|---|---|---|
| MLL1/2 complex (GO:0044665) | −2.47 (0.0046) / −2.01 (0.081) | −1.76 (0.018) / −1.46 (0.18) | −1.46 (0.62) / −0.80 (0.99) |
| MLL1 complex (GO:0071339) | −2.37 (0.0052) / −1.86 (0.13) | −1.73 (0.018) / −1.40 (0.23) | −1.34 (0.74) / −0.66 (0.99) |
| methyltransferase complex (GO:0034708) | −1.71 (0.19) / −1.91 (0.10) | −1.55 (0.037) / −1.49 (0.12) | −1.23 (0.74) / −1.18 (0.82) |
| chromatin (GO:0000785) | −3.91 (5.8e-08) / −1.90 (0.10) | −1.26 (0.032) / −1.25 (0.11) | +1.54 (0.0022) / +1.70 (1.1e-05) |
| protein-containing complex (GO:0032991) | −3.79 (8.3e-08) / −1.65 (0.23) | −1.30 (0.021) / −1.23 (0.11) | +1.33 (0.097) / +1.55 (8.8e-04) |
| transcription regulator complex (GO:0005667) | −3.40 (5.4e-06) / −1.15 (0.59) | −1.30 (0.027) / −1.14 (0.41) | +1.38 (0.17) / +1.80 (1.5e-05) |
| nuclear protein-containing complex (GO:0140513) | −3.39 (2.1e-06) / −2.60 (0.0037) | −1.32 (0.018) / −1.32 (0.040) | +1.00 (0.79) / +1.05 (0.85) |
| protein kinase complex (GO:1902911) | −2.46 (0.0040) / −1.60 (0.23) | −1.69 (0.032) / −1.70 (0.040) | −1.74 (0.17) / −1.16 (0.85) |

- **En el extremo de menor ruido no queda ningún término.** Ninguno tiene NES negativo con p ajustado < 0.1 en ambas réplicas; el más bajo es vesículas intracelulares, 0.056 solo en la Rep 1.
- **Metilación.** Ningún término ligado a metilación o heterocromatina se acerca a la significancia en ningún extremo: methyltransferase complex, heterochromatin, NuRD, PcG, Sin3 y HDAC dan todos p ajustado > 0.65.
- **Chromatin cambia de signo y es robusto.** Pasa de NES negativo con los rankings sin controlar a +1.54 / +1.70 (p ajustado 0.0022 / 1.1e-05). Protein-containing complex también cumple el criterio de p < 0.1 en ambas réplicas (0.097 / 8.8e-04).
  - **Qué significa.** El término incluye 476 de los ~950 factores mapeados, sobre todo factores de unión a secuencia. Lo que dice es que, a igual cobertura, esos factores reducen el ruido algo menos que el resto: cofactores, maquinaria de Pol II y P-TEFb, proteínas de unión a ARN.
  - **Magnitud.** La diferencia es chica: medianas del residuo ~+0.2 a +0.7 contra ~−0.1 a −0.35 puntos de rango.
  - **No lo arrastran los factores de baja cobertura.** Se mantiene en los tres tercios de cobertura con el residuo (p 0.0015, 4e-04 y 0.011 / 4e-05, 0.0037 y 9e-04), y con el percentil en los dos tercios inferiores (`p3_chromatin_por_tercil_cobertura.tsv`).
  - **Lectura.** No es una señal de metilación. Es compatible con que los cofactores sin unión directa al ADN se recluten preferentemente a promotores activos, de modo que su ChIP marca actividad, mientras que los factores de secuencia se unen por motivo también a promotores menos activos. No lo testeé más.
- **MLL1/2 no sobrevive.** Queda con p ajustado 0.62 / 0.99. Es más estricto que la ocupación como covariable de p7 (0.058 / 0.43) y confirma la conclusión del handoff anterior: a nivel complejo, el enriquecimiento de MLL1/2 era cobertura. Lo que queda es lo de WDR5 y KMT2B por miembro.
- **Advertencia.** fgsea pondera por el valor de la métrica, y los residuos de los factores de baja cobertura son más dispersos (SD 1.8 / 1.9 en el tercio inferior contra 1.2 / 1.0 en el superior). Los extremos de esta lista los ocupan muchos factores de pocos promotores (ICE2, ZNF436, ZNF597, TP73…).
- **Versión con la recta.** Rankeando por el residuo lineal sesgado, la primera versión daba MLL1/2 en 0.22 / 0.24 y chromatin en +1.27 (0.24) / +1.37 (0.077). El sesgo de la recta favorecía a los factores de cobertura alta en el extremo de menor ruido.

### p4. Control por isla CpG

**Premisa a corregir.** DNMT3B y ZBTB33 no unen preferentemente promotores sin isla:

| | Unidos con isla | No unidos con isla |
|---|---|---|
| DNMT3B | 0.845 / 0.842 | 0.60 / 0.58 |
| ZBTB33 | 0.844 / 0.833 | 0.54 / 0.53 |
| MAZ | 0.82 / 0.81 | 0.21 / 0.20 |
| SP2 | 0.85 / 0.85 | 0.28 / 0.28 |
| ZBTB14 | 0.955 / 0.953 | 0.51 / 0.48 |
| Mediana de los 982 | 0.865 / 0.856 | |
| Library | 0.71 / 0.69 | |

- Unen más CGI que la library, aunque algo menos que el factor típico.
- Lo que los diferencia de MAZ o SP2 no es el conjunto unido sino el no unido: el grupo no unido de MAZ es 80% sin isla y el de DNMT3B, 40%.

**Percentil en la ventana** (`p4_percentil_ventana_miembros.jpg`, `p4_miembros_grupos_cgi.tsv`):

| Factor | Ajustado | + isla CpG | Solo con isla | Solo sin isla |
|---|---|---|---|---|
| DNMT3B | 0.97 / 0.96 | 0.90 / 0.95 | 0.76 / 0.61 | 0.91 / 0.99 |
| ZBTB33 | 0.94 / 0.95 | 0.92 / 0.95 | 0.78 / 0.91 | 0.95 / 0.82 |
| MAZ | 0.004 / 0.02 | 0.004 / 0.03 | 0.01 / 0.04 | 0.00 / 0.31 |
| SP2 | 0.01 / 0.04 | 0.01 / 0.05 | 0.05 / 0.12 | 0.03 / 0.28 |
| ZBTB14 | 0.04 / 0.02 | 0.09 / 0.05 | 0.08 / 0.02 | 0.14 / 0.50 |
| E2F1 | 0.05 / 0.04 | 0.08 / 0.06 | 0.003 / 0.16 | 0.65 / 0.46 |
| E2F4 | 0.04 / 0.06 | 0.05 / 0.09 | 0.07 / 0.21 | 0.28 / 0.31 |

Coeficiente ajustado dentro de cada estrato, Rep 1 / Rep 2:

| | Con isla | Sin isla |
|---|---|---|
| DNMT3B | 0.96 / 0.90 (3551 / 4107 unidos) | 0.32 / 1.32 (653 / 770) |
| ZBTB33 | 0.87 / 1.62 | 0.30 / −1.53 |

**Tests de grupo** (rank-biserial, p). El residuo loess y el percentil se recalculan dentro de cada versión:

| Versión | Estadístico | Metilación vs. resto | Sin DNMT3B ni ZBTB33 | CpG no metilado vs. resto |
|---|---|---|---|---|
| Ajustado | residuo loess | 0.36 (0.10) / 0.12 (0.60) | 0.18 (0.50) / −0.17 (0.50) | −0.60 (0.0012) / −0.46 (0.012) |
| | percentil | 0.41 (0.063) / 0.14 (0.53) | 0.20 (0.44) / −0.18 (0.49) | −0.68 (2.2e-04) / −0.63 (6.6e-04) |
| + isla CpG | residuo loess | 0.44 (0.046) / 0.16 (0.46) | 0.33 (0.20) / −0.09 (0.73) | −0.55 (0.0027) / −0.42 (0.022) |
| | percentil | 0.48 (0.028) / 0.18 (0.42) | 0.33 (0.20) / −0.12 (0.64) | −0.63 (5.9e-04) / −0.59 (0.0014) |
| Solo con isla | residuo loess | 0.34 (0.12) / 0.10 (0.65) | 0.32 (0.21) / −0.03 (0.92) | −0.58 (0.0016) / −0.28 (0.13) |
| | percentil | 0.38 (0.084) / 0.10 (0.65) | 0.31 (0.24) / −0.07 (0.78) | −0.60 (0.0010) / −0.42 (0.023) |
| Solo sin isla (sin DNMT1) | residuo loess | 0.22 (0.34) / −0.12 (0.60) | −0.05 (0.87) / −0.54 (0.061) | −0.34 (0.062) / −0.17 (0.35) |
| | percentil | 0.26 (0.27) / −0.07 (0.76) | −0.04 (0.88) / −0.51 (0.078) | −0.39 (0.034) / −0.21 (0.25) |

**Veredicto del control de isla CpG**

- **Covariable.** La isla CpG como covariable casi no cambia la posición de ninguno de los dos grupos.
- **Dentro de los promotores con isla:**
  - DNMT3B y ZBTB33 dejan de ser extremos. Quedan en la mitad superior de su ventana (0.76 / 0.61 y 0.78 / 0.91), lejos del 5% superior, con p por permutación 0.49 / 0.79 y 0.45 / 0.19. DNMT3B en particular no se destaca en la Rep 2.
  - El grupo de metilación no difiere del resto (p 0.12 / 0.65 con el residuo y 0.084 / 0.65 con el percentil).
- **Dentro de los promotores sin isla:** DNMT3B queda alto (0.91 / 0.99; p 0.18 / 0.038). Ese resultado lo explica p5: sale casi entero de los picos de HUES-8.
- **El hallazgo de metilación se desinfla.** No se reduce del todo a la clase de promotor: con isla como covariable, la posición casi no cambia. Pero dentro de los promotores con isla, DNMT3B y ZBTB33 no se destacan de los factores de cobertura comparable, y como grupo no hay señal en ninguna versión.
- **El grupo de CpG no metilado sobrevive dentro de los promotores con isla, más débil en la Rep 2.**
  - Con el percentil da p 0.0010 / 0.023; con el residuo loess, 0.0016 / 0.13. La Rep 2 no pasa con el residuo.
  - MAZ, SP2 y ZBTB14 siguen en el 1–12% de mayor reducción entre factores de cobertura comparable; E2F1 solo en la Rep 1.
  - Dentro de los promotores sin isla el efecto se pierde (p 0.034–0.062 / 0.25–0.35). Esos factores unen pocos promotores sin isla (189–1805), así que la potencia baja mucho.

### p5. Procedencia de los picos

Líneas celulares en ReMap 2022 (bed no redundante), promotores de la library unidos y efecto por línea (percentil en la ventana de cobertura con los picos de esa línea sola):

| Factor | Líneas en ReMap | Promotores por línea (fracción del conjunto) | Efecto por línea: percentil Rep 1 / Rep 2 | Experimentos en ChIP-Atlas (hg38) |
|---|---|---|---|---|
| DNMT3B | 2: Hep-G2, HUES-8 | Hep-G2 6773 (0.93); HUES-8 552 (0.08) | Hep-G2 0.73 / 0.68 (coef −0.14 / −0.30); HUES-8 1.00 / 1.00 (coef +5.6 / +6.9; 7–8% con isla) | 35 en 9 tipos celulares (hESC H1 12, MCF10CA1a 6…) |
| ZBTB33 | 8: A-549, hígado, GM12878, HCT-116, K-562, SK-N-SH, Hep-G2, MCF-7 | 311–5642 por línea; ninguna supera el 60% | hígado 1.00 / 0.95; Hep-G2 0.97 / 0.76; SK-N-SH 0.97 / 0.86; GM12878 0.54 / 0.95; A-549 0.78 / 0.74; K-562 0.63 / 0.68; HCT-116 0.50 / 0.52 | 43 en 12 (293: 14, GM12878: 8…) |
| MAZ | 8: HEK293, IMR-90, K-562, HeLa-S3, Hep-G2, MCF-7, A-549, GM12878 | 5412–13864 por línea | HEK293 0.02 / 0.08; HeLa-S3 0.03 / 0.29; IMR-90 0.07 / 0.13; K-562 0.10 / 0.43; otras 0.10–0.31 | 36 en 12 |
| SP2 | 5: HEK293, HEK293T, K-562, WA01, Hep-G2 | HEK293 13961 (0.99) | HEK293 0.02 / 0.07; HEK293T 0.25 / 0.40; K-562, WA01 y Hep-G2 0.31–0.68 | 7 en 4 |
| ZBTB14 | 2: HEK293, Hep-G2 | HEK293 8415 (0.91); Hep-G2 5684 (0.61) | HEK293 0.03 / 0.02; Hep-G2 0.11 / 0.11 | 2 (solo Hep G2) |

Advertencias:

- **DNMT3B.** Tiene picos de dos líneas, pero en la práctica su efecto viene de una sola.
  - Los picos de Hep-G2 dan un factor sin nada particular.
  - Todo el extremo alto lo producen 552 promotores con picos de HUES-8 (célula madre embrionaria), el 92% sin isla.
  - Es un conjunto chico y sesgado a promotores sin isla, que ya sabemos que son más ruidosos. Encaja con la conocida unión de DNMT3B a sitios metilados en ESC, pero como evidencia del eje de metilación se reduce a «los promotores sin isla unidos en HUES-8 no bajan el ruido».
- **ZBTB14.** Tiene dos líneas y el 91% de sus promotores viene de HEK293. Solo 2 experimentos en ChIP-Atlas.
- **SP2.** El 99% viene de HEK293, y el efecto fuerte es solo con esos picos; en las otras líneas queda en el medio.
- **ZBTB33.** Es el único de los cinco con una base amplia (8 líneas). Su posición alta es consistente pero moderada en casi todas: 6 de 8 líneas por encima de 0.5 en ambas réplicas. Las otras dos son HCT-116 (0.50 / 0.52) y MCF-7 (0.38 / 0.69, con solo 161 / 178 promotores). Ninguna línea sola lo explica.
- **Lado de CpG no metilado.** Para los tres factores (MAZ, SP2, ZBTB14), el efecto es más fuerte con los picos de HEK293, la misma línea que la library, aun ajustando por la actividad CAGE en HEK293 y la DNase. Tiene sentido: la ocupación en el contexto endógeno de la propia línea es la más informativa.
- **Número de experimentos en ReMap.** No pude obtenerlo: el bed no redundante no lo guarda y el sitio de ReMap respondía con error. Los conteos de ChIP-Atlas son una aproximación.

### p6. Solapamiento DNMT3B / ZBTB33

| | Rep 1 | Rep 2 |
|---|---|---|
| Unidos por DNMT3B / ZBTB33 | 4204 / 5231 | 4877 / 6135 |
| Unidos por ambos (esperado si fueran independientes) | 3287 (2365) | 3753 (2625) |
| Fracción de DNMT3B también en ZBTB33 | 0.78 | 0.77 |
| Fracción de ZBTB33 también en DNMT3B | 0.63 | 0.61 |
| Jaccard, unidos / no unidos | 0.54 / 0.52 | 0.52 / 0.54 |
| Odds ratio (Fisher) | 5.8 | 5.8 |
| Puesto de ZBTB33 entre los 981 por Jaccard con DNMT3B | 140 | 152 |

- **El solapamiento es grande pero no particular.** Entre factores de cobertura alta casi todos se solapan así. ZBTB33 no está entre los más parecidos a DNMT3B (puesto 140 / 152).
- **El más parecido a DNMT3B es MBD1** (Jaccard 0.80 / 0.80). MBD1 queda en 0.82 / 0.54, así que tampoco es independiente de DNMT3B.

Efecto por categoría, en el modelo ajustado y frente a los promotores sin ninguno de los dos (`p6_efecto_por_categoria.jpg`):

| Categoría | n | Fracción con isla | Coef. (IC 95%) | p |
|---|---|---|---|---|
| Solo DNMT3B | 917 / 1124 | 0.65 / 0.66 | 1.36 (−0.78, 3.50) / 1.50 (−0.45, 3.44) | 0.21 / 0.13 |
| Solo ZBTB33 | 1944 / 2382 | 0.75 / 0.73 | 0.29 (−1.40, 1.99) / 0.27 (−1.26, 1.80) | 0.73 / 0.73 |
| Ambos | 3287 / 3753 | 0.90 / 0.90 | 0.73 (−0.97, 2.43) / 0.85 (−0.72, 2.42) | 0.40 / 0.29 |
| DNMT3B, sin los promotores de ZBTB33 | 917 / 1124 | | 3.13 (0.85, 5.41) / 2.68 (0.62, 4.74) | 0.0071 / 0.011 |
| ZBTB33, sin los promotores de DNMT3B | 1944 / 2382 | | 1.12 (−0.68, 2.93) / 0.91 (−0.71, 2.53) | 0.22 / 0.27 |

- **Los dos efectos no los sostienen los mismos promotores.**
- **DNMT3B.** Lo sostienen sobre todo los promotores que une solo DNMT3B, que son los más pobres en isla (65%) e incluyen los de HUES-8. Excluyendo los de ZBTB33, su efecto sube a +3.1 / +2.7.
- **ZBTB33.** Su efecto es débil en todas las particiones. Su percentil alto viene de ser apenas no negativo en una franja de cobertura donde la mediana es −1.7 / −1.9.
- **No son una réplica.** No son dos marcas del mismo grupo de promotores, pero tampoco dos evidencias de un mismo mecanismo. Cada uno tiene su explicación más simple: DNMT3B, un subconjunto de promotores sin isla de una línea; ZBTB33, un efecto chico y difuso.

## Veredicto

1. **Eje de metilación: no se sostiene como hallazgo.**
   - La observación preliminar se reproduce (DNMT3B y ZBTB33 en el 3–6% superior de su ventana), pero:
     - ningún factor es significativo individualmente con BH;
     - el grupo de metilación no difiere del resto sin esos dos factores;
     - dentro de los promotores con isla, ninguno de los dos se destaca;
     - el extremo de DNMT3B sale de 552 promotores sin isla con picos de HUES-8;
     - el extremo alto lo ocupan sobre todo factores de linaje sin relación con la metilación.
   - Cómo redactarlo si se menciona: DNMT3B y ZBTB33 «están entre los pocos factores de su cobertura cuya unión no se asocia a una reducción del ruido», nunca «se asocian a mayor ruido».
2. **Lado de CpG no metilado: se sostiene con fuerza moderada.**
   - Como grupo predefinido, MAZ, SP1, SP2, SP3, ZBTB14, E2F1, E2F4, NRF1, GABPA y KLF16 reducen el ruido más que los factores de cobertura comparable, en ambas réplicas (residuo loess p 0.0012 / 0.012; percentil p 2.2e-04 / 6.6e-04).
   - Sobrevive a la isla CpG como covariable (residuo 0.0027 / 0.022). Dentro de los promotores con isla sobrevive con el percentil (0.0010 / 0.023), pero con el residuo solo en la Rep 1 (0.0016 / 0.13). No es solo la clase de promotor, pero en la Rep 2 se debilita.
   - Límites:
     - es un efecto de grupo; individualmente nadie pasa BH;
     - no hay control positivo CXXC;
     - el extremo bajo lo comparten con maquinaria general de promotor activo (CDK9, CCNT2, TAF1, PHF8), así que no se puede atribuir específicamente a «leer CpG no metilado»;
     - para MAZ, SP2 y ZBTB14, la señal es más fuerte con los picos de HEK293.
3. **GSEA residualizado.**
   - No aparecen términos de metilación en ningún extremo, y MLL1/2 desaparece (0.62 / 0.99).
   - Lo único robusto es «chromatin» con signo positivo (0.0022 / 1.1e-05): a igual cobertura, los factores de unión a secuencia reducen el ruido algo menos que los cofactores. Es lo opuesto a «cromatina → menos ruido», chico en magnitud, y no tiene relación con la metilación.

## Recomendación para R3.7

- **No hacer un panel.** El eje de metilación, que era la motivación, no sobrevive a los controles. Un panel con el lado de CpG no metilado solo sería un efecto de grupo modesto y difícil de separar de «maquinaria de promotor activo».
- **Mención en el texto, opcional y con cautela.** Una frase posible, siempre con los números Rep 1 / Rep 2: «Entre factores de cobertura comparable, los que unen CpG no metilado (MAZ, SP2, ZBTB14, E2F1/4, entre otros) están entre los de mayor reducción del ruido, también dentro de los promotores con isla CpG (Wilcoxon de grupo sobre el percentil de cobertura, p = 0.0010 / 0.023)». Esa frase es consistente con el resultado de islas CpG del capítulo, pero no lo confirma de forma independiente: la diferencia entre factores se explica en parte por cuáles promotores no unen, y dentro del estrato con isla el efecto se debilita en la Rep 2 (con el residuo, p 0.13).
- **Descartar** la mención de DNMT3B y ZBTB33 como evidencia de metilación.
- **Para el resto de R3.7.** Es un argumento más para controlar la cobertura por defecto en el GSEA: con el ranking crudo, «chromatin» y «nuclear protein-containing complex» salen enriquecidos en el extremo de menor ruido solo por cobertura, y «chromatin» incluso se invierte al controlarla. El control tiene que ser no lineal: la ventana de cobertura o un loess, no una recta en log(n).

## Archivos

- `p1_residuo_ventana_por_factor.tsv`: para los 982 factores y ambas réplicas, coeficiente, ajuste loess, residuo loess (y en z), residuo rlm de referencia, n de la ventana, percentil, mediana de la ventana, p por permutación y p BH.
- `p1_resumen_cobertura.tsv`: Spearman con la cobertura antes y después de residualizar, pendiente de la recta de referencia y conteo BH.
- `p1_miembros_grupos.tsv`: lo mismo para los 17 miembros de los grupos.
- `p1_ranking_residuo_ambas_replicas.tsv`: los 982 en formato ancho por réplica, ordenados por residuo medio.
- `p1_coef_vs_cobertura.{jpg,pdf}`: coeficiente ajustado contra n unidos (log10), con los dos grupos destacados y la curva loess, por réplica.
- `p1_residuo_vs_cobertura.{jpg,pdf}`: residuo contra n unidos.
- `p2_grupos_definicion.tsv`: miembros de cada grupo, n por réplica, si están incluidos y por qué no.
- `p2_test_grupos.tsv`: Wilcoxon de los grupos con las tres versiones del estadístico (residuo loess, percentil y residuo rlm de referencia).
- `p2_residuo_por_grupo.{jpg,pdf}`: distribución del residuo y del percentil por grupo.
- `p3_gsea_residuo_todos_los_terminos.tsv`, `p3_gsea_comparacion.tsv`, `p3_gsea_nes_comparacion.{jpg,pdf}`: GSEA sobre el residuo loess y su comparación con la AUC original y con el coeficiente ajustado.
- `p3_chromatin_por_tercil_cobertura.tsv`: chequeo del NES de «chromatin», anotados contra no anotados dentro de cada tercio de cobertura.
- `p4_residuo_ventana_por_factor_cgi.tsv`: residuo y ventana de los 982 en las cuatro versiones (ajustado, + isla, con isla, sin isla), más la fracción con isla de los unidos y de los no unidos.
- `p4_miembros_grupos_cgi.tsv`, `p4_test_grupos_cgi.tsv`, `p4_n_por_estrato.tsv`: miembros, tests de grupo y n de cada estrato.
- `p4_percentil_ventana_miembros.{jpg,pdf}`: percentil de cada miembro en las cuatro versiones.
- `p4_coef_vs_cobertura_estratos.{jpg,pdf}`: coeficiente contra cobertura dentro de cada estrato.
- `p5_procedencia_lineas_remap.tsv`: líneas celulares, promotores por línea y exclusivos de cada línea.
- `p5_efecto_por_linea_celular.tsv`: efecto y percentil por línea.
- `p5_experimentos_chipatlas.tsv`: experimentos de ChIP-Atlas por tipo celular.
- `p6_solapamiento.tsv`, `p6_jaccard_top15.tsv`: solapamiento de los conjuntos unidos y no unidos, y factores más parecidos a cada uno.
- `p6_efecto_por_categoria.{tsv,jpg,pdf}`: efecto por categoría de unión.

---

# p7. P-TEFb y la maquinaria de pausa

## Método

Usa la misma infraestructura que p1–p4: coeficiente ajustado, residuo loess, percentil en la ventana de ±25% y p por permutación con BH. Cada variante de modelo o estrato recalcula la ventana contra los 982 factores bajo esa misma variante.

**Grupos a priori**, con la tabla completa en `p7_grupos_definicion.tsv`:

- **Pausa / elongación (13).** Son los del pedido presentes: CDK9, CCNT2, NELFA, NELFCD, NELFE, SUPT5H, AFF4, BRD4 y MLLT3. Agregué los del mismo complejo que encontré:
  - AFF1, ELL2 y MLLT1, del SEC;
  - HEXIM1, del 7SK snRNP.
- **Ausentes del catálogo de ReMap sobre la library:** CCNT1, NELFB, SUPT4H1, CDK12, CDK13, ELL3, EAF1/2, HEXIM2, LARP7 y MEPCE. ELL está pero queda bajo el umbral (49 / 58 promotores).
- **Preiniciación (10).**
  - TFIID: TBP, TAF1, TAF2, TAF3, TAF7 y TAF9B.
  - TFIIA: GTF2A2.
  - TFIIB: GTF2B.
  - TFIIF: GTF2F1.
  - TFIIH: CDK7.
  - Excluí TAF15: es una proteína FET y no una subunidad de TFIID, aunque comparta el nombre.
- **Mediador (5):** MED, MED1, MED12, MED25 y MED26. Va solo como referencia.

**Control de «maquinaria en general».** Usé tres versiones:

1. **Ocupación general.** El modelo ajustado más log10(n de los otros 981 factores unidos al promotor + 1), para los 982 factores.
2. **Modelo conjunto a nivel de promotor.** Un solo lm con los indicadores de unión de CDK9, CCNT2, NELFA, NELFE, SUPT5H, AFF4, BRD4, TAF1, TBP, GTF2B y GTF2F1, más las covariables. Da el efecto parcial de cada uno.
3. **CDK9 dentro y fuera del núcleo de preiniciación.** El núcleo es tener unidos a la vez TAF1, TBP y GTF2B.

## Resultados

**Miembros.** Coeficiente ajustado y percentil en la ventana (Rep 1 / Rep 2):

| Factor | Coeficiente | Percentil | p perm. | Ocupación general: percentil | Con isla: percentil | Sin isla: percentil |
|---|---|---|---|---|---|---|
| CDK9 | −6.1 / −6.3 | 0.02 / 0.00 | 0.038 / 0.0066 | 0.01 / 0.00 | 0.06 / 0.00 | 0.06 / 0.23 |
| CCNT2 | −4.2 / −3.0 | 0.00 / 0.05 | 0.0093 / 0.10 | 0.00 / 0.04 | 0.00 / 0.03 | 0.36 / 0.83 |
| NELFA | −7.2 / −4.2 | 0.00 / 0.15 | 0.0056 / 0.30 | 0.00 / 0.18 | 0.00 / 0.35 | 0.00 / 0.09 |
| NELFE | −5.5 / −4.5 | 0.08 / 0.13 | | 0.09 / 0.27 | 0.03 / 0.24 | 0.31 / 0.39 |
| SUPT5H | −5.4 / −3.4 | 0.08 / 0.43 | | 0.11 / 0.77 | 0.02 / 0.51 | 0.29 / 0.91 |
| BRD4 | −3.2 / −2.9 | 0.46 / 0.63 | | 0.26 / 0.28 | 0.10 / 0.01 | 0.23 / 0.70 |
| AFF4 | −3.3 / −2.9 | 0.36 / 0.48 | | 0.48 / 0.64 | 0.44 / 0.33 | 0.51 / 0.93 |
| TAF1 | −5.9 / −5.4 | 0.04 / 0.03 | | 0.08 / 0.10 | 0.57 / 0.82 | 0.01 / 0.01 |
| TBP | −6.5 / −3.5 | 0.01 / 0.35 | | 0.01 / 0.71 | 0.01 / 0.93 | 0.06 / 0.46 |
| GTF2B | −4.7 / −3.6 | 0.07 / 0.19 | | 0.07 / 0.21 | 0.07 / 0.18 | 0.54 / 0.72 |

Ningún factor pasa BH sobre 982: el mínimo de CDK9 es p = 0.038 / 0.0066 con BH = 1 / 1.

CDK9 es el único que se mantiene en el extremo inferior en ambas réplicas en todas las versiones, salvo el estrato sin isla de la Rep 2 (0.23). Ese estrato tiene ~1000–1300 promotores unidos. CCNT2 se sostiene en el estrato con isla, pero casi no une promotores sin isla (267 / 318). NELFA es extrema solo en la Rep 1.

**Test de grupo.** Wilcoxon sobre el percentil de ventana, en el modelo ajustado:

| Comparación | Medianas | p |
|---|---|---|
| Pausa vs. resto | 0.36 / 0.32 vs. 0.49 | 0.049 / 0.10 |
| Pausa sin CDK9 ni CCNT2 vs. resto | 0.41 / 0.43 vs. 0.49 | 0.26 / 0.43 |
| Preiniciación vs. resto | 0.16 / 0.40 vs. 0.49 | 0.17 / 0.97 |
| Preiniciación sin GTF2A2 vs. resto | 0.09 / 0.35 vs. 0.49 | 0.043 / 0.60 |
| **Pausa vs. preiniciación** | 0.36 / 0.32 vs. 0.16 / 0.40 | **0.83 / 0.31** |
| Pausa vs. preiniciación sin GTF2A2 | | 0.89 / 0.50 |
| Mediador vs. resto | 0.51 / 0.58 vs. 0.49 | 0.98 / 0.44 |

- **Pausa vs. preiniciación:** la diferencia de Hodges-Lehmann es −0.02 (IC −0.32 a 0.33) / −0.13 (−0.42 a 0.17).
- **Con el residuo loess,** pausa vs. resto da 0.13 / 0.44 y pausa vs. preiniciación 0.88 / 0.48.
- **Con la covariable de ocupación** el cuadro es el mismo: pausa vs. resto 0.039 / 0.20, pausa vs. preiniciación 0.88 / 0.40.
- **Estrato con isla:** pausa vs. resto 0.024 / 0.029 y pausa vs. preiniciación 0.44 / 0.044. En la Rep 2 la preiniciación sube a 0.74, sobre todo por TAF1 y TBP.
- **Estrato sin isla:** pausa vs. resto 0.067 / 0.76.

Lectura de los tests de grupo:

- **Como grupo, la pausa no se distingue del resto** de forma reproducible. Sin sus dos miembros extremos no queda nada (0.26 / 0.43). El resto del grupo tira hacia abajo con magnitudes chicas, pero ese desplazamiento también lo tiene la preiniciación.
- **La comparación directa pausa vs. preiniciación no los distingue** (0.83 / 0.31). Tampoco tiene potencia para hacerlo: son 13 contra 10 factores y el IC de la diferencia abarca ±0.3 de percentil.
- **Ni siquiera la preiniciación como grupo queda en el extremo** en la Rep 2. Así que la premisa de que toda la maquinaria basal está en el extremo de menor ruido vale para TAF1 y GTF2B, para TBP solo en la Rep 1, y no vale como grupo.

**Modelo conjunto a nivel de promotor.** Efectos parciales con IC 95%, Rep 1 / Rep 2:

| Indicador | Efecto parcial | p |
|---|---|---|
| CDK9 | −3.1 (−5.1 a −1.0) / −4.7 (−6.6 a −2.9) | 0.0036 / 4.1e-07 |
| CCNT2 | −2.6 / −1.8 | 3.3e-04 / 0.0081 |
| BRD4 | −2.1 / −1.9 | 0.0079 / 0.0067 |
| NELFA | −4.5 / −1.0 | 1.6e-04 / 0.33 |
| TAF1 | −0.4 / −2.9 | 0.76 / 0.010 |
| TBP | −3.1 / +0.8 | 0.012 / 0.43 |
| GTF2B | −0.4 / −0.2 | 0.73 / 0.87 |
| GTF2F1 | +1.5 / +0.8 | 0.17 / 0.43 |
| NELFE, SUPT5H, AFF4 | entre −1.2 y +2.2, inconsistentes | |

- **CDK9 dentro del núcleo de preiniciación:** −5.6 / −5.2 (p 0.022 / 0.016). Ojo con este número: el 97% de los promotores con TAF1 + TBP + GTF2B también une CDK9, así que la comparación es contra solo 138 / 183 promotores. De ahí el IC tan amplio (−10.4 a −0.8 / −9.5 a −1.0).
- **CDK9 fuera del núcleo:** −4.3 / −5.5 (p 8.6e-05 / 1.1e-08).
- **Solapamiento:** el 97% de los promotores de CDK9 une TAF1, el 93% TBP y el 82% GTF2B. El Jaccard con todos ellos va de 0.79 a 0.85.

Esta es la evidencia más directa a favor de algo específico. Con las demás uniones en el mismo modelo, los indicadores de P-TEFb (CDK9, CCNT2) y BRD4 conservan un efecto parcial negativo en ambas réplicas. Los de preiniciación no: TAF1, TBP y GTF2B no son consistentes entre réplicas. Dos límites:

- los indicadores son muy colineales (Jaccard ~0.8), así que los efectos parciales son frágiles;
- «unido por CDK9» puede estar marcando los picos de ChIP más intensos más que una función de P-TEFb.

**Estratos por isla CpG.**

- CDK9 en el estrato con isla: 0.06 / 0.00. Sin isla: 0.06 / 0.23.
- Su efecto no se reduce a que une promotores con isla: el 85% de sus promotores tiene isla, contra 71% de la library.
- En el estrato sin isla de la Rep 2 se debilita, con ~1000–1300 unidos.

**Procedencia** (`p7_p8_procedencia_lineas_remap.tsv`, `p7_p8_efecto_por_linea_celular.tsv`):

- **CCNT2: una sola línea (K-562).** En ChIP-Atlas figura 1 experimento. ⚠ Todo su efecto viene de un único experimento en una línea que no es HEK293.
- **CDK9: 10 líneas,** con 83 experimentos en ChIP-Atlas en 21 tipos celulares. La unión de los conjuntos es de 14697 promotores; la de la library, de 7067 / 8456. Las líneas que más cubren son MM1-S (86%), HEK293T (69%) y HCT-116 (50%).
- **El efecto de CDK9 no es homogéneo entre líneas.** Coeficiente y percentil por línea, Rep 1 / Rep 2:

| Línea | Coeficiente | Percentil |
|---|---|---|
| MM1-S | −5.9 / −4.9 | 0.02 / 0.05 |
| HEK293T | −5.3 / −4.0 | 0.006 / 0.06 |
| HCT-116 | −1.5 / −2.8 | 0.18 / 0.05 |
| A-375 | 0.0 / −1.5 | 0.47 / 0.09 |
| BT-474 | +0.6 / +1.0 | |
| MV4-11 | +1.0 / +0.4 | |
| P493-6 | +1.3 / +1.2 | |
| Kelly | +3.5 / +1.6 | 0.99 / 0.90 |

  El efecto lo sostienen los dos conjuntos más amplios, uno de ellos en HEK293T (la línea de la library). En las líneas con conjuntos chicos (~1000–2500 promotores) es nulo o positivo.

## Veredicto p7

- **El grupo de pausa no tiene un efecto propio.** No se sostiene sin CDK9 y CCNT2 (0.26 / 0.43). Tampoco se distingue de la preiniciación (0.83 / 0.31), y esa comparación no tiene potencia para distinguirlos.
- **Lo que queda es CDK9 como factor individual.** Es el extremo de su ventana en ambas réplicas y resiste:
  - la ocupación general;
  - la TATA-box: con la covariable, 0.02 / 0.00 (p8);
  - el estrato con isla;
  - el modelo conjunto con la preiniciación.
- **Límites de CDK9:**
  - no pasa BH;
  - viene sobre todo de dos experimentos, MM1-S y HEK293T;
  - su efecto parcial depende de indicadores muy colineales.
- **CCNT2 viene de un solo experimento** en K-562.
- **Lectura honesta:** no hay evidencia de grupo para «la maquinaria de pausa reduce el ruido». Tampoco alcanza para separarla de «los promotores ocupados por la maquinaria transcripcional son menos ruidosos». CDK9 es el factor individual que mejor sobrevive a ese control, pero de un factor no se puede inferir un mecanismo de pausa.
- **Recomendación: mención en el texto, no panel.** Una frase posible: «Entre los factores de la maquinaria basal, CDK9 (P-TEFb) es el de mayor reducción del ruido respecto de factores de cobertura comparable (percentil 0.02 / 0.00), y el efecto se mantiene al controlar la ocupación general del promotor y la unión del complejo de preiniciación; como grupo, en cambio, los factores de pausa y elongación no se distinguen de los de preiniciación». No vincularlo a la frecuencia de ráfagas como resultado, a lo sumo como especulación en la discusión.

---

# p8. GTF2A2 (TFIIA) y la TATA-box; caso lateral SNAPC1

## Método

- **TATA-box canónica:** `TATA_EPD` de `prom_df`, idéntica a la clase «Canónica» de `tata_scan.R` (el script de `tata_ruido_especificidad` lo verifica con un `stopifnot`).
- **Actividad alta:** percentil 75 de `mean` por réplica, igual que en `ROC_TATA_expr`.
- **Odds ratio de TATA:** Fisher, unidos contra no unidos, calculado para los 982 factores para ubicar a GTF2A2.
- **Modelos:** el ajustado más `TATA_EPD`, y los estratos de actividad con y sin TATA, todos para los 982 factores y con ventana. Los estratos sin y con TATA solo para los factores foco, porque con TATA hay 668 / 820 promotores.
- **Controles:** TAF1, TBP y GTF2B.

## Resultados

**1. Enriquecimiento en TATA.**

- **GTF2A2:** el 29.4% / 29.8% de sus promotores unidos tiene TATA canónica, contra 6.2% / 6.2% del resto.
- **Odds ratio:** 6.4 (IC 5.0–8.0) / 6.4 (5.2–8.0), Fisher p 4.2e-45 / 8.6e-53. **Es el odds ratio más alto de los 982 factores en ambas réplicas.** Le siguen GPS2 (4.1 / 3.7) y CHAF1B (3.2 / 3.2).
- **Controles:** TAF1 0.70 / 0.69, TBP 0.92 / 0.92 y GTF2B 0.63 / 0.68. Los controles unen, si acaso, menos TATA que el resto.

**2. Atenuación al agregar TATA.** Coeficiente e IC 95%, Rep 1 / Rep 2:

| Factor | Ajustado | + TATA-box | Solo sin TATA | Solo con TATA |
|---|---|---|---|---|
| GTF2A2 | +4.3 (1.4 a 7.1) / +3.7 (1.0 a 6.4) | +3.6 (0.6 a 6.5) / +2.6 (−0.1 a 5.4) | +2.3 (p 0.18) / +0.8 (p 0.62) | +1.2 (p 0.75) / +6.3 (p 0.070) |
| TAF1 | −5.9 / −5.4 | −5.9 / −5.4 | −5.5 / −5.3 | −9.5 / −5.7 |
| TBP | −6.5 / −3.5 | −6.6 / −3.6 | −6.2 / −3.7 | −9.3 / −2.4 |
| GTF2B | −4.7 / −3.6 | −4.6 / −3.5 | −4.4 / −3.4 | −6.9 / −4.0 |

- **GTF2A2 se atenúa un 17% / 29%** con la covariable de TATA.
- **En los promotores sin TATA (n = 293 / 334 unidos) deja de ser significativo:** +2.3 (p 0.18) / +0.8 (p 0.62).
- **Los controles no se mueven** (cambios ≤ 0.1).
- **Dentro de los promotores con TATA** el efecto de GTF2A2 es inconsistente (+1.2 / +6.3), con n chico (122 / 142).
- **Percentil de GTF2A2 en su ventana con la covariable:** 1.00 / 0.92.

**3. Concentración en actividad alta.** Coeficiente y percentil de GTF2A2, Rep 1 / Rep 2:

| Estrato | n unidos | Coeficiente | Percentil |
|---|---|---|---|
| Actividad alta | 199 / 226 | +4.8 (0.5 a 9.1) / +5.8 (1.8 a 9.8) | 0.98 / 0.97 |
| Actividad alta + TATA-box | 199 / 226 | +1.2 (−3.4 a 5.8) / +1.0 (−3.3 a 5.2) | 0.63 / 0.64 |
| Actividad baja | 216 / 250 | +2.2 (p 0.27) / +0.8 (p 0.67) | 0.86 / 0.64 |

- **El efecto de GTF2A2 está en el cuartil más activo,** que es donde la TATA tiene efecto sobre el ruido.
- **Ahí desaparece casi por completo al ajustar por TATA:** de +4.8 / +5.8 a +1.2 / +1.0, y el percentil cae de 0.98 / 0.97 a 0.63 / 0.64.
- **Los controles no cambian:** en actividad alta, TAF1 pasa de −6.8 / −5.6 a −6.2 / −5.0 y TBP de −5.7 / +0.6 a −5.7 / +0.5.
- **Interacción unión × actividad alta:** es positiva para GTF2A2 (+5.2 / +8.3). No la uso como evidencia porque los controles también dan interacciones positivas grandes (TAF1 +5.3 / +8.2, TBP +6.3 / +11.9) aunque en los estratos su efecto es igual o más negativo en actividad alta. En el cuartil alto casi todo está unido por la maquinaria (TAF1: 2162 de 2325), así que el término de interacción con pendientes comunes no es confiable.
- **Isla CpG:** en el estrato con isla, el percentil de GTF2A2 es 1.00 / 0.96. Sin isla no es evaluable (menos de 50 unidos).

**Procedencia.**

- **GTF2A2: una sola línea (K-562),** con 732 promotores en el catálogo y 415 / 476 en el análisis.
- **Experimentos:** ChIP-Atlas lista 6, todos en K-562. El bed no redundante de ReMap no dice cuántos de esos entran en el conjunto fusionado.
- ⚠ **Todo viene de una única línea celular y ninguna es HEK293.**

**SNAPC1.**

- **La premisa no se cumple: la library no tiene promotores de genes de snRNA.** Por `GENETYPE` de org.Hs.eg.db son 0 de 11590 genes; la library es casi toda codificante (11580).
- **Cercanía a genes de snRNA:** solo 4 promotores de la library están a menos de 2 kb de un gen de snRNA (CHRLOC hg38). Tres están unidos por SNAPC1, cerca de RNU12, RNU6-1 y RNU7-1.
- **Excluir esos 3 no cambia nada:** +2.58 / +2.97 contra +2.57 / +2.98.
- **Los promotores que une SNAPC1** son codificantes, con 91% / 90% de isla y un enriquecimiento modesto en TATA (OR 1.5 / 1.5).
- **El efecto no depende de la TATA:** con la covariable es +2.4 / +2.8 y percentil 0.99 / 1.00.
- **Procedencia:** una sola línea (MCF-10A), 1 experimento en ChIP-Atlas.

## Veredicto p8

- **GTF2A2: se cumplen las tres predicciones.**
  1. Sus promotores son los más enriquecidos en TATA canónica de los 982 factores (OR 6.4 / 6.4).
  2. El efecto se atenúa con TATA como covariable y desaparece en los promotores sin TATA. Los controles TAF1, TBP y GTF2B no cambian.
  3. El efecto está en el cuartil más activo y ahí lo explica la TATA (percentil de 0.98 / 0.97 a 0.63 / 0.64).
- **Qué dice y qué no.** GTF2A2 no aporta información sobre el ruido más allá de la TATA: es un marcador de los promotores con TATA. Por eso no es un segundo mecanismo, sino una confirmación desde un dato independiente (ocupación de TFIIA en K-562) de que la asociación TATA–ruido sigue la arquitectura del promotor.
- **Límites:**
  - un solo origen (K-562);
  - n chico (415 / 476);
  - no pasa BH (p perm. 0.039 / 0.082).
- **Recomendación: mención en el texto, no panel.** Encaja como una oración de cierre en la sección de la TATA (`ruido_tata` / `tata_ruido_especificidad`). Una frase posible: «La única subunidad de la maquinaria basal cuya unión se asocia a mayor ruido es GTF2A2 (TFIIA): sus promotores son los más enriquecidos en TATA-box de los 982 factores (OR 6.4 / 6.4), y su efecto desaparece al ajustar por la presencia de TATA en los promotores de alta actividad». Siempre aclarando que viene de una sola línea celular (K-562). Si se quiere un panel suplementario, `p8_efecto_focal_por_modelo` es el candidato, pero no lo recomiendo como panel de R3.7 por el origen único.
- **SNAPC1: descartar.** La library no tiene promotores de snRNA, así que la hipótesis de «clase de promotor PSE» no aplica. El efecto queda como observación aislada de un solo experimento (MCF-10A).

## Archivos p7 / p8

- `p7_grupos_definicion.tsv`: miembros de pausa, preiniciación y mediador, n por réplica, incluidos y motivo de exclusión.
- `p7_p8_residuo_ventana_por_factor_variantes.tsv`: para los 982, coeficiente, residuo loess, ventana y p por permutación en cada variante:
  - ajustado;
  - + ocupación;
  - + TATA;
  - con isla / sin isla;
  - actividad alta / baja, con y sin TATA.
- `p7_miembros_grupos.tsv`, `p7_test_grupos.tsv`: miembros y Wilcoxon de grupo, con HL e IC, en ajustado, + ocupación y los estratos de isla.
- `p7_modelo_conjunto_maquinaria.tsv`: efectos parciales del modelo conjunto, CDK9, CCNT2, NELFA y TAF1 dentro y fuera del núcleo de preiniciación, y solapamiento con CDK9.
- `p7_coef_vs_cobertura.{jpg,pdf}`: coeficiente contra n unidos (log10) con los tres grupos destacados y la curva loess, por réplica.
- `p7_percentil_por_grupo.{jpg,pdf}`: percentil por grupo, ajustado y + ocupación.
- `p7_percentil_ventana_miembros.{jpg,pdf}`: percentil de cada miembro de pausa y preiniciación en ajustado, + ocupación y los estratos de isla.
- `p7_p8_procedencia_lineas_remap.tsv`, `p7_p8_experimentos_chipatlas.tsv`, `p7_p8_efecto_por_linea_celular.tsv`: procedencia y efecto por línea de CDK9, CCNT2, GTF2A2 y SNAPC1.
- `p8_tata_por_factor.tsv`: fracción con TATA en unidos y no unidos, OR, Fisher y rango entre los 982.
- `p8_efecto_focal_por_modelo.{tsv,jpg,pdf}`: coeficiente con IC de GTF2A2, TAF1, TBP, GTF2B y SNAPC1 (y CDK9 y CCNT2 en el .tsv) en cada modelo o subconjunto, con su percentil cuando hay ventana.
- `p8_tata_odds_ratio_vs_cobertura.{jpg,pdf}`: OR de TATA contra n unidos para los 982.
- `p8_snapc1_snrna.tsv`: genes de snRNA en la library y cerca de los unidos por SNAPC1, y efecto sin ellos.
