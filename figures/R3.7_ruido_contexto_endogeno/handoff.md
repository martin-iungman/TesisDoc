# ruido_contexto_endogeno: fusión de remap_tf_noise_gsea + chipatlas_histonas + gradiente de ocupación

Script: `R/56_ruido_contexto_endogeno/ruido_contexto_endogeno.R` (~2 min, desde la raíz del repo; salida vía `fig_dir("ruido_contexto_endogeno")`).
Resultados siempre Rep 1 / Rep 2. AUC = P(ruido con la marca o factor > ruido sin ella); AUC < 0.5 = menos ruido.
No se tocó `docs/tesis.md`.

## Resumen

- **Hay dos discrepancias con los valores que me pasaste.**
  - **AUC de H3K4me3 en la Rep 1: da 0.369, no 0.38.** La tabla vieja de `chipatlas_histonas` estaba desactualizada (ver punto 5c).
  - **Orden después de H3K4me3: en la Rep 2 no es H3K27ac y H3K9ac.** H3K4me1, H3K9ac y H3K27ac quedan prácticamente empatados (ver punto 5c).
- **El resto coincide:** 982 / 514, NES y p ajustado de MLL1/2, 63.1 / 59.6 con 115 / 144 promotores frente a 50.5, y la meseta. Hay matices en el p ajustado 0.081 de la Rep 2 y en la meseta de la Rep 2 (ver punto 5).
- **La asociación de H3K4me3 no se puede evaluar en los promotores con isla.** El 99.3% de ellos tiene la marca: solo 43 / 56 no la tienen.
  - Así que no es que «se concentre en los sin isla»: **solo en los sin isla hay variación para medirla.**
  - Ahí es significativa en ambas réplicas: 0.422 / 0.468.
- **Mapping:** la figura fusionada pasa a ser R3.7. `remap_tf_noise_gsea` y `chipatlas_histonas` quedan como `fusionada`, sin número. **Ninguna otra figura se renumera:** `chipatlas_histonas` (R3.8) era la última de Resultados 3 y ya no estaba en `docs/Fig R3.pptx`, que tiene 7 slides.
- **Hay que actualizar el deck de Drive.** El slide 7 de `Fig R3` todavía tiene la figura vieja de ReMap: reemplazalo por `ruido_contexto_endogeno.jpg`. Si no, la rutina diaria `resync-figuras-mapping` va a ver ese slide como `remap_tf_noise_gsea`.

## Qué se hizo

- **Paneles.** Cada uno se exporta suelto (`panel_a` a `panel_e`) y en la compuesta `ruido_contexto_endogeno`, todo en JPG y PDF.
  - **A, B y C van en el formato de las figuras originales**, a pedido tuyo: barras de AUC − 0.5 de `plot_noise_auc_summary()`, con Rep 1 opaca y Rep 2 más clara.
    - Las réplicas van por transparencia en A/B/C y en facetas en D/E. Es la única diferencia de codificación entre paneles.
    - Usé el mismo tamaño de fuente (base 16) en los cinco paneles.
  - **A** (`panel_a_remap_auc`): los 514 factores consistentes, todos de ruido bajo, sin etiquetas ni IC, como `Remap_noise`. El texto que acompaña la figura dice que se evaluaron 982.
  - **B** (`panel_b_gsea_mll`): como `GSEA_MLL1_2`.
    - Curva de enriquecimiento de la Rep 1 y, debajo, todos los factores en el mismo orden, con los 19 del complejo resaltados.
    - El NES y el p ajustado de las dos réplicas van en el texto o en la leyenda (ver punto 5).
  - **C** (`panel_c_chipatlas_auc`): como `chipatlas_noise`, con las 23 marcas consistentes de las 49 evaluadas.
    - Barras con IC 95%. H3K4me3 resaltada en color, tanto en la etiqueta como en su prevalencia.
    - Columna adjunta con el % de promotores con la marca (Rep 1 / Rep 2).
    - El orden es el de `plot_noise_auc_summary()`, por AUC de la Rep 2. El panel C ya no muestra las 26 marcas n.s. (están en el TSV).
    - Cambio respecto del original: la significancia es p de Wald DeLong con BH, como `ruido_summary`; el script viejo usaba el criterio de IC. Por eso la cantidad de marcas que se muestran puede diferir del panel viejo.
  - **D** (`panel_d_h3k4me3_sin_isla`): scatter sin isla CpG, con bins y rango recalculados dentro del subconjunto, como el original.
    - Loess con IC 95%.
    - En cada faceta: n con y sin H3K4me3, y el AUC.
    - Leyenda «Con H3K4me3» / «Sin H3K4me3».
  - **E** (`panel_e_ocupacion`): bin 0 + deciles del resto (por réplica), por isla CpG.
    - IC 95% t, n debajo de cada punto y ρ de Spearman por estrato.
    - Los rangos de factores de cada decil están en `panel_e_ocupacion.tsv`.
    - **Diferencia con el exploratorio:** acá va sobre el set completo. El exploratorio excluía 4 / 5 promotores sin clase de FANTOM5 porque usaba FANTOM5. Los valores no cambian a la precisión reportada.
- **Estilo.**
  - Réplicas por transparencia en A/B/C (formato original) y en facetas en D/E.
  - `theme_pubr` con base 16 en todos.
  - Figura compuesta: A | B, C | D, E abajo a lo ancho, en orden de lectura.
  - Ejes con punto decimal y sin separador de miles.
  - Ninguna leyenda dice TRUE/FALSE.
  - Límites de eje calculados para incluir los IC completos.
- **Colores.** Orange `#D6741F` es «sin isla» en E y «ruido alto» en A/C. Ningún factor ni marca da ruido alto, así que esa categoría no aparece en la figura y no hay ambigüedad. Si preferís, lo cambio.
- **Tablas:**
  - `numeros_texto.tsv`: los números de abajo.
  - `panel_a_remap_auc.tsv`, `panel_b_gsea_mll.tsv`, `panel_c_chipatlas_auc.tsv`, `panel_e_ocupacion.tsv`.
  - `h3k4me3_auc_por_isla.tsv`, `meseta_tramos_fijos.tsv`, `meseta_quiebre.tsv`.

## Números pedidos

**1. Marcas de histona evaluadas: 49.** Criterio: >100 promotores con la marca en ambas réplicas, como el original. De esas, 23 tienen efecto consistente, todas de ruido bajo; 26 son n.s.

**2. Prevalencia de H3K4me3 en la library:** 90.7% (8438 de 9300) / 90.3% (10299 de 11400).

**3. H3K4me3 dentro de los promotores sin isla CpG.** Con bins y rango recalculados en el estrato, que es lo que muestra el panel D:

| | AUC [IC 95%] | p | Con H3K4me3 | Sin H3K4me3 |
|---|---|---|---|---|
| Rep 1 | 0.422 [0.398, 0.446] | 2.2e-10 | 1808 | 792 |
| Rep 2 | 0.468 [0.447, 0.489] | 0.0027 | 2455 | 1045 |

- **Los n del panel D (2600 / 3500) son menores que los del estrato completo (2692 / 3500).** Al rebinear en bins de 100, `add_mean_sw_bins()` descarta el resto de la división.
- **Sensibilidad:** con el rango global restringido al estrato y el estrato completo, da 0.425 / 0.463 (n 1873 con / 819 sin, y 2455 / 1045). Es lo mismo.
- En los sin isla, el 69.5% / 70.1% tiene H3K4me3.

**4. H3K4me3 dentro de los promotores con isla CpG:**

| | AUC [IC 95%] | p | Con H3K4me3 | Sin H3K4me3 |
|---|---|---|---|---|
| Rep 1 | 0.418 [0.327, 0.509] | 0.076 | 6557 | 43 |
| Rep 2 | 0.526 [0.450, 0.602] | 0.51 | 7844 | 56 |

- **No hay comparación posible.** El 99.3% de los promotores con isla tiene H3K4me3, y los 43 / 56 sin la marca dan IC de ±0.09 y estimaciones de sentido opuesto entre réplicas.
- **Por eso no corresponde decir que la asociación «se concentra en los sin isla».** Lo correcto es que se evalúa en los sin isla, porque solo ahí hay variación.
- **Redacción sugerida:** «la ausencia de H3K4me3 es casi exclusiva de los promotores sin isla CpG; dentro de ese grupo, los promotores sin la marca son más ruidosos (AUC 0.42 / 0.47)».
- El AUC global (0.369 / 0.435) mezcla ese efecto con el de la isla. Casi todos los promotores sin H3K4me3 son sin isla, que ya son más ruidosos.

## 5. Confirmación de los valores del texto

| Valor en el texto | Recalculado | ¿Coincide? |
|---|---|---|
| 982 factores evaluados | 982 | ✓ |
| 514 con efecto consistente, ninguno de ruido alto | 514 ruido bajo, 0 ruido alto, 468 n.s. | ✓ |
| NES MLL1/2 −2.47 / −2.01 | −2.47 / −2.01 | ✓ |
| p ajustado 0.005 / 0.081 | 0.005 / 0.081 | ✓ (ver nota a) |
| AUC H3K4me3 0.38 / 0.44 | **0.369 / 0.435** | **✗ en Rep 1** (ver nota c) |
| seguido por H3K27ac y H3K9ac | Rep 1: H3K27ac 0.372, H3K9ac 0.394. **Rep 2: H3K4me1 0.443, H3K9ac 0.445, H3K27ac 0.448** | **Solo en Rep 1** (ver nota c) |
| 63.1 / 59.6 en 115 / 144 promotores sin ningún pico de ReMap | 63.1 [58.0, 68.2] / 59.6 [55.3, 63.9], n 115 / 144 | ✓ |
| frente a 50.5 en el conjunto | 50.5 / 50.5 | ✓ (es 50.5 por construcción: rango 1–100) |
| meseta entre 300 y 400 factores | quiebre 400 / 375 | ✓ en Rep 1; en Rep 2 con matiz (ver nota e) |

**a. p ajustado 0.081 en la Rep 2.**
- Es significativo solo con el corte 0.1 que usa el GSEA (`pvalueCutoff = 0.1`), no con 0.05.
- Si el texto dice «significativo en ambas réplicas», conviene aclarar el umbral, o decir «0.005 en la Rep 1 y 0.081 en la Rep 2».
- El conjunto tiene 19 factores en ambas réplicas.

**c. H3K4me3 y orden de las marcas.**
- **De dónde sale la diferencia.** El 0.38 de la Rep 1 viene de `chipatlas_noise_auc.tsv` (0.378). Ese archivo es del 29/07, anterior al `activity_stats_highconf.tsv` actual (22/10), y no se había regenerado.
- **Lo verifiqué** corriendo `R/22_chipatlas_histonas/chipatlas_histonas.R` sin cambios sobre los datos actuales, con salida al scratchpad: da 0.369 / 0.435, igual que el script nuevo.
- **Valor correcto: 0.37 / 0.44.**
- **Orden.**
  - Rep 1: H3K4me3 (0.369), H3K27ac (0.372), H2BK15ac (0.383), H3K9ac (0.394), H3K4me2 (0.399).
  - Rep 2: H3K4me3 (0.435), H3K4me1 (0.443), H3K9ac (0.445), H3K27ac (0.448), H3K4me2 (0.451).
  - H2BK15ac es n.s. (IC ancho) y solo el 1.4% de los promotores la tiene.
  - Por AUC medio entre réplicas, que es el orden del panel C: H3K4me3 (0.402), H3K27ac (0.410), H3K9ac (0.420), H3K4me1 (0.423), H3K4me2 (0.425).
- **Sugerencia.** «Seguido por H3K27ac y H3K9ac» vale por AUC medio y en la Rep 1. En la Rep 2, las marcas de segundo a quinto lugar están dentro de 0.01 de AUC.
- **Prevalencias (Rep 1 / Rep 2):** H3K4me3 90.7 / 90.3%, H3K27ac 93.4 / 92.9%, H3K9ac 83.0 / 81.8%, H3K4me1 86.3 / 85.3%.
  - Las marcas con AUC más extremo son las más prevalentes. Ese es el motivo de la columna de prevalencia del panel C: la AUC de una marca presente en el 90% compara ese 90% contra un 10% muy particular.

**e. Meseta.**

Medias de `var_rank_sw` por tramo fijo de factores (`meseta_tramos_fijos.tsv`):

| Factores | 0 | 1–100 | 101–200 | 201–300 | 301–400 | 401–500 | 501–600 | >600 |
|---|---|---|---|---|---|---|---|---|
| Rep 1 | 63.1 | 59.4 | 56.2 | 51.5 | 48.7 | 47.1 | 47.4 | 46.6 |
| Rep 2 | 59.6 | 54.2 | 51.9 | 49.8 | 48.3 | 49.5 | 49.2 | 51.1 |

- **Recta quebrada** (quiebre por grilla de 50 a 600, sin bootstrap):
  - Rep 1: quiebre en 400. Pendiente −3.4 antes y ≈0 después (por cada 100 factores).
  - Rep 2: quiebre en 375. Pendiente −2.1 antes y +1.2 después.
  - Los IC bootstrap del exploratorio (`exploratorio_ausencia_senal/p2d`) son 250–500 / 50–400.
- **Lectura.**
  - Rep 1: meseta clara desde ~300–400.
  - Rep 2: mínimo en 300–400 y una suba leve después (51.1 en >600).
- **Sugerencia.** «La reducción del ruido satura alrededor de 300–400 factores» vale para ambas réplicas. «Meseta» es exacto para la Rep 1; para la Rep 2 sería «sin reducción adicional».

**Para el texto, sobre el panel E.**
- En los sin isla el gradiente es claro: ρ −0.13 / −0.08, p 8e-12 / 1e-6.
- En los con isla es mínimo y de signo opuesto entre réplicas: ρ −0.03 / +0.04, p 0.005 / 1e-4.
- Por la n grande, el +0.04 de la Rep 2 es significativo: en los con isla de la Rep 2, más factores van con apenas más ruido.
- No contradice `exploratorio_ausencia_senal`. El ρ −0.003 de ese handoff es del conteo de ausencias, no de la ocupación.

## Cambios en el repo

- **`docs/mapping_figuras.csv`:**
  - **Nueva fila `ruido_contexto_endogeno`:** Fig. R3.7, `listo`, `leyenda_lista = no` (la actualiza la rutina `sync-tesis-md`). Archivos y notas completos.
  - **`remap_tf_noise_gsea` y `chipatlas_histonas`:**
    - `numero_actual = -` y `estado = fusionada`, con la nota «Fusionada en ruido_contexto_endogeno (paneles …)».
    - A `chipatlas_histonas` se le quitó el «Revisar:» (preguntaba si la figura había salido del deck: ya está resuelto) y se le agregó la advertencia de la tabla desactualizada.
    - **Elegí marcarlas en vez de borrarlas:** sus scripts siguen existiendo y otros exploratorios (`R/51`, `R/52`) leen sus tablas vía `fig_dir()`. Con la fila, siguen resolviendo a `figures/sinnum_<slug>/`.
    - `fusionada` es un valor nuevo de `estado`. No había precedente en el CSV.
- **`R/maintenance/sync_figure_folders.R`:**
  - `R3.7_remap_tf_noise_gsea` → `sinnum_remap_tf_noise_gsea`.
  - `R3.8_chipatlas_histonas` → `sinnum_chipatlas_histonas`.
  - `figures/R3.7_ruido_contexto_endogeno/` creada por el script.
  - No se renumeró ninguna otra figura.
- **Pendiente tuyo:** reemplazar el slide 7 de `Fig R3` en Drive por `ruido_contexto_endogeno.jpg` (ver Resumen).
- **Nada commiteado.**
