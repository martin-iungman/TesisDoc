# Panorama de ruido (`ruido_summary`) estratificado por isla CpG

Script: `R/55_exploratorio_ruido_summary_cgi/ruido_summary_cgi.R` (unos segundos, desde la raíz del repo).
Salidas en esta carpeta con prefijo por punto (`p1`…`p7`).
Resultados siempre Rep 1 / Rep 2. AUC = P(ruido con la característica > ruido sin ella). AUC > 0.5 significa más ruido.
No se tocó `docs/tesis.md` ni `docs/mapping_figuras.csv`, y no se promovió nada a figura de tesis.

## Veredicto

**Las asociaciones del panorama no están restringidas a los promotores sin isla CpG. Pero tampoco son generales en el sentido de independientes de la isla.**

- **Ninguna característica difiere entre estratos.**
  - El test de diferencia con isla − sin isla da BH > 0.05 en las 32 características evaluables, en las dos réplicas.
  - La única excepción, en el límite, es «Alta actividad en HEK293», Rep 1 (p corregido 0.045).
  - Por eso no hay base para decir que un efecto «vive en los sin isla». El patrón de `exploratorio_ausencia_senal` (gradiente solo en sin isla) no se generaliza al panorama.
- **Lo que pasa es que buena parte del efecto global es composición entre estratos.**
  - Dentro de cada estrato, casi todos los AUC se acercan a 0.5, en los dos estratos por igual.
  - Las características de contexto endógeno (FANTOM5, especificidad tisular, forma del TSS, elementos transponibles) se concentran en los promotores sin isla. Así heredan el efecto de la isla en el AUC global.
- **El residuo dentro de estratos depende de la réplica.**
  - En la Rep 1 varias características mantienen un efecto chico dentro de ambos estratos (AUC 0.53–0.58).
  - En la Rep 2 desaparece (AUC 0.49–0.51) en los dos estratos.
  - La Rep 2 ya tenía estos efectos más débiles en el global. `figures/exploratorio_replicas_ruido/handoff.md` lo atribuye a la escala de ruido comprimida de la Rep 2, no a peor precisión.
- **Criterio de consistencia:**
  - Global: 15 características lo cumplen.
  - Con isla: 3 (Alto contenido G+C, GC-box, TATA-box).
  - Sin isla: 1 (Sin módulo cis-regulatorio).
  - AUC combinado dentro de estratos («ajustado por isla», ver Método): 5 (G+C, GC-box, CRE, YY1, TATA-box).
- **Para el texto:** las asociaciones con el contexto endógeno no deberían presentarse como asociaciones generales. Lo que se sostiene en las dos réplicas una vez que se compara dentro de la misma clase de isla es:
  - la TATA-box (efecto en los promotores con isla; en los sin isla, potencia insuficiente para decidir);
  - la ausencia de módulo cis-regulatorio (en los sin isla, el único estrato donde se puede evaluar);
  - y efectos chicos de GC-box, G+C, CRE y YY1.

### Características pedidas

AUC [IC 95%], `*` = BH < 0.05 dentro del estrato. «Compatible» = el IC del estrato contiene el AUC global de esa réplica.

**Sin módulo cis-regulatorio**
- Global: 0.628* [0.577, 0.678] / 0.592* [0.549, 0.635].
- Sin isla: 0.599* [0.531, 0.668] / 0.601* [0.548, 0.654]. Consistente y del mismo tamaño que el global.
- Con isla:
  - Rep 1: sin evaluar (46 promotores).
  - Rep 2: 0.526 [0.454, 0.599], n = 54. El IC es compatible tanto con el global como con 0.5.
- **Veredicto:** se sostiene en los sin isla. En los con isla no se puede afirmar ni descartar por falta de potencia: solo el 0.7% de los promotores con isla no tiene CRM, frente al 2.6% de los sin isla. La diferencia entre estratos no es significativa (Rep 2, p corregido 0.43).
- Coincide con `exploratorio_ausencia_senal`.

**No detectado (FANTOM5)**
- Global: 0.597* / 0.543*.
- Con isla: 0.576* [0.535, 0.618] / 0.497 [0.462, 0.533].
- Sin isla: 0.542* [0.515, 0.568] / 0.528 [0.505, 0.552].
- Combinado: 0.552* / 0.519 [0.499, 0.538].
- **Veredicto:** atenuado en los dos estratos. En la Rep 1 el efecto residual está en ambos estratos, y es mayor en los con isla. En la Rep 2 no está en los con isla; el IC excluye el global, así que es ausencia de efecto y no falta de potencia. En los sin isla, Rep 2, es compatible con el global pero no significativo.
- No es un efecto «de los sin isla»: el 73% de los no detectados son sin isla, y de ahí sale gran parte del AUC global.

**Alta especificidad tisular**
- Global: 0.579* / 0.519*.
- Con isla: 0.555* [0.537, 0.573] / 0.494 [0.477, 0.511].
- Sin isla: 0.543* [0.517, 0.569] / 0.506 [0.483, 0.529].
- **Veredicto:** mismo patrón que FANTOM5. Se mantiene en ambos estratos en la Rep 1 y desaparece en ambos en la Rep 2 (con isla: el IC excluye el global).

**Baja especificidad tisular**
- Global: 0.431* / 0.472*.
- Con isla: 0.457* [0.442, 0.471] / 0.490 [0.477, 0.503].
- Sin isla: 0.416* [0.375, 0.457] / 0.456 [0.420, 0.492], con n = 170 / 223.
- Combinado: 0.452* / 0.486 [0.473, 0.498].
- **Veredicto:** en la Rep 2 el IC sin isla excluye 0.5 y es compatible con el global; no llega a BH por potencia. Dentro de los con isla la Rep 2 es casi nula.
- Ojo: solo el 6% de los sin isla tiene baja especificidad, así que es casi una característica de isla.

**Promotores angostos**
- Global: 0.557* / 0.520*.
- Con isla: 0.533* [0.516, 0.550] / 0.510 [0.494, 0.525].
- Sin isla: 0.537* [0.510, 0.563] / 0.505 [0.481, 0.528].
- **Veredicto:** igual a especificidad. Atenuado, presente en ambos estratos en la Rep 1, ausente en la Rep 2.

**Elementos transponibles**
- Global: 0.573* / 0.529*.
- Con isla: 0.539* [0.508, 0.570] / 0.506 [0.478, 0.534].
- Sin isla: 0.541* [0.515, 0.568] / 0.511 [0.488, 0.535].
- **Veredicto:** igual. En la Rep 2 los IC son compatibles con el global (0.529), así que ahí no se puede descartar un efecto chico. En la Rep 1 el global (0.573) cae fuera del IC sin isla: hay atenuación real.

**TSS no canónico**
- Global: 0.551* / 0.520*.
- Con isla: 0.537* [0.517, 0.558] / 0.499 [0.481, 0.518].
- Sin isla: 0.525 [0.501, 0.549] / 0.524 [0.503, 0.545].
- **Veredicto:** en los sin isla, las dos réplicas dan el mismo AUC (~0.525), con IC que excluyen 0.5 en la Rep 2 y apenas tocan 0.5 en la Rep 1. Ninguna llega a BH. Es lo más parecido a «un efecto chico en los sin isla», pero sin potencia para afirmarlo. En los con isla, solo en la Rep 1.

**TATA-box**
- Global: 0.547* / 0.551*.
- Con isla: 0.551* [0.520, 0.583] / 0.561* [0.532, 0.590]. Consistente.
- Sin isla: 0.504 [0.469, 0.540] / 0.519 [0.487, 0.552].
- Combinado: 0.531* / 0.542*.
- **Veredicto:** al revés que la hipótesis de partida, el efecto está en los promotores con isla.
  - En los sin isla, Rep 1, el IC excluye el global.
  - En la Rep 2 es compatible con el global y con 0.5.
  - La diferencia entre estratos no es significativa (p corregido 0.38 / 0.38).
- Reproduce exactamente `exploratorio_tata_ruido/p0_auc_referencia.tsv`. Ahí también se vio que el efecto se concentra en actividad alta, incluso dentro de los sin isla.

## Método

- **Promotores.** Set highconf, 9300 / 11400 tras el descarte de `add_mean_sw_bins()`.
  - Con isla (`CGI_frag`): 6608 / 7900.
  - Sin isla: 2692 / 3500.
  - Ningún promotor sin `CGI_frag`.
- **Características.** `build_tidy_features()` sobre el set completo. Los cortes de tercil (G+C, especificidad, forma, conservación, DNase) y la mediana de HEK293 son los globales, así que cada característica significa lo mismo que en `ruido_summary`. Después se estratifica.
- **AUC.** Igual que `ruido_summary`: `pROC::roc`, IC de DeLong, p de Wald con la varianza de DeLong.
  - BH dentro de cada estrato sobre características × réplicas (63 tests con isla, 65 sin isla), así como el global usa BH sobre características × réplicas (68).
  - El global recalculado coincide exactamente con `ruido_summary_auc.tsv` (el script lo verifica).
- **Mínimo por grupo:** 50 promotores con y 50 sin la característica dentro del estrato y la réplica. Las celdas que no llegan quedan sin evaluar, con el n reportado igual.
- **Diferencia entre estratos.** z = (AUC con isla − AUC sin isla) / √(SE² + SE²), con estratos independientes, BH sobre características × réplicas.
- **AUC combinado dentro de estratos.** Promedio de los dos AUC de estrato ponderado por varianza inversa, BH sobre características × réplicas.
  - Equivale a comparar siempre dentro de la misma clase de isla, con la potencia de los dos estratos juntos.
  - No es parte del pedido original. Lo agregué porque es la forma directa de separar «no significativo por n chico» de «no hay efecto al controlar por isla».
  - Solo se calcula si la característica se evaluó en ambos estratos: «Sin módulo cis-regulatorio» no lo tiene en la Rep 1.
- **Consistencia:** BH < 0.05 y mismo sentido en ambas réplicas, dentro del conjunto (global, estrato o combinado).

## Punto 2 y punto 6. n, celdas sin evaluar y características atadas a la isla

Archivos: `p2_sin_evaluar.tsv`, `p6_prevalencia_por_estrato.tsv`, y n con / sin de cada celda en `p1_p2_auc_por_estrato.tsv`.

Sin evaluar:

| Característica | Estrato | n con (Rep 1 / Rep 2) | Motivo |
|---|---|---|---|
| Islas CpG | ambos | 6608 / 7900 y 0 / 0 | define el estrato |
| Sin módulo cis-regulatorio | con isla, solo Rep 1 | 46 / 54 | n < 50 (en la Rep 2 se evalúa con 54) |
| Tracto poli(dA:dT) proximal | con isla, ambas | 34 / 49 | n < 50 |
| Tracto poli(dA:dT) proximal | sin isla, Rep 1 | 49 / 82 | n < 50 |

- El tracto poli(dA:dT) no tiene que ver con la isla; es raro en general (83 / 131). Ya era n.s. en el global.
- **Características ligadas a la definición de isla:**

| Característica | Prevalencia con isla | Prevalencia sin isla | Decisión |
|---|---|---|---|
| Islas CpG | 100% | 0% | Excluida: sin variación por definición |
| Alto contenido G+C | 44.7% / 44.7% | 5.2% / 4.1% (139 / 142) | Evaluada, con cautela |
| CG en TSS | 20.2% / 20.1% | 5.5% / 5.3% (149 / 187) | Evaluada, con cautela |
| GC-box | 54.3% / 53.7% | 40.8% / 38.9% | Evaluada: no está determinada por la isla |

- **Por qué G+C y CG en TSS quedan.**
  - Pasan el mínimo de n en ambos estratos.
  - Dentro de los con isla miden un gradiente de composición residual, no la isla misma. Ahí G+C alto es consistente (0.467* / 0.453*).
  - En los sin isla son el 4–5% de promotores en el borde de la definición. Su IC es ancho y compatible con el global (G+C 0.451 [0.400, 0.502] / 0.454 [0.407, 0.500]): ahí no hay potencia para decir nada.
- **Casi determinadas por la isla sin ser parte de su definición.** En los sin isla quedan con prevalencia del 5–7%: baja especificidad tisular (170 / 223), alta actividad en HEK293 (195 / 240) y alta accesibilidad DNase (150 / 188). Se evalúan, pero con IC anchos.

## Punto 3. Tabla comparativa

AUC Rep 1 / Rep 2, `*` = BH < 0.05 en su conjunto. «Combinado» = dentro de estratos ponderado. «n» = promotores con la característica (los n sin están en el TSV).
Estado: alto / bajo = consistente con ese sentido; no = no consistente; s/e = sin evaluar.
IC 95%, diferencias entre estratos y compatibilidad con el global: `p3_comparativa_global_estratos.tsv`.

### Elementos del promotor (seq)

| Característica | Global | Con isla | Sin isla | Combinado | n con isla | n sin isla | Estado G / C / S / Comb |
|---|---|---|---|---|---|---|---|
| Islas CpG | 0.399* / 0.444* | — / — | — / — | — | 6608 / 7900 | 0 / 0 | bajo / s/e / s/e / — |
| Alto contenido G+C | 0.435* / 0.440* | 0.467* / 0.453* | 0.451 / 0.454 | 0.465* / 0.453* | 2952 / 3528 | 139 / 142 | bajo / bajo / no / bajo |
| CG en TSS | 0.445* / 0.482* | 0.467* / 0.497 | 0.448 / 0.481 | 0.464* / 0.495 | 1335 / 1588 | 149 / 187 | bajo / no / no / no |
| GC-box | 0.464* / 0.476* | 0.472* / 0.474* | 0.482 / 0.500 | 0.475* / 0.482* | 3590 / 4245 | 1097 / 1362 | bajo / bajo / no / bajo |
| CCAAT | 0.487 / 0.501 | 0.491 / 0.505 | 0.483 / 0.493 | 0.489 / 0.502 | 1275 / 1451 | 479 / 594 | no / no / no / no |
| TG en TSS | 0.491 / 0.493 | 0.496 / 0.494 | 0.486 / 0.494 | 0.494 / 0.494 | 915 / 1070 | 335 / 418 | no / no / no / no |
| CA en TSS | 0.496 / 0.497 | 0.500 / 0.502 | 0.483 / 0.482 | 0.495 / 0.495 | 2323 / 2805 | 1022 / 1356 | no / no / no / no |
| TCT | 0.501 / 0.503 | 0.505 / 0.502 | 0.513 / 0.522 | 0.507 / 0.507 | 1285 / 1484 | 400 / 475 | no / no / no / no |
| Tracto poli(dA:dT) proximal | 0.509 / 0.513 | — / — | — / 0.492 | — | 34 / 49 | 49 / 82 | no / s/e / s/e / — |
| TSS fuerte | 0.514 / 0.505 | 0.516 / 0.510 | 0.492 / 0.486 | 0.508 / 0.502 | 889 / 1048 | 466 / 609 | no / no / no / no |
| TA en TSS | 0.520 / 0.511 | 0.517 / 0.516 | 0.537 / 0.503 | 0.522 / 0.512 | 530 / 613 | 207 / 276 | no / no / no / no |
| CRE (CREB/ATF) | 0.531* / 0.529* | 0.519 / 0.529* | 0.526 / 0.514 | 0.521* / 0.524* | 821 / 944 | 503 / 676 | alto / no / no / alto |
| YY1 | 0.532* / 0.526* | 0.521 / 0.521 | 0.520 / 0.517 | 0.521* / 0.519* | 1023 / 1207 | 640 / 835 | alto / no / no / alto |
| E-box CAGGTR | 0.536* / 0.499 | 0.514 / 0.483 | 0.546* / 0.510 | 0.525* / 0.493 | 1286 / 1551 | 763 / 1003 | no / no / no / no |
| TATA-box | 0.547* / 0.551* | 0.551* / 0.561* | 0.504 / 0.519 | 0.531* / 0.542* | 396 / 469 | 272 / 351 | alto / alto / no / alto |
| SINE | 0.548* / 0.536* | 0.527 / 0.519 | 0.522 / 0.525 | 0.525 / 0.523 | 224 / 271 | 250 / 327 | alto / no / no / no |
| TSS no canónico | 0.551* / 0.520* | 0.537* / 0.499 | 0.525 / 0.524 | 0.532* / 0.510 | 925 / 1165 | 742 / 994 | alto / no / no / no |
| LINE | 0.560* / 0.495 | 0.522 / 0.464 | 0.522 / 0.478 | 0.522 / 0.474 | 71 / 90 | 146 / 196 | no / no / no / no |
| Elementos transponibles | 0.573* / 0.529* | 0.539* / 0.506 | 0.541* / 0.511 | 0.540* / 0.509 | 360 / 436 | 569 / 775 | alto / no / no / no |

### Características del promotor (endo)

| Característica | Global | Con isla | Sin isla | Combinado | n con isla | n sin isla | Estado G / C / S / Comb |
|---|---|---|---|---|---|---|---|
| Alta accesibilidad de cromatina (DNase-seq) | 0.430* / 0.491 | 0.467* / 0.518* | 0.401* / 0.480 | 0.460* / 0.514 | 2950 / 3612 | 150 / 188 | no / no / no / no |
| Baja especificidad tisular | 0.431* / 0.472* | 0.457* / 0.490 | 0.416* / 0.456 | 0.452* / 0.486 | 2669 / 3238 | 170 / 223 | bajo / no / no / no |
| Alta actividad en HEK293 | 0.445* / 0.501 | 0.480* / 0.527* | 0.411* / 0.470 | 0.472* / 0.520* | 2517 / 2857 | 195 / 240 | no / no / no / no |
| Promotores anchos | 0.470* / 0.503 | 0.490 / 0.516 | 0.466 / 0.492 | 0.487 / 0.512 | 2357 / 2747 | 271 / 362 | no / no / no / no |
| Enhancers a 50kb | 0.497 / 0.510 | 0.492 / 0.502 | 0.485 / 0.516 | 0.490 / 0.506 | 4474 / 5378 | 2018 / 2605 | no / no / no / no |
| Promotor unidireccional | 0.500 / 0.519* | 0.490 / 0.527* | 0.521 / 0.478 | 0.495 / 0.517 | 1146 / 1345 | 276 / 375 | no / no / no / no |
| Alta conservación (16 a -50pb) | 0.510 / 0.512 | 0.509 / 0.514 | 0.505 / 0.501 | 0.508 / 0.510 | 2138 / 2548 | 961 / 1251 | no / no / no / no |
| Alta conservación (-50 a -150) | 0.512 / 0.513* | 0.509 / 0.510 | 0.499 / 0.507 | 0.506 / 0.509 | 2059 / 2451 | 1040 / 1348 | no / no / no / no |
| Insertado en humanos | 0.513 / 0.520 | 0.504 / 0.510 | 0.511 / 0.532 | 0.506 / 0.516 | 132 / 148 | 50 / 61 | no / no / no / no |
| Alta conservación (-150 a -235) | 0.515* / 0.496 | 0.515 / 0.491 | 0.498 / 0.499 | 0.510 / 0.493 | 2093 / 2502 | 1006 / 1297 | no / no / no / no |
| Sin actividad en ratón | 0.529* / 0.504 | 0.524 / 0.489 | 0.515 / 0.512 | 0.521 / 0.498 | 354 / 431 | 218 / 309 | no / no / no / no |
| Promotores angostos | 0.557* / 0.520* | 0.533* / 0.510 | 0.537* / 0.505 | 0.534* / 0.508 | 1565 / 1817 | 1120 / 1457 | alto / no / no / no |
| Alta especificidad tisular | 0.579* / 0.519* | 0.555* / 0.494 | 0.543* / 0.506 | 0.551* / 0.498 | 1380 / 1611 | 1459 / 1850 | alto / no / no / no |
| No detectado (FANTOM5) | 0.597* / 0.543* | 0.576* / 0.497 | 0.542* / 0.528 | 0.552* / 0.519 | 209 / 271 | 571 / 741 | alto / no / no / no |
| Sin módulo cis-regulatorio | 0.628* / 0.592* | — / 0.526 | 0.599* / 0.601* | — / 0.575* | 46 / 54 | 69 / 90 | alto / s/e / alto / no |

## Punto 4. Consistencia y cambios de estado

Archivos: `p4_consistencia_resumen.tsv`, `p4_consistencia_estado.tsv`.

| Conjunto | Evaluadas | Consistentes | Ruido alto | Ruido bajo |
|---|---|---|---|---|
| Global | 34 | 15 | 10 | 5 |
| Con isla | 31 | 3 | 1 (TATA-box) | 2 (G+C, GC-box) |
| Sin isla | 32 | 1 | 1 (Sin módulo cis-regulatorio) | 0 |
| Combinado dentro de estratos | 32 (Sin módulo cis-regulatorio solo en Rep 2) | 5 | 3 (TATA-box, CRE, YY1) | 2 (G+C, GC-box) |

- **Las 15 consistentes globales**, y en qué estrato se mantienen:
  - Se mantienen en con isla: G+C, GC-box, TATA-box.
  - Se mantiene en sin isla: Sin módulo cis-regulatorio.
  - No se mantienen en ningún estrato: CG en TSS, CRE, YY1, SINE, TSS no canónico, Elementos transponibles, No detectado (FANTOM5), Alta y Baja especificidad tisular, Promotores angostos.
  - Islas CpG: sin evaluar por definición.
- **Ninguna característica gana consistencia en un estrato sin tenerla en el global.**
- **Ninguna cambia de sentido de forma consistente.** Lo que sí aparece es discordancia entre réplicas dentro de los con isla: ambas significativas, pero de sentido opuesto.
  - Alta actividad en HEK293: 0.480* / 0.527*.
  - Alta accesibilidad DNase: 0.467* / 0.518*.
  - En la Rep 2, dentro de los con isla, estas dos van hacia más ruido; en el global y en la Rep 1, hacia menos. Además, «Promotor unidireccional» da 0.490 / 0.527*.
  - Ninguna de las tres era consistente en el global, así que no cambia el panorama, pero conviene no apoyarse en ellas.

## Punto 5. Rango global restringido vs. rango recalculado dentro del estrato

Archivos: `p5_sensibilidad_rango_por_estrato.tsv`, `p5_sensibilidad_resumen.tsv`.

| Estrato | Celdas | r de Pearson (AUC) | Mediana \|ΔAUC\| | Máx. \|ΔAUC\| | Celdas que cambian de significancia | Características que cambian de estado |
|---|---|---|---|---|---|---|
| Con isla | 63 | 0.983 | 0.002 | 0.016 | 3 | 0 |
| Sin isla | 65 | 0.983 | 0.004 | 0.017 | 1 | 0 |

- **Las celdas que cambian están todas en el borde de BH:**
  - Con isla: Elementos transponibles Rep 1 (p corregido 0.048 → 0.059), DNase Rep 2 (0.018 → 0.100), Promotor unidireccional Rep 2 (0.008 → 0.059).
  - Sin isla: TSS no canónico Rep 1 (0.196 → 0.039).
- **Igual que con la TATA, las dos versiones dan lo mismo:** ninguna característica cambia de estado de consistencia.

## Punto 7. Figura

`p7_summary_noise_cgi.jpg` / `.pdf`, con todos los valores en `p7_summary_noise_cgi.tsv`.

- **Filas:** las características consistentes en al menos uno de los tres conjuntos, separadas en A = elementos del promotor y B = características del promotor, como en `ruido_summary`. Orden por AUC global medio.
- **Columnas:** Global / Con isla CpG / Sin isla CpG, con el mismo eje x en las tres.
- **Colores:**
  - Barras naranja / verde si la característica es consistente en ese conjunto, gris si no.
  - Rep 1 opaca y Rep 2 más clara, igual que en `ruido_summary`.
- **n:** a la derecha de cada columna, el n con la característica (Rep 1 / Rep 2). «n < 50 (46 / 54)» marca la celda no evaluada; en ese caso igual se dibuja la réplica que sí se evaluó.
- **Mejoras posibles si se promueve:**
  - Agregar la columna combinada.
  - Marcar con un símbolo las celdas cuyo IC es compatible con el global.
