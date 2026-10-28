# Handoff — TATA-box y ruido transcripcional (2026-10-28)

Exploratorio para la sección TATA-box de R3 (ruido). Pregunta central: la
TATA-box se asocia a mayor ruido solo en promotores de alta actividad; ¿es
porque el efecto depende del nivel de actividad (**lectura A**) o porque las
TATA de los promotores poco activos son más débiles o están mal posicionadas
(**lectura B**)?

- **Script:** `R/44_exploratorio_tata_ruido/tata_ruido.R` (corre entero en
  unos minutos desde la raíz del repo; recalcula la búsqueda del motivo, que
  tarda segundos). Puntos p0–p11, cada uno con prefijo en los archivos.
- **Salidas:** `figures/exploratorio_tata_ruido/` (sin `fig_dir()`, sin número).
- **No es figura de la tesis:** no hay fila en `docs/mapping_figuras.csv` y no
  hay que agregarla.
- Commits: `d7b4d1d` (p0–p5), `26b0c09` (p6), `ad2a3df` (p7), `c70c5bb` (p8), `17440d3` (p9), `147e5d8` (p10), p11 en el commit que agrega esa sección.

Números siempre como Rep 1 / Rep 2.

---

## Método (Paso 0 — importante para comparar con R3)

- **Rango de ruido (`var_rank_sw`):** `add_noise_rank()` sobre todos los
  promotores high-conf ∩ `prom_df` de cada réplica, igual que ruido_cgi:
  intervalos de 100 promotores por `mean` (`row_number`, se descartan los
  `n %% 100` menos activos → 9300 / 11400) y `row_number(var)` dentro de cada
  intervalo (empates rotos por orden de fila; 35 / 416 promotores). Después
  se restringe a los FP sin N (los que tienen búsqueda del motivo).
- **AUC:** exacta (Mann-Whitney, pROC) con IC 95% y p de Wald de DeLong
  (`auc_vs()` en el script). AUC > 0.5 = la feature se asocia a más ruido.
- **Actividad alta:** percentil 75 de `mean` dentro de cada réplica.
- **TATA:** `R/functions/tata_scan.R` (Bucher + pesos W de EPD, corte −8.17,
  canónica −31 a −26 misma hebra = `TATA_EPD` exacto; ventanas de 10 pb
  alineadas en −33 a −24; 10 matrices permutadas, seed 2024).
- **"Sin motivo":** sin ninguna coincidencia ≥ −8.17 en el fragmento, en
  ninguna hebra.
- **TATA fuerte (p4 en adelante):** canónica con puntaje ≥ −5.67 (mediana de
  las canónicas; −5.68 / −5.67 por réplica). El umbral −5 da lo mismo.

### Problemas encontrados en las figuras de R3 (ya corregidos)
1. `noise_auc()` (trapecio sobre umbrales 1–99) subestimaba ~0.01 y fallaba
   en estratos chicos. Ahora `noise_roc()` y `noise_roc_by_expr()` usan
   `noise_auc_exact()` (pROC) — `R/functions/plot_helpers.R`.
2. `noise_roc_by_expr()` cortaba en `0.75 * max(mean_sw)` con las dos réplicas
   juntas (en Rep 1 "alta" era el 8.6% superior). Ahora percentil 75 de `mean`
   por réplica.
3. `ruido_summary_auc.tsv` / `summary_noise.jpg` estaban desactualizados.

Regenerados: ruido_cgi (isla CpG 0.399 / 0.444, antes 0.390 / 0.434),
ruido_tata (TATA 0.547 / 0.551; ROC_TATA_expr alta 0.595 / 0.606, baja
0.515 / 0.519) y ruido_summary (15 features consistentes; entran YY1 y CRE).
Filas del mapping actualizadas. Commits `20a61a7`, `e09c267`, `35222f7`.

---

## Resultados

### p0 — valores de referencia (TATA_EPD vs. el resto)
- Global 0.547 / 0.551; alta actividad 0.595 / 0.606; baja 0.515 / 0.519.
- Con isla 0.551 / 0.561; sin isla 0.504 / 0.519.
- `noise_auc()` global, para comparar con versiones viejas: isla CpG
  0.390 / 0.434, TATA 0.538 / 0.541.

### p1 — posición
- Ventana canónica (−33 a −24, misma hebra): AUC 0.554 / 0.554 (n 689 / 850),
  BH < 0.05.
- El pico **no es exclusivo** (a diferencia de la actividad): el perfil queda
  en ~0.52–0.56 en casi todo el fragmento.
- Por clase vs. sin motivo: canónica 0.556 / 0.556; solo no canónica
  0.525 / 0.510; solo hebra opuesta 0.532 / 0.529; hebra opuesta en la
  ventana canónica sin la canónica 0.559 / 0.527 (n 97 / 122).

### p2 — composición
- Permutadas: algunas reproducen AUC > 0.5 en la ventana canónica (perm4
  0.583 / 0.606, perm6 0.594 / 0.548).
- Modelo `var_rank_sw ~ clase + A/T proximal + A/T del fragmento`: la
  canónica baja de +5.7 / +5.7 rangos a +2.7 / +3.4 pero se mantiene
  (p 0.024 / 0.002); las no canónicas desaparecen.
- **≈ la mitad del efecto sobre el ruido es composición A/T**; queda un
  efecto propio del motivo canónico.

### p3 — puntaje
- Plano (~0.50–0.52) hasta −7; salto en ≥ −5 (0.572 / 0.571). Bins
  intermedios inconsistentes entre réplicas. Ventana control ~0.5.

### p4 — puntaje × actividad (A vs. B) — **resultado clave**

| | Baja actividad | Alta actividad |
|---|---|---|
| Puntaje alto | 0.520 / 0.508 (n 172 / 234) | **0.648 / 0.642** (n 162 / 177) |
| Puntaje bajo | 0.529 / 0.542 (n 213 / 271) | 0.538 / 0.562 (n 121 / 138) |

- Alto vs. bajo: p 0.004 / 0.023 en alta actividad; n.s. en baja.
- Spearman puntaje–ruido en las canónicas: 0.16 / 0.17 en alta, ~0 en baja.
- Las TATA de baja actividad son apenas más débiles (mediana −5.9 vs. −5.2).
- Potencia: IC superior de "alto × baja" = 0.56 / 0.55, excluye ~0.64.
- **Apoya la lectura A.**

### p5 — dentro de islas CpG
- Mismo patrón: puntaje alto, alta actividad 0.626 / 0.647; baja actividad
  0.498 / 0.480. Sin isla: n 28–49, no interpretable.

### p6 — qué distingue a las TATA fuertes de baja actividad
TATA fuertes alta vs. baja (Fisher, BH) + control `variable ~ actividad_alta *
tata_fuerte` sobre todos los promotores. 88% de las TATA fuertes caen en el
mismo grupo de actividad en ambas réplicas (133 alta / 145 baja); con las
concordantes da lo mismo.

- **Específico de las TATA (interacción):** alta especificidad tisular
  (82 / 81% en baja vs. 45 / 38% en alta; OR 0.18 / 0.14 vs. 0.40 / 0.41 en
  el resto), alta actividad en HEK293 (OR 8.1 / 7.6 vs. 3.7 / 3.4),
  unidireccional (interacción clara solo en Rep 1), PyPu en −1/0 (más
  frecuente en las de **baja** actividad, al revés que en el resto).
- **Diferencia general alta vs. baja, no específica:** isla CpG, CCAAT, G+C.
- **Sin efecto:** Polycomb (cualquiera y cada uno), TSS no canónico, E-box,
  poli(dA:dT), Inr, TSS fuerte.
- Lectura: H1 "regulados" sí; H1 Polycomb/angostos no; H2 no; H3 no (PyPu al
  revés); H4 no específica.

### p7 — especificidad tisular continua (Gini de FANTOM5)
- Gini mediano TATA fuertes: 0.51 / 0.50 (alta) vs. 0.80 / 0.81 (baja);
  resto 0.34 / 0.35 vs. 0.44 / 0.45. Interacción (percentiles) p 0.05 / 0.012.
- Muestras de FANTOM5 activas: misma dirección, interacción n.s.
- Actividad vs. Gini: Spearman −0.41 / −0.44 en TATA fuertes vs. −0.26 / −0.28
  en el resto; pendiente −47 / −54 vs. −32 / −34 (p 0.026 / 0.001).
- Las TATA fuertes son más activas que el resto **en todos los quintiles de
  Gini** (q1 87 / 91 vs. 60 / 61; q5 55 / 50 vs. 39 / 39). El 46% de las TATA
  fuertes está en el quintil más específico.
- Lectura: las TATA fuertes de baja actividad son sobre todo TATA de
  promotores tejido-específicos; el motivo sigue sumando actividad, pero no
  alcanza el cuartil superior.

### p8 — ruido de las TATA fuertes según especificidad
- Todas las actividades juntas: AUC baja con la especificidad (tercil bajo
  0.64 / 0.66 → alto 0.50 / 0.53; p 0.02 / 0.01).
- **Separando por actividad desaparece:** alta actividad elevado en todos los
  terciles (0.63 / 0.67, 0.66 / 0.61, 0.59 / 0.68; p extremos 0.56 / 0.84);
  baja actividad ~0.5 en todos (tercil alto 0.465 / 0.486, n 128 / 170).
- Sin correlación ruido–Gini dentro de cada estrato (rho 0.03–0.07).
- Lectura: el efecto sobre el ruido sigue a la actividad, no a la clase de
  promotor → refuerza A.
- Límite: "TATA fuerte × expresión amplia × baja actividad" tiene n 6 / 12.

### p9 — perfil posicional y controles de composición dentro de cada estrato de actividad
Mismo análisis que p1–p2, con casos y referencia "sin motivo" del mismo
estrato (percentil 75 por réplica).

- **Actividad alta:** la ventana canónica sube a 0.596 / 0.603 (n 289 / 322)
  y el resto del perfil queda en ~0.5 (mediana fuera de la canónica 0.56 /
  0.51 misma hebra; las otras ventanas BH < 0.05 son de una sola réplica, n
  24–47 y en direcciones opuestas). Muchas ventanas no canónicas tienen n < 20
  y no se grafican.
- Clases, actividad alta: canónica 0.601 / 0.607; solo no canónica 0.529 /
  0.494; solo hebra opuesta 0.513 / 0.524; hebra opuesta en −31 a −26 sin
  canónica 0.515 / 0.488 (n 31 / 45). Todas las no canónicas n.s.
- Permutadas, actividad alta: ninguna reproduce el AUC en la ventana canónica
  (todas BH n.s.; n 20–56).
- Modelo A/T, actividad alta: canónica +10.2 / +10.8 rangos → +8.25 / +8.47
  con A/T proximal y del fragmento (se mantiene ~80%).
- **Actividad baja:** canónica 0.525 / 0.526 n.s. y +2.6 → ~0 con A/T; las
  clases no canónicas y las permutadas 4 y 6 sí dan AUC > 0.5 (0.58–0.63).
- Lectura: el componente de composición A/T sobre el ruido vive en los
  promotores de baja actividad; en los de alta actividad el efecto es
  específico del motivo y de la posición canónica, como para la actividad.

### p10 — `var_rank_sw ~ clase * actividad_alta + A/T` en un solo modelo
Por réplica, sobre todos los promotores. Efecto de la canónica vs. sin motivo
(rangos), baja / alta actividad:

| Modelo | Baja | Alta | p interacción |
|---|---|---|---|
| Sin A/T | +2.6 / +2.6 | +10.2 / +10.8 | 0.002 / <0.001 |
| + A/T proximal y del fragmento (pendiente común) | −0.9 / +0.1 (n.s.) | +7.6 / +8.7 | <0.001 / <0.001 |
| + A/T × actividad (sensibilidad) | −1.1 / +0.1 (n.s.) | +8.3 / +8.5 | <0.001 / <0.001 |

- El ajuste resta una cantidad parecida en ambos estratos (~2.5–3.5 rangos en
  baja, ~2–2.6 en alta), pero en baja era todo el efecto y en alta es ~20–25%
  (se retiene 75–81%).
- Con A/T la diferencia entre estratos no se achica; la interacción se
  mantiene o aumenta.
- Clases no canónicas: con A/T quedan en ~0 en ambos estratos.
- Lectura: A/T explica el pequeño efecto de baja actividad y una fracción
  menor del de alta; el efecto específico de la TATA canónica existe solo en
  alta actividad.

### p11 — ¿es razonable el corte en el percentil 75?
TATA canónica vs. sin motivo, efecto sobre el ruido en función de la
actividad (percentil dentro de cada réplica).

- **Por decil:** AUC ~0.49–0.54 en los deciles 1–6, 0.56 / 0.55 en el 7
  (n.s.), 0.60–0.63 en los deciles 8–10 (significativos). n 36–198 por decil.
- **Spline (`ns(act_pct, 4) * canonica`):** efecto ~0–3 rangos hasta P55;
  sube entre ~P60 y ~P80 (+7.3 / +6.2 en P65, ya con IC > 0; +12.0 / +10.2 en
  P75) y queda en meseta ~10–12 por encima de P80. Transición gradual, no
  escalón.
- **Barrido del corte (P50–P90):** el AUC por encima casi no cambia
  (0.59–0.61 con cualquier corte); el de abajo sube de 0.50–0.51 (cortes
  ≤ P65) a 0.54 (P90). Con P75 el grupo bajo queda en 0.525 / 0.526 (p 0.10 /
  0.06) porque incluye el decil 7.
- **Bisagra (efecto plano hasta c y lineal después):** perfil plano; un
  efecto lineal en toda la actividad ajusta igual que cualquier c ≤ ~60
  (Rep 1) / ~68 (Rep 2), y los quiebres más altos ajustan peor. Los datos no
  identifican un umbral; el modelo bisagra está mal especificado para una
  curva con meseta, así que la descripción buena es el spline.
- Lectura: P75 es coherente (no cambia ninguna conclusión y es un corte
  elegido a priori, sin sobreajuste), pero no es "el" umbral: el efecto
  aparece de forma gradual entre ~P60 y ~P80. Un corte en P65–P70 separa un
  poco más limpio (abajo ~0.51–0.52). Recomendado para la tesis: mantener P75
  y mostrar la curva continua (deciles o spline) como respaldo.
- El panel de deciles se guarda además como respaldo (no es panel de la
  figura ni va a la tesis) en la carpeta de `tata_ruido_especificidad`:
  `respaldo_auc_por_decil_actividad.jpg/.pdf/.tsv`, vía `fig_dir()`. No está
  en la columna `archivos` del mapping.

---

## Pendientes

- **Tesis (Google Doc, no `docs/tesis.md`, que se sobreescribe al sincronizar):**
  - Resultados R3.2: "AUC de 0,39 y 0,43" → 0,40 y 0,44.
  - Métodos "Asociación de las características del promotor con el ruido":
    agregar que el binning es por réplica y descarta el resto de n/100, cómo
    se resuelven empates (orden de fila) y el criterio de actividad alta
    (percentil 75 de la media dentro de cada réplica). Texto sugerido en la
    conversación del 2026-10-27.
- **Panel propuesto para ruido_tata:** AUC por celda puntaje × actividad
  (p4_puntaje_x_actividad), como forest horizontal en el formato de
  `tata_controles_modelo_AT`. No implementado.
- `R/42_exploratorio_ruido_composicion`: `p5_ROC_gc_alto.jpg` quedó con el
  AUC por trapecio (no se volvió a correr).
- `summary_noise.jpg` (ruido_summary): la barra de error de Rep 1 se ve
  corrida respecto de su barra; ya pasaba antes, no se tocó.
- Posible siguiente paso: posición fina de la TATA dentro de −31 a −26 por
  estrato de actividad (la parte "mal posicionadas" de la lectura B no se
  evaluó más allá de canónica vs. no canónica).
