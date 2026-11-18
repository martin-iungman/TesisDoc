# Exploratorio: embeddings de modelos de lenguaje de ADN (2026-11-17)

## Pregunta

Cuánto de la actividad media explicable de cada promotor puede predecirse desde la secuencia de 252 pb
sola, y si una representación preentrenada sobre el genoma (DNABERT-2, Nucleotide Transformer v2,
HyenaDNA) agrega algo sobre bases simples. Destino posible: el techo de varianza explicable de la
Discusión.

## Archivos

- `preparar_dataset.py` → `data/processed/embeddings_adn/embeddings_adn_dataset_v1.tsv`: una fila por promotor medido
  en al menos una réplica (11955; 9397 en Rep 1, 11447 en Rep 2, 8889 en ambas), con la secuencia
  orientada (copiada de `transcriptional_library/Library_data/data/promoters_wo_dupl.fa`; la TATA cae
  en −31 a −27, así que está en el sentido de la transcripción), actividad media por réplica,
  las 8 características del modelo aditivo de R/38, G+C, CpG o/e y el fold.
- `evaluar.py`: bases y función `evaluar_embeddings()`. Salida en `figures/exploratorio_embeddings_adn/`.
- `extraer_embeddings_colab.ipynb`: extracción de embeddings en Colab con GPU y evaluación con
  `evaluar.py` (escala cruda). Ya corrido: no repetir. Los embeddings quedaron en Drive
  (`embeddings_adn/emb_*.npz`) y se bajan a `data/processed/embeddings_adn/`.
- `mejoras.py`: corre local, sin GPU, sobre los embeddings guardados. Todo en percentiles. Guarda las
  predicciones fuera de fold de cada modelo en `figures/exploratorio_embeddings_adn/mejoras/predicciones/`
  y no recalcula un modelo si su archivo ya existe (`--rehacer` para forzar, `--listar` para ver el
  catálogo, `--solo a,b` para correr algunos). Escribe `resultados.tsv` y `bootstrap.tsv`.
  Con `--objetivo promedio` entrena y evalúa contra el promedio de percentiles de las dos réplicas
  (ver "Objetivo: promedio de percentiles") y escribe en `mejoras_promedio/`, con el mismo formato.
  Requiere `scipy` y `scikit-learn` (instalados con pip el 2026-11-17).

- `finetuning.py` + `finetuning_colab.ipynb` (2026-11-18, corrido en Colab con GPU; ver "Fine-tuning"): fine-tuning de
  DNABERT-2 contra una CNN ajustada (diseño en el encabezado de `finetuning.py`). Multitarea con tres
  salidas (percentil Rep 1, Rep 2 y promedio), así una corrida sirve para los dos modos de
  `mejoras.py`. Hiperparámetros elegidos en una validación interna por cromosoma dentro del
  entrenamiento (CNN: 12 configuraciones al azar; DNABERT-2: lr 2e-5 o 5e-5), 3 semillas. Guarda cada
  (modelo, fold, semilla) en Drive y retoma si Colab se desconecta. Las predicciones salen con el
  formato de `mejoras.py` (`cnn_ajustada`, `dnabert2_ft`, más `_s0` = una semilla); `mejoras.py` las
  lee como modelos externos si están en `predicciones/` (ver la última celda del notebook). Probado
  local en CPU con un subconjunto (`--prueba`); en local, `transformers` tiene que ir en una ruta corta
  por el límite de 260 caracteres de Windows.

**Regla de organización (Martín, 2026-11-17):** Colab solo para lo que necesita GPU (extracción de
embeddings, fine-tuning). Todo lo demás corre local y deja sus salidas en TesisDoc.

## Diseño

- **Split:** 5 folds por cromosoma, fijos en el dataset (fold 0: chr1, 5, 22, X; 1: chr7, 13, 16, 19, 20;
  2: chr2, 6, 10, 14, 21; 3: chr4, 8, 12, 17; 4: chr3, 9, 11, 15, 18). Agrupa también por gen.
- **Evaluación:** predicción fuera de fold, cada réplica por separado. `pearson_otra_replica` correlaciona
  la predicción entrenada con una réplica contra la otra, en los 8889 medidos en ambas.
- **Techo:** r entre réplicas = 0.865. Con un predictor perfecto, el R² esperado contra una réplica es
  ~0.865 y la correlación ~0.93 (no 0.865). `frac_techo` = R² / 0.865.
- **Ridge:** variables estandarizadas, α elegido por GCV dentro del entrenamiento. Para embeddings, la
  capa también se elige dentro del entrenamiento de cada fold.

## Escala

Las dos réplicas tienen distribuciones distintas (Rep 2 corrida hacia abajo: 35% de los promotores
bajo 1.5 contra 3% en Rep 1) y ambas una cola larga de promotores muy activos. Pearson entre réplicas
sobre la media cruda: 0.865; sobre percentiles (= Spearman): 0.688. El 0.865 está dominado por la
cola. **Escala principal decidida: percentil dentro de cada réplica** (representa al promotor típico,
no depende de la calibración de cada réplica, coherente con el resto de la tesis). La escala cruda
queda como sensibilidad. Entrenar con el promedio de percentiles de las dos réplicas pero evaluar
contra cada réplica por separado no mejora nada (aditivo 0.148 contra 0.148, k-meros por ventana
0.215 contra 0.212 de fracción del techo). Entrenar **y evaluar** contra el promedio es otra cosa
(cambia el techo): ver "Objetivo: promedio de percentiles".

En percentiles (R² / fracción del techo 0.688): aditivo 0.102 / 0.096 (15% / 14%); aditivo + composición
0.121 / 0.112 (17% / 16%); k-meros por ventana 0.146 / 0.135 (21% / 20%).

DNABERT-2 y los otros modelos en percentiles (`embeddings_percentil.tsv`, celda agregada al notebook de
Colab, capas restringidas a las ganadoras en escala cruda): DNABERT-2 ventanas 0.190 / 0.179 (28% / 26%
de lo explicable); DNABERT-2 media 0.151 / 0.148 (22% / 22%); NT v2 ventanas 0.145 / 0.142 (21% / 21%);
HyenaDNA ventanas 0.142 / 0.135 (21% / 20%). El orden es el mismo que en escala cruda: NT y HyenaDNA
quedan a la par de los k-meros por ventana y DNABERT-2 con ventanas los supera por ~6 puntos de
fracción del techo. Ese notebook no guardó predicciones fuera de fold, así que `mejoras.py` vuelve a
ajustar estos cuatro modelos (necesita las predicciones para el bootstrap); sus números deberían coincidir
con estos y sirven de chequeo.

## Resultados de `mejoras.py` por réplica (`mejoras/resultados.tsv`, `mejoras/bootstrap.tsv`)

Percentil dentro de cada réplica, techo 0.688. Corrido el 2026-11-17 (~85 min local). Los cuatro
modelos ya corridos en Colab reproducen `embeddings_percentil.tsv` (DNABERT-2 ventanas 0.192 / 0.178
contra 0.190 / 0.179).

| Modelo | R² Rep 1 / Rep 2 | Fracción del techo |
|---|---|---|
| DNABERT-2 ventanas, capas 3 + 6 | 0.200 / 0.188 | 0.29 / 0.27 |
| DNABERT-2 c3 + k-meros 1 a 4 por ventana + aditivo_comp | 0.196 / 0.182 | 0.29 / 0.27 |
| DNABERT-2 ventanas c3 | 0.192 / 0.178 | 0.28 / 0.26 |
| HyenaDNA + NT v2 + DNABERT-2, ventanas | 0.191 / 0.180 | 0.28 / 0.26 |
| Gradient boosting: PCA(128) de DNABERT-2 c3 + aditivo_comp | 0.165 / 0.155 | 0.24 / 0.23 |
| DNABERT-2 media c6 | 0.155 / 0.150 | 0.23 / 0.22 |
| k-meros 1 a 5 por ventana (4092 variables) | 0.155 / 0.147 | 0.23 / 0.21 |
| k-meros 1 a 4 por ventana (1020) | 0.146 / 0.135 | 0.21 / 0.20 |
| NT v2 ventanas c11 | 0.145 / 0.142 | 0.21 / 0.21 |
| HyenaDNA ventanas c2 | 0.142 / 0.135 | 0.21 / 0.20 |
| Aditivo + composición | 0.121 / 0.112 | 0.18 / 0.16 |
| Aditivo | 0.102 / 0.096 | 0.15 / 0.14 |

Lecturas (IC 95% por bootstrap de la diferencia de R² contra k-meros 1 a 4 por ventana):
1. La ventaja de DNABERT-2 no es por cantidad de variables: k ≤ 5 suma solo 0.008 / 0.013 a la base
   (significativo), y DNABERT-2 c3 le sigue sacando ~0.035. DNABERT-2 c3 contra k-meros 1 a 4:
   +0.046 [0.036, 0.055] / +0.043 [0.036, 0.051]. (Contra k ≤ 5 no hay IC: el bootstrap usa solo
   k-meros 1 a 4 y aditivo como referencias.)
2. Sumar k-meros y aditivo a DNABERT-2 agrega ~0.004: DNABERT-2 ya contiene esa información. Sumar
   HyenaDNA y NT tampoco agrega.
3. HyenaDNA y NT v2 no se distinguen de los k-meros 1 a 4 (IC incluyen o tocan el 0).
4. El gradient boosting rinde menos que ridge sobre el mismo embedding; probablemente el PCA pierde
   información. No es evidencia de no linealidad aprovechable.
5. Capas 3 + 6 es lo mejor, apenas sobre la capa 3 sola. Nada pasa del ~29% de lo explicable.

## Objetivo: promedio de percentiles (`mejoras_promedio/`)

Solo los 8889 promotores medidos en las dos réplicas. Percentil dentro de cada réplica, calculado en
ese subconjunto; se entrena y evalúa contra el promedio de los dos. Techo = confiabilidad del promedio
por Spearman-Brown, 2r/(1+r) = 0.812 (r = 0.684 entre percentiles en el subconjunto). Los R² absolutos
suben porque el objetivo tiene menos ruido; la fracción del techo es la comparable.

Corrido el 2026-11-17 (~40 min local).

| Modelo | R² | Fracción del techo | Dif. R² contra k-meros 1 a 4 [IC 95%] |
|---|---|---|---|
| DNABERT-2 ventanas, capas 3 + 6 | 0.264 | 0.33 | +0.069 [0.058, 0.079] |
| DNABERT-2 c3 + k-meros 1 a 4 por ventana + aditivo_comp | 0.263 | 0.32 | +0.068 [0.059, 0.077] |
| DNABERT-2 ventanas c3 | 0.257 | 0.32 | +0.062 [0.052, 0.072] |
| HyenaDNA + NT v2 + DNABERT-2, ventanas | 0.255 | 0.31 | +0.060 [0.048, 0.070] |
| Gradient boosting: PCA(128) de DNABERT-2 c3 + aditivo_comp | 0.222 | 0.27 | +0.027 [0.016, 0.038] |
| k-meros 1 a 5 por ventana | 0.208 | 0.26 | +0.013 [0.009, 0.017] |
| DNABERT-2 media c6 | 0.207 | 0.26 | +0.012 [0.000, 0.023] |
| NT v2 ventanas c11 | 0.202 | 0.25 | +0.007 [−0.003, 0.017] |
| k-meros 1 a 4 por ventana | 0.195 | 0.24 | — |
| HyenaDNA ventanas c2 | 0.195 | 0.24 | −0.000 [−0.011, 0.009] |
| Aditivo + composición | 0.159 | 0.20 | −0.036 [−0.048, −0.025] |
| Aditivo | 0.135 | 0.17 | −0.060 [−0.073, −0.048] |

Lecturas:
1. El orden y las conclusiones son los mismos que por réplica: DNABERT-2 con ventanas supera a todas
   las bases (también a k ≤ 5, por ~0.05 de R²), sumarle k-meros, aditivo u otros modelos casi no
   agrega, HyenaDNA y NT v2 no se distinguen de los k-meros, gradient boosting rinde menos que ridge.
2. Las diferencias son algo más grandes y los IC algo más angostos que por réplica (objetivo con menos
   ruido).
3. La fracción del techo sube en todos los modelos (DNABERT-2 c3: 0.32 contra 0.28 / 0.26 por
   réplica). Si Spearman-Brown valiera exacto, debería coincidir. No está resuelto por qué no
   coincide; candidatos: el subconjunto medido en ambas réplicas no es el mismo conjunto de
   promotores, o los errores de las réplicas no son independientes ni de igual confiabilidad. Por eso
   la fracción del techo no se compara directamente entre los dos modos: el techo de "lo explicable
   desde la secuencia" queda entre ~29% (por réplica) y ~33% (promedio).

## Fine-tuning: DNABERT-2 contra CNN ajustada (2026-11-18)

Predicciones de `finetuning_colab.ipynb` en `mejoras/predicciones/` y `mejoras_promedio/predicciones/`;
resumen y bootstrap rehechos con `mejoras.py` en los dos modos (referencias del bootstrap ampliadas a
`cnn_ajustada` y `dnabert2_vent_c3c6`). Sin ensamblar = promedio de 3 semillas; `_s0` = una semilla.

| Modelo | R² Rep 1 / Rep 2 | Fracción del techo | R² promedio | Fracción del techo |
|---|---|---|---|---|
| DNABERT-2 fine-tuning, 3 semillas | 0.220 / 0.210 | 0.32 / 0.31 | 0.294 | 0.36 |
| DNABERT-2 fine-tuning, 1 semilla | 0.206 / 0.192 | 0.30 / 0.28 | 0.277 | 0.34 |
| DNABERT-2 congelado, ventanas c3 + c6 (ridge) | 0.200 / 0.188 | 0.29 / 0.27 | 0.264 | 0.33 |
| DNABERT-2 congelado, ventanas c3 (ridge) | 0.192 / 0.178 | 0.28 / 0.26 | 0.257 | 0.32 |
| CNN ajustada, 3 semillas | 0.180 / 0.170 | 0.26 / 0.25 | 0.240 | 0.30 |
| CNN ajustada, 1 semilla | 0.171 / 0.155 | 0.25 / 0.23 | 0.223 | 0.27 |
| k-meros 1 a 4 por ventana (ridge) | 0.146 / 0.135 | 0.21 / 0.20 | 0.195 | 0.24 |

Diferencias de R² con IC 95% por bootstrap (Rep 1 / Rep 2 / promedio):
- DNABERT-2 fine-tuning contra CNN, 3 semillas cada uno: +0.040 [0.028, 0.053] / +0.040 [0.030, 0.050]
  / +0.054 [0.040, 0.068]. Con 1 semilla cada uno, +0.026 / +0.022 / +0.037 (IC sin incluir el 0;
  `_s0` contra `_s0` sale de restar, la referencia del bootstrap es la CNN de 3 semillas).
- DNABERT-2 fine-tuning contra DNABERT-2 congelado c3 + c6: con 3 semillas +0.020 [0.011, 0.030] /
  +0.022 [0.014, 0.030] / +0.030 [0.018, 0.041]; con 1 semilla +0.006 [−0.005, 0.016] / +0.004
  [−0.006, 0.013] / +0.013 [−0.001, 0.025].
- CNN ajustada contra k-meros 1 a 4: +0.034 / +0.035 / +0.045 (3 semillas), IC sin incluir el 0.
- DNABERT-2 congelado c3 + c6 contra CNN ajustada de 3 semillas: +0.020 / +0.018 / +0.024, IC sin
  incluir el 0.

Lecturas:
1. **El preentrenamiento aporta frente a una base supervisada seria.** DNABERT-2 supera a la CNN
   ajustada en los tres objetivos, comparando igual contra igual (3 contra 3 semillas, o 1 contra 1).
   Más aún: DNABERT-2 **congelado** con ridge ya le gana a la CNN ajustada de 3 semillas.
2. **El fine-tuning suma poco sobre los embeddings congelados.** Con una semilla no se distingue de
   ridge sobre capas 3 + 6 (IC incluyen el 0). La ventaja de la versión de 3 semillas (~0.02 a 0.03)
   es en buena parte el efecto de promediar modelos: promediar 3 semillas de la CNN también le suma
   0.009 a 0.017. Para comparar contra ridge (un único modelo determinista), la cifra justa es `_s0`.
3. La CNN ajustada queda bien por encima de la CNN sin ajustar de `bases.tsv` y de los k-meros, así
   que es una base razonable, no un hombre de paja.
4. Hiperparámetros elegidos (`finetuning/<modelo>/busqueda_fold*.json`, bajados de Drive; los `.npy`
   y los `.json` por semilla quedaron solo en Drive, las predicciones ya están ensambladas):
   - CNN: en los 5 folds ganan 4 bloques (el máximo del espacio) y casi siempre el kernel inicial más
     chico (7, el mínimo); 14 a 27 épocas. Que gane el borde del espacio sugiere que una CNN más
     profunda o con kernels más chicos podría rendir algo más: la CNN es una base razonable, no
     necesariamente la mejor posible.
   - DNABERT-2: lr 2e-5 en 3 folds y 5e-5 en 2; paró en 3 a 5 épocas (de 8), así que la parada
     temprana actuó. En el fold 4, lr 5e-5 se desestabilizó (val 0.98 contra 0.73).
5. Techo: lo mejor (fine-tuning, 3 semillas) llega a ~31 a 32% de lo explicable por réplica y ~36%
   contra el promedio. Con un único modelo, ~29% y ~34%. La conclusión para la Discusión no cambia:
   desde los 252 pb solos se explica alrededor de un tercio de la varianza explicable.

## Resultados de las bases en escala cruda (`bases.tsv`)

| Modelo | Variables | R² Rep 1 | R² Rep 2 | Fracción del techo |
|---|---|---|---|---|
| Aditivo (8 elementos) | 8 | 0.109 | 0.099 | 0.13 / 0.12 |
| Aditivo + G+C y CpG o/e | 10 | 0.132 | 0.119 | 0.15 / 0.14 |
| k-meros 1 a 5 | 1364 | 0.133 | 0.126 | 0.15 / 0.15 |
| k-meros 1 a 4 en tres ventanas | 1020 | 0.164 | 0.155 | 0.19 / 0.18 |
| CNN desde cero (sin ajuste fino) | | 0.115 | 0.114 | 0.13 / 0.13 |

Lecturas:
1. El R² del modelo aditivo fuera de fold (0.109 / 0.099) es casi igual al ajustado sobre todos los
   datos en R/38 (0.106 / 0.097): ese número no estaba inflado por sobreajuste.
2. La composición continua alcanza para igualar a todos los k-meros sin posición: sin información
   posicional, la secuencia aporta poco más allá de la composición y los elementos conocidos.
3. La posición sí agrega: los k-meros por ventana (−236 a −61, −60 a −11, −10 a +15) son lo mejor,
   alrededor de un 50% más de R² que el modelo aditivo. Aun así, todo queda por debajo del 20% de lo
   explicable.
4. La CNN no está ajustada (arquitectura y tasa de aprendizaje fijas, para con 7 a 11 épocas porque
   sobreajusta rápido). Con ~6000 a 7000 secuencias de entrenamiento por fold es un mal punto de
   comparación todavía: tomarla como cota inferior, no como "lo que da una CNN".

## Resultados de los embeddings en escala cruda (`embeddings.tsv`, `embeddings_por_capa.tsv`, corrido en Colab)

Ridge sobre embeddings congelados; capa elegida dentro del entrenamiento de cada fold.

| Modelo | R² Rep 1 | R² Rep 2 | Fracción del techo | Capa elegida |
|---|---|---|---|---|
| HyenaDNA, media | 0.116 | 0.126 | 0.13 / 0.15 | 1 de 5 |
| HyenaDNA, ventanas | 0.151 | 0.144 | 0.17 / 0.17 | 2 |
| NT v2 100M, media | 0.131 | 0.120 | 0.15 / 0.14 | 11 de 22 |
| NT v2 100M, ventanas | 0.154 | 0.147 | 0.18 / 0.17 | 11 |
| DNABERT-2, media | 0.157 | 0.156 | 0.18 / 0.18 | 3 a 6 de 12 |
| **DNABERT-2, ventanas** | **0.216** | **0.202** | **0.25 / 0.23** | 3 |

Lecturas:
1. **DNABERT-2 con ventanas es lo mejor de todo**: ~30% más de R² que los k-meros por ventana
   (0.164 / 0.155) y el doble que el modelo aditivo, consistente en las dos réplicas y contra la otra
   réplica (r 0.50 / 0.47). Aun así llega a ~25% de lo explicable.
2. HyenaDNA y NT v2 no superan a los k-meros por ventana: rinden como una base de composición y
   posición gruesa.
3. En los tres modelos las capas tempranas o intermedias rinden más que la última (la última capa está
   especializada en la tarea de preentrenamiento). Conservar la posición (ventanas) mejora siempre.
4. Comparación no del todo pareja: DNABERT-2 con ventanas tiene 2304 variables contra 1020 de los
   k-meros (k ≤ 4). Ridge regulariza, pero conviene probar k-meros por ventana con k ≤ 5 antes de
   atribuir la ganancia al preentrenamiento.

## Pendiente

- Ventanas más finas o más variantes de fine-tuning: baja prioridad. El fine-tuning de un único modelo
  no se distingue de ridge sobre embeddings congelados; la representación parece saturada.
- Para la Discusión: Tang et al. 2025 (Genome Biology) y DART-Eval (Patel et al. 2025) encuentran que
  los embeddings congelados de modelos de lenguaje de ADN no aportan ventajas sustanciales sobre bases
  convencionales en tareas regulatorias. Acá eso vale para HyenaDNA y NT v2; DNABERT-2 sí supera a las
  bases, incluida una CNN ajustada (ver "Fine-tuning"), y congelado rinde casi igual que con
  fine-tuning.
- Opcional: agregar a las secuencias los flancos fijos del vector, como se hizo para el predscore de PUFFIN.
