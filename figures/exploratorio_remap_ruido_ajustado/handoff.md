# Control por actividad endógena y accesibilidad — unión de factores de ReMap y ruido

Script: `R/48_exploratorio_remap_ruido_ajustado/remap_ruido_ajustado.R`.
Salidas: esta carpeta, con prefijo por punto del pedido (p1–p4).
Todos los resultados se dan como Rep 1 / Rep 2.

## Método

- **Factores.** Son los 982 de `remap_tf_noise_auc.tsv` (remap_tf_noise_gsea).
  - La unión se toma de `data/processed/remap_tf_hits.tsv`, con picos en cualquier línea celular, igual que en el análisis original.
- **Promotores.** Es el set highconf con `var_rank_sw`, el rango de varianza dentro del bin de actividad (1–100). Son 9300 / 11400 promotores.
- **Covariables.** Usé las versiones continuas que tiene `prom_df`, no las binarizaciones del panorama:
  - `hek_tpm`: CAGE en HEK293, transformada como log10(TPM + 1).
  - `mean_dnase`: DNase-seq.
  - No hay NAs, así que los tres modelos usan exactamente los mismos promotores.
- **Modelos.** Ajusté tres modelos lineales por factor y por réplica:
  - crudo: `var_rank_sw ~ unido`
  - solo_actividad: `var_rank_sw ~ unido + log10(hek_tpm + 1)`
  - ajustado: `var_rank_sw ~ unido + log10(hek_tpm + 1) + mean_dnase`
- **Estadística.** El coeficiente de `unido` se reporta con su IC 95% (t). El p está corregido por BH dentro de cada réplica y de cada modelo.
- **Consistencia.** Un factor es consistente si tiene p corregido < 0.05 en ambas réplicas y el coeficiente tiene el mismo signo en las dos.
- **GSEA.** Usé los mismos parámetros que remap_tf_noise_gsea: `gseGO(ont = "CC")`, `set.seed(2024)` antes de cada tanda de réplicas, `seed = TRUE` y el resto por defecto.
  - La lista se ordena de mayor a menor métrica. NES < 0 quiere decir que el término se concentra en el extremo de menor ruido.
  - Lo corrí sobre cuatro rankings: la AUC original y el coeficiente de cada uno de los tres modelos.
  - Usé `pvalueCutoff = 1` solo para poder listar todos los términos. El p ajustado no cambia por eso, porque BH se calcula igual sobre todos los términos testeados.
  - El ranking por AUC reproduce exactamente el GSEA actual. Por ejemplo, MLL1/2 da 0.0046 / 0.081 y nuclear protein-containing complex da 2.09e-06 / 0.0037.

### Advertencia metodológica: la accesibilidad es en parte mediador

- **Por qué el modelo ajustado es conservador.** La accesibilidad de la cromatina es en parte consecuencia de la unión de factores, no solo un confusor previo. Al ajustar por `mean_dnase` se descuenta también parte del efecto que se quiere medir. Por eso el modelo ajustado es conservador por construcción. Se reportan siempre las dos versiones y la ajustada no es «la correcta».
- **Por qué agregué el modelo solo_actividad.** La actividad endógena en HEK293 (CAGE) también es en parte consecuencia de la unión. Aun así, está más lejos del sitio de unión que la accesibilidad, porque integra todo lo que regula al gen en su locus y no solo la apertura local. Por eso el modelo solo_actividad es un punto intermedio razonable: controla que «unido» signifique «gen activo en HEK293» sin descontar directamente la apertura que genera el factor.
- **Correlación entre covariables.** `log_hek` y `mean_dnase` están bastante correlacionadas (Spearman 0.66 / 0.67). Las dos se asocian débilmente al ruido:
  - actividad endógena: ρ = −0.14 / −0.02
  - accesibilidad: ρ = −0.18 / −0.05

  En la réplica 2 la asociación es mucho más débil, y eso explica casi todas las diferencias entre réplicas que siguen (`p1_covariables_correlacion.tsv`).

## Resultados

### p2. Cuántos factores mantienen el efecto

Factores con efecto significativo y del mismo signo en ambas réplicas, sobre 982:

| Ranking / modelo | Menor ruido | Mayor ruido |
|---|---|---|
| AUC original | 514 | 0 |
| crudo (lm) | 502 | 0 |
| solo_actividad | 443 | 1 |
| ajustado | 276 | 4 |

- **Mayor ruido con solo_actividad:** el único es GTF2A2.
- **Mayor ruido con ajustado:** BATF, GTF2A2, SNAPC1 y ZC3H8. Los efectos son de +2 a +4 posiciones y todos tienen pocos promotores unidos, salvo SNAPC1.

**Reducción del efecto.** Se calcula entre los 502 factores consistentes en el modelo crudo (`p2_reduccion_efecto.tsv`):

| | Rep 1 | Rep 2 |
|---|---|---|
| Coeficiente mediano, crudo | −6.80 | −2.70 |
| Coeficiente mediano, solo_actividad | −4.37 | −2.85 |
| Coeficiente mediano, ajustado | −2.19 | −2.40 |
| Reducción mediana (1 − ajustado/crudo), solo_actividad | 35% | −5% (no se reduce) |
| Reducción mediana, ajustado | 66% | 11% |

- En la réplica 1 el ajuste se lleva una parte grande del efecto, pero no todo.
- En la réplica 2 el efecto casi no cambia, porque las covariables apenas se asocian al ruido en esa réplica.
- **Conclusión de p2.** «La unión se asocia a menor ruido» es en parte, pero no del todo, una forma indirecta de decir «promotores activos/accesibles». Después del ajuste completo, 276 factores mantienen el efecto en ambas réplicas, todavía con una asimetría casi total hacia menor ruido.

### p3. GSEA crudo vs. ajustado

La tabla completa está en `p3_gsea_comparacion.tsv` y todos los términos en `p3_gsea_todos_los_terminos.tsv`. Cada celda es NES (p ajustado), Rep 1 / Rep 2:

| Término | AUC original | crudo (lm) | solo_actividad | ajustado |
|---|---|---|---|---|
| MLL1/2 complex (GO:0044665) | −2.47 (0.0046) / −2.01 (0.081) | −1.64 (0.046) / −1.46 (0.19) | −1.56 (0.035) / −1.47 (0.11) | −1.76 (0.018) / −1.46 (0.18) |
| MLL1 complex (GO:0071339) | −2.37 (0.0052) / −1.86 (0.13) | −1.58 (0.069) / −1.39 (0.21) | −1.50 (0.095) / −1.41 (0.19) | −1.73 (0.018) / −1.40 (0.23) |
| SWI/SNF superfamily-type complex (GO:0070603) | −1.55 (0.20) / −2.11 (0.057) | −1.42 (0.10) / −1.37 (0.19) | −1.24 (0.29) / −1.36 (0.19) | −1.39 (0.15) / −1.42 (0.14) |
| ATPase complex (GO:1904949) | idéntico a SWI/SNF superfamily (mismos 44 factores en la lista) | | | |
| SWI/SNF complex (GO:0016514) | −1.71 (0.14) / −1.57 (0.28) | −1.82 (0.012) / −1.46 (0.19) | −1.69 (0.0099) / −1.49 (0.11) | −1.62 (0.042) / −1.42 (0.19) |
| nuclear protein-containing complex (GO:0140513) | −3.39 (2.1e-06) / −2.60 (0.0037) | −1.75 (1.2e-07) / −1.37 (0.030) | −1.50 (3.7e-06) / −1.34 (0.021) | −1.32 (0.018) / −1.32 (0.040) |
| catalytic complex (GO:1902494) | −2.29 (0.0061) / −2.68 (0.0022) | −1.52 (0.0035) / −1.37 (0.11) | −1.34 (0.024) / −1.34 (0.063) | −1.31 (0.12) / −1.39 (0.040) |
| protein-containing complex (GO:0032991) | −3.79 (8.3e-08) / −1.65 (0.23) | −1.79 (5.8e-09) / −1.25 (0.12) | −1.53 (4.7e-07) / −1.28 (0.062) | −1.30 (0.021) / −1.23 (0.11) |
| chromatin (GO:0000785) | −3.91 (5.8e-08) / −1.90 (0.10) | −1.78 (5.8e-09) / −1.34 (0.030) | −1.46 (4.2e-06) / −1.35 (0.021) | −1.26 (0.032) / −1.25 (0.11) |
| transferase complex (GO:1990234) | −2.28 (0.0081) / −1.92 (0.095) | −1.63 (0.0027) / −1.44 (0.12) | −1.52 (0.0046) / −1.43 (0.068) | −1.54 (0.018) / −1.47 (0.048) |

**Lectura de la tabla**

- **Lo que más cambia es la métrica, no el ajuste.** La diferencia grande está entre la columna AUC y las tres columnas de coeficiente. El paso de crudo a ajustado cambia mucho menos. Para MLL1/2 en la réplica 1, el ajuste incluso mejora el p: 0.046 → 0.018.
- **Por qué cambia tanto la métrica.** fgsea pondera cada factor por el valor de su métrica (gseaParam = 1).
  - La AUC es siempre positiva y está concentrada entre 0.38 y 0.53, así que los pesos son casi uniformes y el test se parece a un KS sin ponderar. fgsea lo advierte: «all values in the stats vector are greater than zero».
  - El coeficiente lineal tiene signo y un rango amplio, así que pesa mucho más a los factores extremos.
  - Por eso la comparación limpia del efecto del ajuste es crudo (lm) vs. ajustado, no AUC vs. ajustado.
- **Términos genéricos.** Comparados con el criterio de p ajustado < 0.1 en ambas réplicas:
  - Nuclear protein-containing complex lo cumple en las cuatro versiones, con NES más chico después del ajuste.
  - Transferase complex lo cumple en las cuatro, y con el ajuste completo da 0.018 / 0.048.
  - Catalytic complex lo cumple con AUC y con solo_actividad, pero no con crudo (lm) ni con ajustado.
  - Chromatin lo cumple con crudo (lm) y con solo_actividad.
  - Protein-containing complex lo cumple solo con solo_actividad.
- **SWI/SNF.** No es robusto: no llega a p < 0.1 en ninguna réplica en ninguna de las versiones ajustadas, salvo el término chico SWI/SNF complex (GO:0016514) en la réplica 1.

### p4. Miembros de MLL1/2 complex

El término GO:0044665 tiene 19 miembros con datos en ReMap (`p4_mll1_2_miembros.tsv`, `p4_mll1_2_miembros.jpg`). El término de GO es amplio: además del núcleo COMPASS/MLL incluye factores asociados de función general, como TAF1, MAX, E2F6, MGA y RNF2.

- **Sostienen el efecto en ambas réplicas incluso con el ajuste completo (13 de 19).** Coeficiente ajustado, Rep 1 / Rep 2:
  - Núcleo COMPASS/MLL: KMT2A −3.6 / −3.2, KMT2B −5.0 / −4.9, WDR5 −3.7 / −4.2, RBBP5 −4.1 / −4.4, ASH2L −3.3 / −2.8.
  - Asociados: HCFC1 −3.0 / −3.7, PHF20 −3.3 / −2.5, E2F6 −5.4 / −4.2, MAX −4.8 / −5.0, MGA −3.5 / −2.1, RNF2 −3.6 / −3.8, TAF1 −5.9 / −5.4, RUVBL2 −2.2 / −2.6.
- **Pierden el efecto con el ajuste completo.**
  - MEN1 (+0.3 / −0.9), CHD8 (−0.2 / −1.3) y TAF7 (−0.5 / −0.1): eran significativos en crudo y siguen siéndolo con solo_actividad, salvo TAF7 en la réplica 2.
  - MCRS1 (−1.6 / −1.9): queda significativo solo en la réplica 2.
  - RUVBL1 y KAT8: no eran consistentes ni en crudo. KAT8 tiene 143 / 150 promotores unidos y cambia de signo entre réplicas.
- **Conteo de miembros consistentes:** 16 de 19 en crudo, 16 de 19 con solo_actividad y 13 de 19 con el ajuste completo.
- **Posición frente al resto de los factores.** Con el ajuste completo, la mediana de los miembros está en el percentil 0.14 / 0.22 del ranking (0 = menor ruido). Su coeficiente mediano es −3.34 / −2.83, contra −1.08 / −1.05 del total de factores.
- **Quién sostiene el enriquecimiento.** No son unos pocos miembros. Lo sostiene el núcleo catalítico COMPASS completo (KMT2A, KMT2B, WDR5, RBBP5, ASH2L) más los asociados de unión amplia.
- **MEN1 cae.** Es el componente específico de MLL1/2 frente a los otros COMPASS (SET1A/B). Pierde el efecto con el ajuste completo, aunque lo mantiene con solo_actividad.

## Lectura: ¿se sostiene el resultado de MLL1/2?

1. **Réplica 1.** Se sostiene. MLL1/2 sigue enriquecido en el extremo de menor ruido con las tres versiones del coeficiente: p ajustado 0.046 crudo, 0.035 solo_actividad y 0.018 ajustado.
2. **Réplica 2.** No llega a significativo con ninguna versión del coeficiente (0.11–0.19), tampoco sin ajustar. Ya con la AUC era borderline (0.081). El ajuste no es lo que lo debilita: en la réplica 2 el ajuste casi no cambia nada. Lo que cambia es la métrica de ranking.
3. **Por miembro.** Es la evidencia más sólida. El núcleo COMPASS/MLL completo (KMT2A, KMT2B, WDR5, RBBP5, ASH2L) mantiene, en ambas réplicas y con el ajuste completo, un efecto de 2 a 3 veces el efecto mediano del resto de los factores.
4. **Confusor posible no controlado: amplitud de unión.** El efecto de cada factor se correlaciona fuerte con cuántos promotores une, incluso con el modelo ajustado: Spearman entre n unidos y coeficiente = −0.63 / −0.66. Los miembros de MLL1/2 unen muchos más promotores que el factor típico (mediana 6191 / 7376 contra 3554 / 4164, sobre 9300 / 11400). Con factores que unen ~70–80% de los promotores, la comparación es en la práctica «el 20–30% no unido contra el resto», y ese grupo no unido es particular.

   Parte del enriquecimiento de MLL1/2 puede reflejar que son factores de unión amplia (regiones HOT). Si la figura siguiente depende de este resultado, sugiero un control adicional:
   - agregar el número total de factores unidos por promotor como covariable, o
   - comparar MLL1/2 contra factores de amplitud de unión similar.

**Recomendación para la tesis.**
- Presentar el resultado de MLL1/2 como «enriquecido en el extremo de menor ruido», con significancia clara en la réplica 1 y la misma dirección en la réplica 2.
- Apoyarlo en el análisis por miembro del núcleo COMPASS, que es consistente en ambas réplicas y sobrevive al ajuste.
- Mencionar que la versión ajustada es conservadora y que no controla por amplitud de unión.

## Criterio de consistencia para reportar el GSEA

- **Confirmo el número.** En el GSEA actual (ranking por AUC), MLL1/2 complex da p ajustado 0.0046 / 0.081. NES −2.47 / −2.01. Se reproduce exactamente con este script.
- **Criterio recomendado:** p ajustado < 0.1 en ambas réplicas y el mismo signo de NES.
  - 0.1 es el `pvalueCutoff` que ya usa el script y el corte habitual de FDR en GSEA. El GSEA no combina réplicas por sí solo, así que conviene exigir el corte en cada una.
  - Con este criterio, MLL1/2 en el GSEA original cumple (0.0046 / 0.081), pero ninguna de las versiones por coeficiente cumple.
- **Criterio más estricto:** p ajustado < 0.05 en ambas réplicas. Lo cumplen solo los términos genéricos, como nuclear protein-containing complex y catalytic complex.
- **Lo que no conviene hacer:** reportar el p de una sola réplica, o presentar 0.081 como significativo a 0.05.
- **Redacción posible:** «enriquecido en el extremo de menor ruido (NES −2.47 / −2.01; p ajustado 0.0046 / 0.081, umbral 0.1)».

## Archivos

- `p1_coeficientes_por_factor.tsv`: coeficiente, IC 95% y p corregido de `unido` en los tres modelos, más la AUC original, por factor y réplica.
- `p1_covariables_correlacion.tsv`: correlaciones de las covariables con el ruido y entre sí.
- `p1_crudo_vs_ajustado.jpg`: dispersión del efecto crudo contra los dos ajustados, con los miembros de MLL1/2 en rojo.
- `p2_conteo_consistentes.tsv`, `p2_reduccion_efecto.tsv`: conteos y reducción del efecto.
- `p3_gsea_comparacion.tsv`: NES y p ajustado lado a lado para los términos pedidos y todos los que dan p ajustado < 0.1 en algún ranking.
- `p3_gsea_todos_los_terminos.tsv`: todos los términos de los cuatro rankings.
- `p3_gsea_nes_comparacion.jpg`: NES de los términos pedidos en los cuatro rankings.
- `p3_gsea_mll1_2_ajustado_Rep1.jpg`, `p3_gsea_mll1_2_ajustado_Rep2.jpg`: curvas de enriquecimiento de MLL1/2 sobre el ranking ajustado.
- `p4_mll1_2_miembros.tsv`, `p4_mll1_2_miembros.jpg`: efecto crudo y ajustados de cada miembro, con n de promotores unidos y si está en el leading edge.

---

# Puntos 5 y 6

Script: `R/48_exploratorio_remap_ruido_ajustado/remap_ocupacion_actividad.R`, que corre en ~1 min. Las salidas tienen prefijo p5_ y p6_ y están en esta misma carpeta.

El pedido anunciaba «tres puntos nuevos», pero solo llegaron dos (5 y 6). Si había un tercero, no lo recibí.

## p5. Ocupación cuantitativa y ruido

**Conteo.** Para cada promotor conté cuántos factores distintos tienen al menos un pico que lo solapa, en dos versiones:

- con el catálogo no redundante completo de ReMap (1205 factores; no hay marcas de histona en `remap_tf_hits.tsv`);
- con los 982 factores del análisis principal.

**Las dos versiones dan lo mismo.**

- La correlación de Spearman entre ellas es 1.000 / 1.000. La mediana es 431 / 412 factores con el catálogo completo y 430 / 411 con los 982.
- Los 223 factores que quedan afuera son los que unen 100 promotores o menos en alguna réplica, así que casi no cambian el conteo.
- Todo lo que sigue usa el catálogo completo; los resultados con los 982 coinciden hasta el tercer decimal (`p5_ocupacion_versiones.tsv`).

**Relación con la binaria «sin módulo cis-regulatorio».** Es exactamente «0 factores» (verificado: sale del mismo bed de ReMap).

- Solo el 1.2% / 1.3% de los promotores está en esa categoría: 115 promotores en la Rep 1.
- La binaria es, en la práctica, el extremo inferior de una variable continua donde el promotor típico está unido por más de 400 factores.

**Relación con el ruido.**

- Spearman crudo entre el conteo y el rango de varianza: −0.14 / −0.05.
- Medias por bin, Rep 1 (`p5_ruido_por_bin_ocupacion.jpg`; bins de 0 factores y deciles de los unidos):

  | Bin (factores unidos) | Rango de varianza medio |
  |---|---|
  | 0 | 63.1 |
  | 1–53 | 59.8 |
  | 54–153 | 57.7 |
  | 153–252 | 53.5 |
  | 252–352 | 50.3 |
  | más de ~350 | 46–48 (meseta) |

- La Rep 2 tiene la misma forma, más plana.
- **La relación es gradual en el tercio inferior (hasta ~300 factores) y después se satura.** En los dos o tres deciles superiores hay incluso una leve suba.

**Primer control: actividad endógena, accesibilidad e isla CpG.** El conteo está muy correlacionado con las tres covariables. Spearman con la Rep 1:

| Covariable | Correlación con el conteo |
|---|---|
| Accesibilidad | 0.79 |
| Actividad en HEK293 | 0.65 |
| Isla CpG | 0.58 |

Lo que queda después de ajustar:

- La Spearman parcial (sobre rangos) cae a 0.02 / 0.003. **En el grueso de la distribución, el gradiente se explica por las covariables.**
- En las medias por bin ajustadas, lo único que sobrevive es el extremo inferior:
  - 0 factores: 57.3 / 58.4.
  - 1–53 factores: 52.4 / 52.5.
  - A partir de ~100 factores, todo queda en 49–52, alrededor de la media de 50.5.

**Segundo control: conteo frente a binaria** (`p5_modelos_binaria_vs_conteo.tsv`). Modelos con `log2(n + 1)`, con la binaria y con ambas:

- Sin covariables:
  - El conteo explica más que la binaria: R² 0.023 contra 0.002 en la Rep 1, 0.005 contra 0.001 en la Rep 2.
  - Con las dos en el modelo, el conteo sigue siendo muy significativo (p 3e-45 / 8e-11). La binaria pierde el efecto: en la Rep 1 cambia de signo; en la Rep 2 queda en p 0.86.
- Con covariables:
  - El conteo agrega información sobre la binaria: −0.75 [−1.26, −0.23] / −0.67 [−1.13, −0.21] posiciones de rango por cada duplicación del número de factores, p 0.004 / 0.004.
  - La binaria no agrega nada sobre el conteo: p 0.55 / 0.22.
  - Los efectos son chicos: R² total 0.035 / 0.010.
  - La discrepancia con la Spearman parcial (≈0) indica que el efecto ajustado está en el extremo de pocos factores, que la escala log2 pondera más.

**Lectura de p5.**

- El conteo subsume a la anotación binaria. Justifica reemplazar «sin módulo cis-regulatorio» por el número de factores unidos, o al menos complementarla.
- Pero la señal es sobre todo de promotores con muy poca o ninguna ocupación. Por encima de ~100–300 factores se satura.
- Ajustada por actividad endógena, accesibilidad e isla CpG, la señal es chica.

## p6. Efecto sobre la actividad frente a efecto sobre el ruido

**Datos.** Para 982 factores usé:

- el efecto sobre la actividad: HL de la unión sobre la actividad media (`remap_tf_activity_wilcoxon.tsv`, remap_tf_activity_gsea);
- el efecto sobre el ruido: AUC (`remap_tf_noise_auc.tsv`).

Usé percentiles dentro de cada réplica en vez de estandarizar. Así quedan en la misma escala y no pesan los factores extremos.

**Correlaciones** (`p6_correlaciones.tsv`):

| | Rep 1 | Rep 2 |
|---|---|---|
| Spearman global | −0.70 | −0.30 |
| Spearman con n unidos, actividad | 0.62 | 0.58 |
| Spearman con n unidos, ruido | −0.84 | −0.70 |
| Spearman parcial controlando por n unidos | −0.42 | +0.18 |

- Los factores que más suben la actividad son los que más bajan el ruido.
- **Buena parte de la correlación es cobertura.** En la Rep 2 la asociación desaparece y hasta cambia de signo al controlar por n unidos.

**Factores que se despegan de la tendencia** (`p6_residuos_por_factor.tsv`). Ajusté un modelo robusto (`MASS::rlm`) del percentil de ruido sobre el de actividad y estandaricé el residuo por su MAD.

- **Criterio pedido, |z| ≥ 2 en ambas réplicas con el mismo signo: ningún factor lo cumple.**
- Los 15 de cada extremo por residuo medio, con el mismo signo en ambas réplicas, están en `p6_top_desviados_consistentes.tsv`:
  - **Menos ruido de lo esperado:** ZNF846, HMGB1, SREBP2, TRIM25, MAF1, BCOR, ZBTB7A, KMT2A, ZFX, YY1, TFAP2C, KDM4C, CTCFL, ZNF783 y GLIS2. Mediana de cobertura: percentil 0.85; 11 de 15 unen más de 4000 promotores.
  - **Más ruido de lo esperado:** GTF2A2, MLXIP, CHAF1B, SNAPC1, IRF9, TAF2, PAX8, PAX6, PALB2, E2F7, DDX20, MYF5, REL, MAFB y SOX11. Mediana de cobertura: percentil 0.18; todos unen menos de 2500 promotores.

**Residuo controlando también por cobertura** (`pct_ruido ~ pct_actividad + pct_n`; `p6_top_desviados_consistentes_ajustado_cobertura.tsv`).

- 23 factores quedan con «menos ruido de lo esperado» (|z| ≥ 2 en ambas réplicas) y 1 con «más ruido».
- Pero los de menos ruido son casi todos ZNF con muy poca cobertura (110–800 promotores; mediana del percentil de cobertura 0.07): ZNF846, ZMYM2, YBX3, ZNF331, ZNF266, ZNF322 y otros.
- Es lo esperable si la relación ruido–cobertura no es lineal y si las estimaciones con pocos promotores son imprecisas.

**Respuesta a la pregunta de cobertura.** Sí: en las dos versiones, los factores que se despegan son de cobertura extrema.

- Sin controlar por cobertura: los de menos ruido de lo esperado tienen cobertura alta y los de más ruido, cobertura baja.
- Controlando por cobertura: aparecen los ZNF de cobertura mínima.
- En ninguna de las dos versiones aparece un grupo de cobertura intermedia que se despegue. **Tal como está, el análisis de desviados no es interpretable.**

**Grupos destacados.** Mediana del residuo z medio entre réplicas:

- **MLL1/2 (18 miembros con ambos efectos):**
  - Sin controlar por cobertura: −0.76.
  - Controlando por cobertura: −0.39.
  - Los del núcleo quedan sistemáticamente por debajo de la tendencia, con menos ruido de lo que predice su efecto sobre la actividad: KMT2A −2.0 / −1.2, RBBP5 −1.7 / −1.3, WDR5 −1.2 / −1.4, ASH2L −1.5 / −0.9 y KMT2B −1.0 / −1.0.
  - Ninguno llega a |z| ≥ 2 en ambas réplicas, y casi todos tienen cobertura alta.
- **Polycomb:**
  - Grupo: RING1, JARID2, CBX8 y PHF19, los que remap_tf_activity_gsea marca en actividad baja, más los componentes canónicos de PRC1/PRC2 presentes.
  - Mediana z: −0.01 sin controlar por cobertura y 0.11 controlando.
  - Están sobre la tendencia. Los de R2 están en el extremo de menor actividad y mayor ruido (percentil de actividad 0.00–0.01; de ruido 0.97–1.00 en la Rep 1), que es justo donde los ubica la tendencia general. Ruido alto y actividad baja no son dos señales independientes.
  - Además, tienen cobertura baja (195–3000 promotores).
- **NF-YA:** el factor con mayor efecto sobre la actividad (percentil 0.99) y de los que más bajan el ruido (0.03 / 0.22). Residuo −0.48 / −0.36, sobre la tendencia.

## ¿Qué entra en la figura R3.7?

1. **Control por actividad endógena y accesibilidad (puntos 1–4): candidato a panel.**
   - Panel recomendado: el forest plot por miembro de MLL1/2 (`p4_mll1_2_miembros.jpg`), crudo frente a ajustado.
   - Es la forma más honesta de mostrar que el resultado de MLL1/2 sobrevive al control.
   - El GSEA ajustado va como texto, con los números Rep 1 / Rep 2 y la advertencia de que en la Rep 2 no llega a 0.1.
2. **Ocupación cuantitativa (p5): texto, no panel.**
   - El resultado útil es metodológico: el conteo subsume a la binaria y justifica reemplazarla en el texto del panorama.
   - Ajustado, el efecto es chico y concentrado en los promotores con muy pocos factores (~1% sin ninguno).
   - Si se quiere mostrar algo, el gráfico de bins encaja mejor como complemento de R3.6 que en R3.7.
3. **Actividad frente a ruido (p6): descartar como panel.**
   - La correlación global se explica en gran parte por la cobertura; en la Rep 2 cambia de signo al controlarla.
   - No hay desviados consistentes con el criterio pedido.
   - Los candidatos que aparecen son todos de cobertura extrema.
   - A lo sumo, una frase: que NF-YA y Polycomb están sobre la tendencia, y que el núcleo COMPASS queda levemente por debajo en ambas réplicas sin llegar al umbral.

## Archivos (puntos 5 y 6)

- `p5_ocupacion_versiones.tsv`: comparación de las dos versiones del conteo, correlaciones con el ruido y con las covariables.
- `p5_rho_parcial_covariables.tsv`: Spearman parcial del conteo con el ruido, controlando por las covariables.
- `p5_ruido_por_bin_ocupacion.tsv`, `p5_ruido_por_bin_ocupacion.jpg`: rango de varianza medio por bin, crudo y ajustado, con IC 95%.
- `p5_modelos_binaria_vs_conteo.tsv`: modelos con la binaria, el conteo y ambas, con y sin covariables, y F parciales.
- `p6_correlaciones.tsv`: Spearman global, parcial y con la cobertura, y el ajuste robusto.
- `p6_residuos_por_factor.tsv`: percentiles, residuos con y sin control por cobertura, y grupo, para los 982 factores.
- `p6_desviados_y_destacados.tsv`: los grupos destacados y los desviados.
- `p6_top_desviados_consistentes.tsv`, `p6_top_desviados_consistentes_ajustado_cobertura.tsv`: los 15 de cada extremo.
- `p6_actividad_vs_ruido.jpg`: dispersión por réplica, coloreada por cobertura, con MLL1/2, Polycomb y NF-YA destacados.
- `p6_residuo_vs_cobertura.jpg`: residuo frente a cobertura.

# Punto 7: control por amplitud de unión (cobertura)

Script: `R/48_exploratorio_remap_ruido_ajustado/remap_ruido_cobertura.R`. Salidas con prefijo `p7`.

## Método

- **Coeficiente:** el de `unido` del modelo ajustado de p1 (`var_rank_sw ~ unido + log10(hek_tpm + 1) + mean_dnase`). Para los modelos nuevos, mismo modelo con un agregado o una restricción.
- **p7.1 Ventana de cobertura:** mismo procedimiento que `nfy_especificidad`.
  - Referencia: los otros 981 factores (se excluye al propio).
  - Ventana: factores que unen un n de promotores a ±25% del focal.
  - Puesto 1 = coeficiente más negativo (mayor reducción del ruido). Percentil = % de factores comparables con coeficiente ≤ al del focal, así que bajo = más reducción.
  - También con una ventana que excluye a todos los miembros de MLL1/2: cambia el percentil en menos de 2 puntos, así que no se reporta aparte (está en la tabla).
- **p7.3 Ocupación total como covariable:** modelo ajustado + log2(n + 1), con n = factores distintos del catálogo ReMap completo (1205) con pico sobre el promotor. Incluye al propio factor; sobre ~400 factores en escala log2 su aporte es despreciable.
- **p7.4 Referencia restringida:** modelo ajustado excluyendo los promotores del decil inferior de ocupación (≤ 46 / ≤ 41 factores). Como versión más estricta, también sin el cuartil inferior (≤ 204 / ≤ 185).
- **GSEA de MLL1/2 complex** sobre los rankings nuevos, con los mismos parámetros que p3.

## Resultados

### p7.1 Posición del núcleo COMPASS (modelo ajustado)

Puesto / n de factores comparados, y percentil. Rep 1 / Rep 2.

| Factor | n unidos | Universo | Ventana ±25% |
|---|---|---|---|
| KMT2B | 6679 / 7969 | 33/982 (p3) / 25/982 (p2) | 33/350 (p9) / 21/345 (p6) |
| WDR5 | 5590 / 6717 | 103/982 (p10) / 66/982 (p7) | 50/369 (p13) / 29/353 (p8) |
| RBBP5 | 7041 / 8478 | 76/982 (p8) / 49/982 (p5) | 70/313 (p22) / 41/300 (p13) |
| KMT2A | 7374 / 8912 | 115/982 (p12) / 160/982 (p16) | 98/282 (p35) / 131/264 (p49) |
| ASH2L | 7173 / 8587 | 141/982 (p14) / 218/982 (p22) | 117/302 (p39) / 173/291 (p59) |

- La ventana de cobertura de estos factores es la región de efecto más fuerte de todo el catálogo: la mediana del coeficiente en la ventana es −2.0 a −3.1, contra −1.08 / −1.05 en los 982.
- **KMT2B y WDR5** quedan en el 10–15% más fuerte entre factores de cobertura comparable en ambas réplicas. **RBBP5** queda en p22 / p13.
- **KMT2A y ASH2L** quedan en la mediana de su ventana (p35 / p49 y p39 / p59). Su efecto es el que tiene cualquier factor que une esa cantidad de promotores.
- Los otros 14 miembros de MLL1/2 están en `p7_1_ventana_cobertura_mll.tsv`, para todos los modelos.

### p7.2 Figura

`p7_2_efecto_vs_cobertura_compass.jpg`: coeficiente ajustado frente a n unidos (log10), una columna por miembro del núcleo y una fila por réplica, con la ventana sombreada. Mismo formato que el panel de cobertura de `nfy_especificidad`.

La forma de la nube es la explicación del confusor. El coeficiente está centrado en 0 hasta ~3000 promotores y cae de forma pronunciada por encima de ~4000, que es donde están todos los miembros del núcleo.

### p7.3 Ocupación total como covariable

| | Ajustado | Ajustado + log2(ocupación) |
|---|---|---|
| Mediana de los 982 | −1.08 / −1.05 | 0.12 / 0.07 |
| Mediana de MLL1/2 (19) | −3.34 / −2.83 | −1.11 / −0.71 |
| Mediana del núcleo COMPASS | −3.70 / −4.17 | −1.96 / −2.45 |
| Consistentes de menor ruido (p corregido < 0.05 en ambas, mismo signo) | 276 | 27 |
| Spearman n unidos – coeficiente | −0.63 / −0.66 | −0.36 / −0.42 |
| GSEA MLL1/2 (NES; p ajustado) | −1.76; 0.018 / −1.46; 0.18 | −1.86; 0.058 / −1.55; 0.43 |

- La ocupación absorbe casi todo el efecto «promedio» de la unión: la mediana de los 982 pasa a 0. El coeficiente de log2(ocupación) es −1.29 / −1.29 por duplicación.
- Por miembro del núcleo, coeficiente y p corregido, Rep 1 / Rep 2:
  - **WDR5:** −2.57 (0.007) / −3.07 (0.0002). Ventana p8 / p3.
  - **KMT2B:** −3.17 (0.021) / −3.04 (0.011). Ventana p10 / p5.
  - **RBBP5:** −1.96 (0.19) / −2.45 (0.042). Ventana p26 / p13.
  - **KMT2A:** −1.37 (0.42) / −0.71 (0.66). Ventana p38 / p58.
  - **ASH2L:** −0.97 (0.59) / −0.30 (0.88). Ventana p47 / p69.
- Queda dependencia residual de la cobertura (Spearman −0.36 / −0.42): la covariable a nivel promotor no absorbe todo el sesgo a nivel factor.

### p7.4 Referencia restringida

| | n promotores | Grupo no unido del núcleo |
|---|---|---|
| Completo | 9300 / 11400 | 1926–3710 / 2488–4683 |
| Sin decil inferior (ocupación > 46 / > 41) | 8358 / 10251 | 1078–2801 / 1442–3571 |
| Sin cuartil inferior (ocupación > 204 / > 185) | 6973 / 8544 | 542–1869 / 715–2368 |

Coeficiente y p corregido, Rep 1 / Rep 2:

| Factor | Sin decil inferior | Sin cuartil inferior |
|---|---|---|
| WDR5 | −2.97 (0.0004) / −3.27 (< 0.0001) | −2.26 (0.060) / −2.22 (0.037) |
| KMT2B | −3.97 (0.0003) / −3.71 (0.0001) | −2.05 (0.38) / −2.86 (0.087) |
| RBBP5 | −2.86 (0.015) / −2.99 (0.003) | −1.19 (0.63) / 0.26 (0.94) |
| KMT2A | −2.46 (0.060) / −1.93 (0.092) | −2.57 (0.23) / −1.66 (0.42) |
| ASH2L | −2.06 (0.11) / −0.94 (0.45) | −1.84 (0.43) / −0.47 (0.89) |

- **Sin el decil inferior:** WDR5, KMT2B y RBBP5 se mantienen en ambas réplicas, y sus percentiles de ventana casi no cambian (p11 / p6, p9 / p5, p25 / p17). KMT2A y ASH2L pierden significación.
  - GSEA MLL1/2: NES −1.70 / −1.48; p ajustado 0.043 / 0.36.
  - La mediana de los 982 baja a −0.48 / −0.53 y quedan 100 consistentes.
- **Sin el cuartil inferior:** el grupo no unido de KMT2A, ASH2L y RBBP5 queda en 540–930 promotores, y se pierde potencia en todo el catálogo (mediana de los 982: 0.31 / −0.06; 10 consistentes).
  - Solo WDR5 se mantiene cerca del umbral en ambas réplicas.
  - GSEA MLL1/2: p ajustado 0.26 / 0.98.
  - Con este corte no se puede separar «el efecto desaparece» de «no hay potencia».

## Lectura: ¿el resultado de MLL1/2 se sostiene entre factores de cobertura comparable?

**A nivel complejo, no: se diluye.**

- Buena parte de la señal de MLL1/2 en el GSEA es cobertura. Sus miembros están entre los factores que más promotores unen, y en esa franja casi cualquier factor da un efecto fuerte.
- Con la ocupación como covariable o sin el decil inferior, el GSEA de MLL1/2 queda en p ajustado 0.06 / 0.43 y 0.04 / 0.36.
  - La Rep 1 sigue cerca o debajo de 0.05.
  - La Rep 2, que ya no pasaba 0.1 con el modelo ajustado (0.18), se aleja.
- Con el criterio recomendado (p ajustado < 0.1 en ambas réplicas con el mismo signo), MLL1/2 complex no pasa en ningún ranking ajustado. Solo pasa con la AUC original.

**A nivel miembro, se sostiene para una parte del núcleo, con fuerza moderada.**

- **WDR5 y KMT2B:**
  - Están en el 3–15% más fuerte entre factores de cobertura comparable, en ambas réplicas y con todos los modelos.
  - Siguen significativos con la ocupación como covariable y sin el decil inferior.
  - WDR5 incluso sin el cuartil inferior.
  - En el modelo ajustado su efecto es −3.7 a −5.0 puntos de rango sobre 100, contra −2.0 a −2.9 de la mediana de su ventana: el exceso propio es de ~1.5 a 2.5 puntos. Ajustado por ocupación, es −2.6 a −3.2. Es claro pero modesto.
- **RBBP5:** intermedio. Se sostiene sin el decil inferior y con la ocupación en la Rep 2; no en la Rep 1.
- **KMT2A (MLL1) y ASH2L:** no se distinguen de los factores de su misma cobertura en ninguna versión.
  - Su efecto en el modelo ajustado era el esperable por la cantidad de promotores que unen.
  - Este es el punto más importante para la redacción: el componente catalítico de MLL1 no muestra un efecto propio.

## Implicancias

- **Para la sección:**
  - No conviene redactar «los promotores unidos por el complejo MLL1/2 tienen menos ruido» como un resultado específico del complejo.
  - Una redacción defendible: «el efecto de la unión sobre el ruido depende fuertemente de cuántos promotores une cada factor; entre factores de cobertura comparable, WDR5 y KMT2B (MLL2) se mantienen entre los de mayor reducción del ruido, pero no KMT2A ni ASH2L».
  - El GSEA de MLL1/2 con la AUC original puede mencionarse como el hallazgo que motivó el control, no como el resultado final.
- **WDR5 no es exclusivo de MLL1/2:** también está en SET1A/B, en NSL (con KAT8, MCRS1 y PHF20) y en ATAC.
  - PHF20, de NSL, también se sostiene en la Rep 1 con la ocupación como covariable y sin el decil inferior.
  - El dato de WDR5 no se puede atribuir a MLL1/2 en particular.
- **Para R3.8 (histonas):**
  - La motivación «MLL1/2 → H3K4me3 → menor ruido» queda debilitada. El único miembro catalítico que se sostiene es KMT2B, y el núcleo como complejo no se sostiene.
  - Si R3.8 se mantiene, conviene motivarla de forma independiente: las marcas de histona como contexto cromatínico endógeno, igual que la accesibilidad. Como mucho, mencionar KMT2B y WDR5 como indicio, no como fundamento.
  - Si el análisis de histonas muestra que H3K4me3 se asocia a menor ruido, conviene aplicar el mismo control por amplitud: H3K4me3 marca casi todos los promotores activos, así que tiene exactamente el mismo problema de cobertura.
- **Metodológico, para todo R3.7:** el número de promotores unidos debería ser parte del modelo o del criterio de comparación por defecto. Sin ese control, cualquier GSEA sobre el ranking de efectos premia a los complejos de cobertura alta (maquinaria general: TAF1, MAX, E2F6, HCFC1 tienen el mismo patrón).

## Archivos (punto 7)

- `p7_1_ventana_cobertura_mll.tsv`: para los 19 miembros de MLL1/2 y los seis modelos, puesto y percentil en el universo, en la ventana ±25% y en la ventana sin miembros de MLL1/2, además de la mediana de la ventana.
- `p7_2_efecto_vs_cobertura_compass.jpg`: figura de cobertura del núcleo COMPASS (modelo ajustado).
- `p7_3_resumen_modelos.tsv`: por modelo y réplica, la mediana de los 982, de MLL1/2 y del núcleo, la Spearman n unidos – coeficiente y el conteo de consistentes.
- `p7_3_4_coeficientes_por_factor.tsv`: coeficientes de los 982 factores en los tres modelos nuevos.
- `p7_3_4_mll1_2_miembros.jpg`: forest plot de los miembros de MLL1/2 en el modelo ajustado y los tres controles, con el núcleo en negrita.
- `p7_4_n_promotores.tsv`: n de promotores y cortes de ocupación por réplica.
- `p7_gsea_mll1_2.tsv`: NES y p ajustado de MLL1/2 complex en los siete rankings (cuatro de p3 y tres nuevos).
