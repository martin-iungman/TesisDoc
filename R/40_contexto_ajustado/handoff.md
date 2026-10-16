# Handoff — características de contexto endógeno: ajustado y estratificado por isla CpG (2026-10-16)

Único handoff de los dos exploratorios sobre las características de contexto
endógeno de `summary_endo`:

- **Modelo ajustado (análisis principal):** `R/40_contexto_ajustado/contexto_ajustado.R`,
  con salidas en `figures/exploratorio_contexto_ajustado/`. Incluye el
  modelo con todas las características de contexto juntas y el ajuste de
  cada característica solo por secuencia.
- **Estratificado por isla CpG (referencia):** `R/39_contexto_cgi/contexto_cgi.R`,
  con salidas en `figures/exploratorio_contexto_cgi/`.

Cada script corre en menos de 1 min desde la raíz del repo. Nada está
commiteado.

## Estado general
- **Figura de tesis R2.13** (slug `summary_features_endogeno`):
  `summary_endo.jpg` pasó a ser figura propia, sin cambios de contenido.
  - `R/15_summary_features/summary_features.R` lo guarda vía `fig_dir("summary_features_endogeno")`.
  - El JPG se movió con `git mv` desde `figures/R2.11_summary_features_secuencia/`.
  - `summary_seq.jpg` y `activity_summary.tsv` no se tocaron.
- Ninguno de los dos exploratorios es figura de la tesis: no usan `fig_dir()`
  y no tienen número.

---

## 1. Modelo multivariable ajustado (principal)

### Diseño
- **Respuesta:** percentil de actividad dentro de cada réplica, calculado sobre
  todos los promotores high-conf. Es el mismo que en `R/38_aditividad`. Los
  coeficientes están en puntos de percentil.
- **Modelo:** un `lm` por réplica, con las características de contexto como
  efectos principales, todas a la vez, más:
  - **razón CpG o/e del fragmento** (`cpg_oe`, continua);
  - **contenido de G+C** (`g_c`, continuo).

  Las dos covariables están estandarizadas: el coeficiente se lee por 1 DE.
  Reemplazan a la estratificación por isla CpG.
- **Efecto crudo:** `lm(percentil ~ característica)` en la misma muestra y con
  la misma codificación que el ajustado.
- **Ajuste por comparaciones múltiples:** BH dentro de cada réplica y modelo.
  La familia son las 13 características; covariables e indicadores quedan
  afuera.
- **"Sobrevive":** p BH < 0.05 en ambas réplicas e IC 95% del mismo signo.

### Características excluidas por redundancia (15 → 13)

| Par redundante | Se deja | Se saca | Por qué |
|---|---|---|---|
| Alta / baja especificidad tisular (terciles extremos del mismo índice de Gini, phi −0.50) | **Baja especificidad tisular** | Alta especificidad tisular | Menos colineal con CpG o/e (r 0.36 vs. −0.48). La referencia pasa a ser el tercil medio + el alto. |
| Promotores angostos / anchos (terciles extremos del mismo ancho interquantil, phi −0.50) | **Promotores anchos** | Promotores angostos | Menos colineal con CpG o/e (0.24 vs. −0.33) y con especificidad tisular. La referencia pasa a ser angostos + intermedios. |

El resto se mantiene:
- Las tres ventanas de conservación son regiones distintas (phi 0.47–0.49
  entre ventanas vecinas).
- HEK293 y DNase-seq miden cosas distintas (phi 0.36).
- "Insertado en humanos" y "Sin actividad en ratón" son categorías excluyentes
  del mismo turnover, no redundantes.

**VIF generalizado:** todos ≤ 2.2 (el máximo es `cpg_oe`), así que no hay
colinealidad problemática (`ajustado_vif.tsv`).

### Datos faltantes
Los NA no son al azar. Borrarlos deja casi solo promotores con expresión
endógena: en casos completos, "No detectado" baja de 8% a 0.2%.

- **Especificidad tisular y forma** no están definidas para los no detectados
  en FANTOM5. El Gini es NA exactamente en los 821 "No detectado". Ahí se
  codifican FALSE, porque "No detectado (FANTOM5)" ya distingue a esos
  promotores. Sin esto, "No detectado" queda constante y no se puede estimar.
- **Modelo principal:** los NA restantes se codifican FALSE y se agrega un
  indicador de faltante por variable:
  - turnover sin anotación de Young et al.: 24%;
  - orientación sin señal PRO-cap suficiente: 22%;
  - forma NA en promotores detectados: 6%.

  n = 9.814 (Rep 1) y 11.604 (Rep 2).
- **Sensibilidad:** el mismo modelo en casos completos, sin indicadores.
  n = 6.464 y 7.415.

### Resultados: crudo vs. ajustado (modelo principal, Rep 1 / Rep 2, puntos de percentil)

| Característica | Prev. | Crudo | Ajustado [IC 95%] | p BH ajust. | Sobrevive ajust. | Casos completos |
|---|---|---|---|---|---|---|
| Alta actividad en HEK293 | 29% | +18.0 / +18.1 | +10.7 [9.3, 12.0] / +10.0 [8.7, 11.3] | 4e-51 / 5e-51 | sí | sí |
| Alta accesibilidad (DNase-seq) | 33% | +17.9 / +19.4 | +10.6 [9.3, 11.9] / +12.3 [11.0, 13.5] | 1e-52 / 2e-80 | sí | sí |
| Promotor unidireccional | 15% | +7.8 / +7.9 | +2.1 [0.5, 3.6] / +1.8 [0.4, 3.3] | 0.016 / 0.023 | sí | sí |
| Sin actividad en ratón | 6% | −9.3 / −7.9 | −4.6 [−6.9, −2.3] / −2.6 [−4.7, −0.6] | 3e-4 / 0.023 | sí | sí |
| Alta conservación (−150 a −235) | 33% | −2.7 / −2.4 | −2.1 [−3.4, −0.8] / −2.0 [−3.2, −0.8] | 0.005 / 0.004 | sí | sí |
| Sin módulo cis-regulatorio | 1.2% | +2.1 / −2.0 (n.s.) | +13.8 [8.6, 18.9] / +8.5 [3.9, 13.1] | 8e-7 / 0.001 | sí (ver ⚠) | no (n = 3) |
| No detectado (FANTOM5) | 8% | −12.6 / −11.1 | −3.6 [−6.0, −1.2] / −1.0 [−3.1, 1.2] | 0.008 / 0.41 | no (solo Rep 1) | no (0.2% de prevalencia, IC muy ancho) |
| Baja especificidad tisular | 31% | +11.6 / +12.5 | −1.3 [−2.7, 0.1] / −1.0 [−2.3, 0.3] | 0.10 / 0.15 | no | **sí, negativo** (−1.8 / −1.9) |
| Promotores anchos | 28% | +8.2 / +9.0 | +0.7 [−0.6, 2.0] / +1.5 [0.3, 2.7] | 0.36 / 0.023 | no | no |
| Enhancers a 50kb | 70% | −2.2 / −2.5 | −0.5 / −0.8 (n.s.) | 0.43 / 0.19 | no | no |
| Alta conservación (16 a −50pb) | 33% | +0.7 / +1.0 (n.s.) | ≈ 0 | 0.98 / 0.78 | no | no |
| Alta conservación (−50 a −150) | 33% | −1.0 / −0.7 (n.s.) | +1.0 / +1.5 | 0.25 / 0.037 | no | no (Rep 2 sí) |
| Insertado en humanos | 2% | +2.4 / +3.8 (n.s.) | +1.5 / +2.9 (n.s.) | 0.47 / 0.15 | no | no |

**Covariables** (por 1 DE, `ajustado_covariables.tsv`):
- CpG o/e: +3.4 / +3.2 (p ≈ 1e-17 en ambas).
- G+C: +0.7 (p 0.03) / +0.1 (n.s.).

Ya ajustado por las características endógenas, la composición del fragmento
aporta poco, y lo que aporta lo lleva la razón CpG o/e, no el G+C.

### R² del modelo (`ajustado_r2.tsv`)

| Modelo | Rep 1 | Rep 2 |
|---|---|---|
| Completo (13 características + indicadores + CpG o/e + G+C) | **0.144** | **0.146** |
| Solo CpG o/e + G+C | 0.071 | 0.069 |
| Solo características (+ indicadores) | 0.134 | 0.139 |
| Completo, casos completos | 0.127 | 0.127 |

Las características de contexto explican casi el doble que la composición del
fragmento. Sumar CpG o/e + G+C a las características agrega ~0.01 de R².

### Lectura
1. **Robustas al ajuste: HEK293 y DNase-seq.** Bajan de ~+18 a ~+10–12
   puntos, pero siguen siendo por lejos los efectos más grandes. Se ajustan
   mutuamente: comparten parte del efecto.
2. **Sobreviven, atenuadas, en ambos modelos:**
   - "Promotor unidireccional": +8 → +2;
   - "Sin actividad en ratón": −9/−8 → −5/−3;
   - "Alta conservación (−150 a −235)": casi sin cambio, ≈ −2.
3. **Desaparecen al ajustar:** baja especificidad tisular, promotores anchos,
   enhancers a 50kb y, en Rep 2, "No detectado".
   - Su efecto crudo se explica por la covariación con HEK293, DNase-seq y
     CpG o/e.
   - Baja especificidad tisular incluso pasa a levemente negativa: −1 a −2
     puntos, significativa en casos completos.
4. ⚠ **"Sin módulo cis-regulatorio":** de ≈ 0 crudo pasa a +8 / +14 ajustado.
   Es un efecto de supresión con solo 120 / 148 promotores (1.2%).
   - Está asociada a "No detectado" (phi 0.30) y a baja expresión endógena.
     Una vez ajustado por eso, los pocos promotores sin CRM resultan más
     activos de lo esperado.
   - En casos completos quedan 3 promotores y el efecto no es estimable.
   - No lo tomaría como hallazgo sin mirarlo aparte.
5. **Comparación con el estratificado** (sección 2):
   - Coinciden en que HEK293 y DNase-seq son independientes de la isla CpG.
   - Coinciden en que unidireccional, sin actividad en ratón y conservación
     −150 a −235 se sostienen.
   - Coinciden en que enhancers a 50kb es composición.
   - Difieren en especificidad tisular, "No detectado" y forma: en el
     estratificado se sostenían dentro de los grupos de isla CpG, pero en el
     ajustado desaparecen. Lo que las absorbe es el ajuste por HEK293 y
     DNase-seq, no la composición CpG. El estratificado no ajustaba por las
     otras características.

### Ajustado solo por secuencia (`secuencia_efectos.tsv`, `secuencia_r2.tsv`, `secuencia_crudo_vs_ajustado.jpg`)

**Diseño.** Cada característica de contexto entra en un modelo propio, de a
una, junto con:
- todas las características de secuencia del summary (grupo `seq` de
  `feature_groups`, 19 en total): TATA, CCAAT, GC-box, isla CpG, YY1, CRE,
  E-box, poli(dA:dT), dinucleótido del TSS, INR fuerte, TCT, LINE/SINE/TE y
  alto G+C;
- CpG o/e y G+C continuas.

No entran las demás características de contexto. Como entran de a una, se
evalúan las 15, incluidas alta especificidad tisular y promotores angostos.
Los datos y la codificación son los mismos que en el modelo principal: n
9.814 / 11.604 e indicador de faltante de la propia característica. BH dentro
de réplica sobre las 15.

- R² de solo secuencia: **0.129 / 0.119**.
- VIF máximo: 7.4, en "CA en TSS". Las categorías de dinucleótido del TSS son
  casi excluyentes entre sí. Esto infla los coeficientes de secuencia, no los
  de contexto.

| Característica | Crudo | Ajust. secuencia [IC 95%] | ΔR² | Ajust. completo | Sobrevive (secuencia) |
|---|---|---|---|---|---|
| Alta actividad en HEK293 | +18.0 / +18.1 | +12.4 [11.1, 13.6] / +12.3 [11.1, 13.5] | 0.032 / 0.030 | +10.7 / +10.0 | sí |
| Alta accesibilidad (DNase-seq) | +17.9 / +19.4 | +11.4 [10.1, 12.6] / +13.1 [11.9, 14.2] | 0.026 / 0.034 | +10.6 / +12.3 | sí |
| Baja especificidad tisular | +11.6 / +12.5 | +6.4 [5.1, 7.6] / +6.8 [5.6, 7.9] | 0.008 / 0.009 | −1.3 / −1.0 (n.s.) | sí |
| Alta especificidad tisular | −9.2 / −10.5 | −4.1 [−5.4, −2.8] / −5.4 [−6.6, −4.2] | 0.003 / 0.006 | (excluida) | sí |
| Promotores anchos | +8.2 / +9.0 | +3.5 [2.2, 4.7] / +4.1 [3.0, 5.3] | 0.003 / 0.004 | +0.7 / +1.5 (n.s.) | sí |
| Promotores angostos | −4.0 / −5.8 | −3.4 [−4.7, −2.1] / −4.8 [−6.0, −3.6] | 0.002 / 0.005 | (excluida) | sí |
| Promotor unidireccional | +7.8 / +7.9 | +3.0 [1.5, 4.5] / +3.0 [1.6, 4.4] | 0.001 / 0.001 | +2.1 / +1.8 | sí |
| Sin actividad en ratón | −9.3 / −7.9 | −6.9 [−9.1, −4.7] / −5.2 [−7.2, −3.2] | 0.003 / 0.002 | −4.6 / −2.6 | sí |
| Alta conservación (−150 a −235) | −2.7 / −2.4 | −2.2 [−3.3, −1.0] / −1.9 [−3.0, −0.9] | 0.001 / 0.001 | −2.1 / −2.0 | sí |
| No detectado (FANTOM5) | −12.6 / −11.1 | −2.4 [−4.5, −0.3] / −1.3 [−3.2, 0.6] | <0.001 | −3.6 / −1.0 | no (solo Rep 1) |
| Enhancers a 50kb | −2.2 / −2.5 | −1.2 / −1.4 (Rep 1 n.s.) | <0.001 | n.s. | no |
| Sin módulo cis-regulatorio | +2.1 / −2.0 (n.s.) | +7.7 [2.7, 12.6] / +3.8 (n.s.) | <0.001 | +13.8 / +8.5 | no (solo Rep 1) |
| Alta conservación (16 a −50pb) | +0.7 / +1.0 | +0.9 / +1.2 (Rep 2 p 0.04) | <0.001 | n.s. | no |
| Alta conservación (−50 a −150) | n.s. | n.s. | <0.001 | n.s. | no |
| Insertado en humanos | n.s. | n.s. | <0.001 | n.s. | no |

**Lectura (crudo → secuencia → completo)**
1. **La secuencia explica parte de todos los efectos grandes, pero no todo.**
   - HEK293 y DNase-seq conservan ~2/3 del efecto crudo (+12/+13) y son los
     que más R² agregan sobre la secuencia (~0.03 cada uno).
   - Especificidad tisular y forma conservan ~40–60% y siguen claramente
     significativas: la secuencia sola no las absorbe.
2. **"No detectado (FANTOM5)" es casi todo secuencia.** Pasa de −12/−11 a
   −2/−1, y no sobrevive en Rep 2. Sus fragmentos son mayormente sin isla
   CpG (72% en Rep 1, vs. 29% del total), y eso explica casi toda su baja
   actividad.
3. **Lo que elimina especificidad tisular y forma en el modelo completo es el
   ajuste por las otras características de contexto, no la secuencia.**
   Sobre todo HEK293 y DNase-seq, que covarían con la especificidad tisular
   (phi 0.39–0.43).
   - Ajustadas por secuencia sola siguen en +6/+7 (baja especificidad) y
     +3/+4 (anchos).
   - En el completo caen a ≈ 0.
   - Es decir, la parte del efecto que no es secuencia la comparten con la
     actividad/accesibilidad endógena en HEK293.
4. **Robustas en las tres versiones:** unidireccional, sin actividad en ratón
   y conservación −150 a −235. Son efectos chicos (2–7 puntos) pero
   independientes de secuencia y del resto del contexto.
5. **"Sin módulo cis-regulatorio"** ya aparece positivo ajustando solo por
   secuencia en Rep 1 (+7.7), no en Rep 2. Con 1.2% de prevalencia sigue
   siendo inestable (ver ⚠ arriba).

### Distribución de la razón CpG o/e (`cpg_oe_distribucion.jpg` / `.tsv`)
Son los 20.851 promotores de la library. Los demás tipos de fragmento no
tienen `cpg_oe`.

**Sí, es bimodal, pero asimétrica.**
- **Modo principal:** ~0.87, el componente CpG-rico, con 74% del peso en una
  mezcla de 2 normales (media 0.86, DE 0.17).
- **Modo secundario:** ~0.25, ancho y bajo, con 26% del peso (media 0.28,
  DE 0.16).
- **Evidencia:**
  - La mezcla de 2 normales es muy superior a 1 normal: ΔBIC = 5.052.
  - D de Ashman = 3.55. Valores > 2 indican componentes bien separados.
  - La densidad kernel tiene dos modos.
- **El valle es poco profundo.** El antimodo está en ~0.44, con una densidad
  del 82% de la altura del modo menor. Es decir, hay un continuo entre los
  dos grupos, no dos poblaciones discretas.
- **Corte de isla CpG.** El corte de 0.6 cae en el hombro del componente
  principal, no en el valle (0.44).
  - 69% de los fragmentos tiene o/e > 0.6.
  - 67% cumple `CGI_frag`, que además exige G+C ≥ 50%.
- **Pico en 0:** 358 fragmentos (1.7%) no tienen ningún CpG.

---

## 2. Estratificado por isla CpG (referencia)

### Diseño
- **Efecto:** HL + IC 95% de `mean ~ característica`, cada réplica por
  separado. Es el mismo criterio que la figura R2.13.
- **p-valor:** Wilcoxon estratificado por réplica.
- **Ajuste:** BH dentro de cada grupo (Global / Con isla CpG / Sin isla CpG),
  sobre las 15 características.
- **Grupos** (`CGI_frag`, high-conf):
  - con isla: 6.990 / 8.056;
  - sin isla: 2.833 / 3.557.
- **Cortes:** los terciles y la mediana son los globales de
  `build_tidy_features()`.
- **Control:** el grupo Global reproduce exactamente R2.13: mismos HL y mismas
  características significativas.

### Salidas
- `contexto_cgi_efectos.tsv`: formato largo.
- `contexto_cgi_resumen.tsv`: una fila por característica.
- `contexto_cgi_forest.jpg`: los tres efectos por característica, con
  facetas Rep 1 / Rep 2.

### Resultados (HL Rep 1 / Rep 2)

| Característica | Prev. con isla | Prev. sin isla | Global | Con isla | Sin isla | Sobrevive en |
|---|---|---|---|---|---|---|
| Alta actividad en HEK293 | 38% | 7% | +0.33 / +0.45 | +0.29 / +0.40 | +0.20 / +0.27 | ambos |
| Alta accesibilidad (DNase-seq) | 45% | 6% | +0.31 / +0.46 | +0.26 / +0.40 | +0.30 / +0.41 | ambos |
| Baja especificidad tisular | 42% | 8% | +0.18 / +0.27 | +0.13 / +0.19 | +0.10 / +0.16 | ambos |
| Alta especificidad tisular | 21% | 68% | −0.18 / −0.26 | −0.14 / −0.22 | −0.04 / −0.09 | ambos |
| No detectado (FANTOM5) | 3% | 21% | −0.18 / −0.21 | −0.17 / −0.21 | −0.05 / −0.03 | ambos |
| Promotores anchos | 38% | 15% | +0.10 / +0.17 | +0.07 / +0.12 | −0.00 / +0.07 | solo con isla |
| Promotor unidireccional | 19% | 26% | +0.08 / +0.11 | +0.12 / +0.16 | +0.01 / +0.03 | solo con isla |
| Promotores angostos | 26% | 60% | −0.12 / −0.20 | −0.08 / −0.14 | −0.03 / −0.07 | solo con isla |
| Sin actividad en ratón | 5% | 8% | −0.13 / −0.15 | −0.18 / −0.21 | −0.02 / −0.03 | solo con isla |
| Alta conservación (−150 a −235) | 32% | 38% | −0.04 / −0.05 | −0.03 / −0.04 | −0.02 / −0.01 | solo con isla |
| Alta conservación (16 a −50pb) | 32% | 36% | +0.01 / +0.02 (n.s.) | +0.05 / +0.07 | −0.03 / −0.02 (n.s.) | solo con isla |
| Enhancers a 50kb | 68% | 75% | −0.04 / −0.05 | −0.02 / −0.03 (n.s.) | −0.03 / −0.03 (n.s.) | ninguno |
| Alta conservación (−50 a −150) | 31% | 39% | n.s. | n.s. | n.s. | ninguno |
| Insertado en humanos | 2% | 4% | n.s. | n.s. | n.s. | ninguno |
| Sin módulo cis-regulatorio | 0.7% | 2.6% | n.s. | n.s. | n.s. | ninguno |

Las prevalencias son de Rep 1 (Rep 2 difiere menos de 2 puntos).

### Lectura
- **Ningún efecto se invierte** de signo dentro de un grupo donde sobrevive.
- **Cinco efectos son independientes de la isla:** HEK293, DNase-seq, las dos
  de especificidad tisular y "No detectado". Se atenúan sin isla, salvo
  DNase-seq.
- **Seis efectos solo se sostienen con isla:** shape, unidireccional, sin
  actividad en ratón y conservación −150 a −235. A esto se suma conservación
  16 a −50pb, que queda enmascarada en el global.
  - "Solo con isla" puede ser falta de poder: el grupo sin isla tiene menos de
    la mitad de promotores. Angostos sin isla es limítrofe.
- **Enhancers a 50kb** pierde la significancia en ambos estratos. Sugiere
  composición: hay más enhancers cerca de promotores sin isla, que son menos
  activos.

## Pendientes / ideas
- **"Sin módulo cis-regulatorio":** mirarlo aparte. Comparar actividad con y
  sin CRM dentro de promotores detectados y con HEK293/DNase similares.
- **Especificidad tisular:** si se quiere decir que su efecto es "mediado"
  por expresión endógena / accesibilidad, conviene un modelo explícito de
  mediación o, al menos, la secuencia de modelos anidados.
- **Interacción con CpG o/e:** el ajustado supone efectos aditivos. Para
  testear lo que sugería el estratificado (efectos solo con isla), habría que
  agregar interacciones característica × CpG o/e.
