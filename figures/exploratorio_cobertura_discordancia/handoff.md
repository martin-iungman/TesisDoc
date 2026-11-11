# Cobertura de ReMap en la discordancia reportero - endógeno

Script: `R/65_exploratorio_cobertura_discordancia/cobertura_discordancia.R` (~1.5 min). Es el diagnóstico de `R/51_exploratorio_remap_cobertura`, aplicado a `discordancia_remap_ora` (`R/60`). La figura de R/60 no se modificó.

- Modelo: `dif_signed ~ avg_rank + unido`, por réplica (MCO + HC1, matricial y verificado contra lm + sandwich). Además el modelo conjunto de la figura, leído de `tables/discordancia_remap_lmm.tsv`.
- 949 factores, promotores activos en HEK293 (5376 / 6298).
- Corrección por cobertura: loess robusto (span 0.5) de β sobre log10(n unidos).
  - Significancia ajustada: z = (β − esperado) / SE, BH.
  - Ranking: residuo / escala local.
  - Percentil del residuo en la ventana ±25%: no llega a FDR < 0.05 para ningún factor. El mínimo p por permutación (~0.005–0.01) no sobrevive BH sobre 949 factores.
- El ORA "original" reproduce el de R/60 (7 términos sobreestima y 54 subestima).

## Resultados (Rep 1 / Rep 2 / conjunto)

1. **Spearman β–cobertura: +0.28 / +0.24 / +0.26.** En el capítulo de ruido era −0.84 / −0.69. La asociación es real pero 3 veces más débil. No es lineal: es plana hasta ~2000 promotores, sube entre 3000 y 5000, y cae en el tramo más alto.
2. **Se reproduce el mecanismo.**
   - Media de los unidos plana: ρ 0.00 / −0.03.
   - Media de los no unidos baja con la cobertura: ρ −0.39 / −0.34. Va de −0.003 a −0.033 en Rep 1, y de −0.024 a −0.045 en Rep 2.
   - Los no unidos de los factores de alta cobertura tienen menor avg_rank (ρ −0.80 / −0.83) y discordancia negativa.
3. **El control por avg_rank no atenúa nada.**
   - ρ crudo 0.28 / 0.22, ajustado 0.28 / 0.24.
   - Sobre el residuo de dif_signed ~ avg_rank, los no unidos siguen moviéndose (ρ −0.40 / −0.37).
4. **Marcadores** (en el conjunto, puesto 1 = mayor sobreestimación sobre 949):

   | Factor | Cobertura | Percentil | Puesto β | Puesto ajustado | FDR ajustado conjunto | FDR Rep 1 / Rep 2 |
   |---|---|---|---|---|---|---|
   | SETDB1 | 3520 (54%) | 0.50 | 43 | 29 | 1.3e-10 | 1e-9 / 9e-9 |
   | L3MBTL2 | 5505 (84%) | 0.79 | 97 | 202 | 0.033 | 0.25 / 0.0096 |
   | NCOR1 | 5020 (77%) | 0.69 | 194 | 358 | 0.37 | 0.53 / 0.35 |

5. **Conjuntos significativos.**
   - Sobreestima: 367 → 175 en el conjunto (317 → 142 / 358 → 190). Todos son un subconjunto de los originales.
   - Subestima: 97 → 177. Se conservan los 97 y aparecen 80 nuevos.
   - Ningún factor invierte el signo.
   - ORA, lado sobreestima: de 7 términos, todos genéricos, se mantienen 2 ("transcription regulator complex", "nuclear protein-containing complex").
   - ORA, lado subestima: de 54 términos se mantienen 32, entre ellos npBAF, nBAF, SWI/SNF, chromatin remodeling y nucleosome. Se pierden 22 y aparecen 5, entre ellos SWI/SNF superfamily-type complex.
   - De los 46 sobreestima concordantes del criterio HEK293, sobreviven 20. Sobreviven los 14 subestima.
6. **Filtro solo-HEK293.**
   - Los 263 factores que pasan tienen menor cobertura total (mediana 2814 frente a 3850 de los que no pasan; Wilcoxon p = 0.0014).
   - Proporción que pasa por quintil de cobertura: 0.37 / 0.27 / 0.25 / 0.27 / 0.23.
   - El sesgo existe pero es moderado, y va hacia la baja cobertura: no selecciona a favor de los factores que el artefacto infla.
   - Dentro de los 263, β en HEK casi no correlaciona con la cobertura en HEK (ρ 0.12).

Archivos: `p1_coef_vs_cobertura`, `p2_media_unidos_no_unidos`, `p2_avg_rank_unidos_no_unidos`, `p4_puesto_original_vs_ajustado`, `p6_cobertura_filtro_hek` (jpg y pdf), y los `.tsv` por punto.
