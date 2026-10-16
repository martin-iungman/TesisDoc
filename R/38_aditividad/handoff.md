# Handoff — análisis exploratorio de aditividad (2026-10-15)

## Estado
- Script: `R/38_aditividad/aditividad.R` (corre en ~1 min desde la raíz del repo).
- **Figura de tesis R2.12** (slug `aditividad_elementos_promotor`, desde
  2026-10-16; slide 12 de `docs/Fig R2.pptx`): `figures/R2.12_aditividad_elementos_promotor/`
  vía `fig_dir()`, dos paneles sin letras en la imagen — A =
  `violines_n_elementos.jpg` (Parte 1), B = `grilla_pos_x_neg.jpg` (Parte 4).
- El resto de las salidas sigue siendo exploratorio, en `figures/exploratorio_aditividad/`
  (no versionado).
- Script, este handoff, la fila del CSV y la carpeta R2.12 commiteados.
- Depende de `YY1_motif`, `CRE_motif`, `Ebox_CAGGTR` y `polyAT_prox` en
  `data/processed/prom_df.tsv` (`R/00_prom_features/build_prom_features.R`,
  commit c92c496).

## Decisiones tomadas
- Carpeta `38_` (no `37_`: ya existía `R/37_g4_representatividad/`).
- **Positivos (5):** isla CpG (`CGI_frag`), `TATA_EPD`, `CCAAT_EPD`, `YY1_motif`, `CRE_motif`.
  - GC-box excluida: casi redundante con la isla CpG.
  - INR fuerte (`INR_strong_TSS`) evaluado y excluido: +0.001 de R², baja el
    Spearman del conteo, efecto solo dentro de isla CpG.
- **Negativos (3):** `other_TSS`, `Ebox_CAGGTR`, `polyAT_prox`.
  - TE (`TE_superclass` no NA) excluido: efecto crudo −0.09/−0.13, pero ≈ 0
    (n.s.) ajustado por isla CpG (37% de isla con TE vs. 74% sin TE).
- Categorías con < 100 promotores (en alguna réplica) se agrupan con la
  vecina; en la grilla se omiten celdas con < 50.
- Efectos: HL (IC 95%) + Wilcoxon, BH dentro de réplica; cada réplica por separado.

## Partes y resultados principales (Rep 1 / Rep 2)
1. **N.º de positivos (0, 1, 2, 3 o más):** gradiente sin saturación; HL vs. 0
   +0.17/+0.30/+0.53 y +0.19/+0.37/+0.67; ρ = 0.25 / 0.23.
2. **Motivos por grupo de isla:** fuerte con isla (HL "2 o más" +0.34 / +0.46),
   débil sin isla (+0.12 / +0.14; 1 motivo n.s. en Rep 2).
3. **Modelo lineal:** R² aditivo 0.084 / 0.073 → con interacciones de a pares
   0.088 / 0.078 (F p = 2e-4 / 1e-8). Sinergias reproducibles: isla × TATA,
   isla × YY1, TATA × CCAAT. Obs. vs. pred. por combinación r = 0.96 / 0.95.
   Con negativos como efectos principales: R² 0.089 / 0.077; TSS no canónico
   −0.06 / −0.10, E-box −0.10 / −0.10, poli(dA:dT) −0.19 / −0.19.
4. **Negativos y puntaje neto:** HL de negativos vs. 0 −0.10/−0.20 y
   −0.12/−0.26; ρ = −0.12. Grilla: los negativos bajan la mediana dentro de
   cada nivel de positivos, más cuanto más positivos. Puntaje neto
   (pos − neg) monótono de −2 o menos a 3 o más; ρ = 0.26 / 0.25 (vs. 0.25 /
   0.23 solo positivos).

## Robustez de la Parte 3: percentil como respuesta y techo entre réplicas
Mismos tres modelos con el percentil de actividad dentro de cada réplica
(rango promedio en empates) como respuesta; coeficientes en puntos de
percentil. Tablas: `p3_modelos_r2_F_media_vs_percentil.tsv`,
`p3_interacciones_media_vs_percentil.tsv`, `p3_interacciones_percentil.tsv`,
`p3_coef_aditivo_pos_neg_percentil.tsv`, `p3_techo_acuerdo_replicas.tsv`.

| Modelo (Rep 1 / Rep 2) | Actividad media | Percentil |
|---|---|---|
| R² aditivo | 0.084 / 0.073 | 0.091 / 0.084 |
| R² con interacciones | 0.088 / 0.078 | 0.092 / 0.086 |
| F interacciones (p) | 3.3 (2e-4) / 5.7 (1e-8) | 2.1 (0.021) / 2.9 (0.001) |
| R² positivos + negativos | 0.089 / 0.077 | 0.097 / 0.090 |
| F negativos (p) | 16.3 (1e-10) / 17.9 (1e-11) | 23.6 (3e-15) / 27.0 (2e-17) |

Interacciones candidatas (estimador [IC 95%], p BH):

| Interacción | Rep | Actividad media | Percentil |
|---|---|---|---|
| Isla × TATA | 1 | +0.18 [0.07, 0.29], 0.010 | +1.8 [−2.6, 6.1], 0.54 n.s. |
| Isla × TATA | 2 | +0.21 [0.08, 0.34], 0.006 | +3.1 [−0.9, 7.1], 0.31 n.s. |
| Isla × YY1 | 1 | +0.13 [0.06, 0.21], 0.008 | +5.4 [2.4, 8.4], 0.005 |
| Isla × YY1 | 2 | +0.14 [0.05, 0.23], 0.006 | +3.1 [0.3, 5.8], 0.092 n.s. |
| TATA × CCAAT | 1 | +0.16 [0.04, 0.28], 0.039 | +2.3 [−2.5, 7.1], 0.53 n.s. |
| TATA × CCAAT | 2 | +0.33 [0.18, 0.48], 1e-4 | +7.7 [3.2, 12.2], 0.008 |
| Isla × CRE | 2 | +0.13, 0.023 | +4.1, 0.034 |
| Isla × CCAAT | 2 | +0.12, 0.029 | +0.7, 0.82 n.s. |

- En percentil la mejora por interacciones se achica (ΔR² 0.002 vs.
  0.003–0.005) y la mayoría de las sinergias pierde significancia: **isla ×
  TATA no sobrevive en ninguna réplica**; isla × YY1 y TATA × CCAAT
  sobreviven solo en una réplica cada una (todas siguen con signo positivo).
  Parte de la "sinergia" en la escala de bins parece ser compresión cerca
  del piso: con los efectos principales solos, el aditivo ya rinde igual o
  mejor en percentil.
- Efectos principales y negativos se mantienen en percentil (isla +14,
  TATA +9 a +13, CCAAT +10, YY1 / CRE +2 a +3; TSS no canónico −4,
  E-box −4, poli(dA:dT) −6 a −9 puntos de percentil).
- **Techo** (9152 promotores medidos en ambas réplicas): Pearson Rep 1 vs.
  Rep 2 r = 0.87 en actividad media (r² = 0.76), r = 0.70 en percentil
  (r² = 0.49; Spearman 0.70). Con un predictor perfecto del valor verdadero
  el R² esperado contra una réplica es ≈ r (no r²), o sea ~0.87 / ~0.70: los
  modelos (R² ~0.08–0.10) explican ~10–13% de la varianza explicable.
- `p3_interacciones.jpg`: la leyenda muestra solo las categorías presentes
  (sin "Subaditiva") y se separaron los paneles para que no se pisen los ejes.
  Las otras figuras no se regeneraron.

## Pendientes posibles
- Caveats a discutir: escala de actividad acotada (bins 1–6); YY1/CRE/E-box
  buscados en todo el fragmento vs. EPD con posición; R² bajo (~9%).
