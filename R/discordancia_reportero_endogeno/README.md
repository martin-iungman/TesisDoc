# discordancia_reportero_endogeno (carpeta de staging, sin numero)

Sacada del pipeline numerado el 2026-07-30 (era R4.6/R4.7 -
R/27_rank_dif_lmm y R/28_remap_rank_dif_lmm_gsea) mientras se decide
que version del analisis usar para la tesis. Los scripts de esta
carpeta no escriben via `fig_dir()`/`docs/mapping_figuras.csv` - cada
uno define su propio `out_dir` dentro de
`figures/discordancia_reportero_endogeno/<misma_subcarpeta>/`.

Las tres subcarpetas analizan lo mismo (¿que features/TFs/marcas de
histona predicen que el reportero se comporte distinto de la actividad
endogena real?) con tres disenos progresivamente mas cuidadosos:

## 1_original_sin_corregir

Portado de `transcriptional_library/Analysis/scripts/rank_dif_lmm.R`.
Controla por `rank_reporter_scaled` unicamente. **Esta confundido**:
`remap_tf.R` incluye el diagnostico (`remap_rank_dif_vs_activity.jpg`)
que lo prueba - el efecto de cada TF sobre esta "diferencia de rango"
correlaciona r~0.8 con su efecto sobre actividad cruda (R2.6/R7), es
decir, el analisis practicamente no aporta nada mas alla de re-detectar
que features/TFs predicen actividad.

## 2_dif_signed_bland_altman

Diseno corregido (Bland-Altman): en vez de controlar por
`rank_reporter_scaled`, descompone `(rank_reporter, rank_endo)` en
`avg_rank` (nivel, promedio de ambos rangos normalizados) y
`dif_signed` (discordancia con signo, `re - rr`). Controlando por
`avg_rank`, la correlacion residual con el efecto de actividad cruda
cae de ~0.8 a ~0.25-0.3. Responde: "¿este feature/TF hace que el
reportero sub o sobreestime la actividad endogena, en promedio?"
(direccion del sesgo, no magnitud).

## 3_doble_glm_dispersion

Mismo `avg_rank`/`dif_signed`, pero con un modelo doble-GLM (`glmmTMB`
con `dispformula`) que ademas de la media ajusta la **dispersion**
(varianza) de `dif_signed` - responde "¿este feature/TF genera mas
incertidumbre/inconsistencia entre reportero y endogeno, en cualquier
sentido (no solo un sesgo direccional)?". Se probo porque la MAGNITUD
de la discordancia (|dif_signed|) resistio 5 intentos previos de
desconfundirla de la actividad cruda (lineal, spline, residuo, matching
exacto, avg_rank+profundidad de lectura) - este modelo de dispersion,
mas riguroso, da el mismo resultado (r~-0.8) para el efecto sobre
DISPERSION: no se pudo desconfundir tampoco, parece ser una propiedad
real de los datos, no un artefacto de metodo. Donde SI aporta algo es
en encontrar excepciones a esa tendencia general (residuo de
dispersion~actividad) - TATA-box es la mayor excepcion entre las
features curadas, consistente con R3.1 (`ruido_cgi_tata`, diseno
totalmente distinto).

## Pendiente

Cuando se decida cual usar, promover esa version de vuelta al pipeline
numerado (nuevo R/NN_.../, figures/R4.N_.../, fila en
`docs/mapping_figuras.csv`) y borrar esta carpeta.
