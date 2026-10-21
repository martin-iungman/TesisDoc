# G4 por hebra a lo largo de la medición (R/37_g4_representatividad)

Todos los promotores juntos, sin estratificar por isla CpG. El G+C se
controla con bins de ancho fijo (0,05):
- Contra Sin G4 se usa una referencia apareada por bin (3:1).
- Codificante vs. molde se compara directo, porque tienen el mismo G+C.
- La retención se testea con CMH estratificado por bin.

El texto para Resultados está en `texto_resultados.md`. La figura principal es
`figures/R1.4_g4_representatividad/g4_retencion_por_gc.jpg` (y .pdf) — Fig. R1.4.

## Resumen (G4-seq K+; el patrón canónico da lo mismo)

**Presort 2022.** Comparado con Sin G4 de igual G+C:
- Solo codificante: 72,6% con lecturas contra 90,3% esperado (OR CMH 0,23) y
  −0,76 log10 lecturas (IC 95% −0,81 a −0,70).
- Solo molde: −0,16 log10 lecturas.
- Codificante vs. molde directo: −0,56.

**Sorting.** Lecturas de la codificante: −0,70 / −0,68 log10. La molde no
pierde (+0,04 / +0,02). Retención de la codificante contra Sin G4 por paso:
- Library → detectado: −21 / −21 pp.
- Detectado → sum_counts ≥ 100: −10 / −6 pp.
- Filtros de consistencia: −1 a −2 pp.
- Library → alta confianza: −19 / −20 pp (OR 0,41 / 0,42). La molde, en
  cambio, +5 / +4 pp (OR 1,25 / 1,20).

**Actividad** (HL contra Sin G4 apareado, alta confianza):

| Clase | Rep 1 | Rep 2 | Robustez |
|---|---|---|---|
| Solo molde | −0,078 | −0,089 | Igual sin filtros de consistencia y con presort ≥ 100 lecturas |
| Solo codificante | −0,017 (n.s.) | 0,000 (n.s.) | Igual en los mismos escenarios |
| Codificante vs. molde | +0,066 | +0,098 | — |

Con el patrón canónico, la codificante da −0,046 en Rep 1 (p ajustado 0,02)
y −0,017 en Rep 2 (n.s.).

**Otros datos.**
- Spearman entre lecturas y actividad: ρ = 0,16 / 0,09.
- Solo el peor caso extremo (todos los promotores perdidos de más con
  actividad mínima) vuelve negativo el efecto de la codificante.
- Prevalencia estandarizada por G+C: codificante 15,4% → 10,3% / 11,0%;
  molde 12,1% → 13,6% / 13,5%.
