# Rutina: resincronizar numeración de figuras

Cuando `docs/Fig MyM.pptx` (o el equivalente de Resultados) cambie —
slides agregados, movidos o reordenados — pegar el siguiente prompt en
una sesión de Claude Code con este repo:

---

> `docs/Fig MyM.pptx` cambió. Necesito resincronizar la numeración de
> figuras en `docs/mapping_figuras.csv`.
>
> 1. Releé el pptx completo: extraé el texto de cada slide (`markitdown`)
>    y las imágenes embebidas (`python-pptx`, iterando shapes de tipo
>    PICTURE). El número de slide = número de figura actual (slide N =
>    Fig. M{N}; para el pptx de Resultados sería R{N} si corresponde).
> 2. Para cada slide, identificá qué slug de `docs/mapping_figuras.csv`
>    corresponde por contenido (paneles, texto de la leyenda, valores
>    numéricos si los hay) - no asumas que el slug que hoy tiene ese
>    número sigue siendo el correcto. Si un slide no tiene contenido
>    (placeholder) o no coincide con ningún slug existente, señalalo en
>    vez de forzar una correspondencia.
> 3. Actualizá `numero_actual` en el CSV para cada slug cuyo número
>    cambió. Si un slug ya no aparece en el pptx o aparece uno nuevo sin
>    slug, decime antes de agregar/eliminar filas.
> 4. Corré `Rscript R/maintenance/sync_figure_folders.R` para renombrar
>    las carpetas de `figures/` según el CSV actualizado. Confirmá el
>    resultado (qué se renombró, qué quedó igual).
> 5. NO toques `R/<numero>_<slug>/` - esos números son de orden de
>    pipeline, no de figura, y no cambian.
> 6. Resumime qué figuras cambiaron de número y cuáles quedaron
>    ambiguas (si las hay).

---

**Por qué hace falta un prompt y no solo un script**: identificar "qué
contenido es ahora la Fig. M9" es una comparación semántica (leyendas,
paneles, valores) que ya nos hizo cambiar de opinión una vez en esta
tesis (ver commit donde `composicion_secuencia_library` pasó de M9 a
M10) — no es mecánico. `sync_figure_folders.R` sí es mecánico y seguro
de automatizar del todo: solo renombra carpetas para que coincidan con
lo que ya está en el CSV.
