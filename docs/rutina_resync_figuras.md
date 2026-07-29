# Rutina: resincronizar numeración de figuras

Corre automáticamente todos los días hábiles a las 8am como tarea
programada (`resync-figuras-mapping`). Este documento describe la misma
rutina para cuando haga falta pegarla a mano en una sesión de Claude
Code con este repo.

Los decks fuente viven en la carpeta de Drive
https://drive.google.com/drive/folders/1yKXenI9vXq5CzfckUxYKU97EUSaK7Zus
(de martin.iungman@gmail.com), como Google Slides nativas cuyo título
empieza con "Fig" (`Fig MyM`, `Fig R1`, `Fig R2`, ... — pueden agregarse
más decks con el tiempo). Cada una corresponde 1:1 a un archivo local
`docs/Fig <sufijo>.pptx`.

---

> Necesito resincronizar la numeración de figuras en
> `docs/mapping_figuras.csv` a partir de los decks de Drive.
>
> 1. Traé los decks desde Drive: buscá en la carpeta
>    `1yKXenI9vXq5CzfckUxYKU97EUSaK7Zus` los archivos cuyo título empiece
>    con "Fig" (son Google Slides nativas, no pptx subidos - hay que
>    descargarlas con `exportMimeType` de PowerPoint) y sobrescribí
>    `docs/Fig <título>.pptx` en el repo con el contenido bajado (si el
>    título es nuevo, creá el archivo).
> 2. Corré `git status --porcelain -- "docs/Fig *.pptx"` para ver cuáles
>    decks cambiaron de verdad (modificados o nuevos). Si no cambió
>    ninguno, no hay nada que hacer.
> 3. Para cada deck que cambió: releelo completo, extrayendo el texto de
>    cada slide (`markitdown`) y las imágenes embebidas (`python-pptx`,
>    iterando shapes de tipo PICTURE). El número de slide = número de
>    figura actual según el deck: en `Fig MyM` slide N = Fig. M{N}; en
>    `Fig R{deck}` slide N = Fig. R{deck}.{N} (mismo esquema que ya usa
>    el CSV, ej. "Fig. R1.2"). Si aparece un deck con un nombre que no
>    encaja en ninguno de los dos esquemas, no inventes uno - señalalo.
> 4. Para cada slide, identificá qué slug de `docs/mapping_figuras.csv`
>    corresponde por contenido (paneles, texto de la leyenda, valores
>    numéricos si los hay), considerando solo las filas de la sección de
>    ese deck - no asumas que el slug que hoy tiene ese número sigue
>    siendo el correcto. Si un slide no tiene contenido (placeholder) o
>    no coincide con ningún slug existente, señalalo en vez de forzar
>    una correspondencia.
> 5. Actualizá `numero_actual` en el CSV para cada slug cuyo número
>    cambió. Si un slug ya no aparece en el deck, o el deck es
>    enteramente nuevo (sin filas existentes de esa sección), no
>    agregues/elimines filas ni inventes slugs - decime antes.
> 6. Corré `Rscript R/maintenance/sync_figure_folders.R` para renombrar
>    las carpetas de `figures/` según el CSV actualizado. Confirmá el
>    resultado (qué se renombró, qué quedó igual).
> 7. NO toques `R/<numero>_<slug>/` - esos números son de orden de
>    pipeline, no de figura, y no cambian.
> 8. Agregá al commit los `docs/Fig *.pptx` que cambiaron, el CSV y las
>    carpetas renombradas de `figures/`. No hagas push - dejá el commit
>    local en main para revisar.
> 9. Resumime qué decks se bajaron y cambiaron, qué figuras cambiaron de
>    número (por deck) y cuáles quedaron ambiguas o sin resolver (si las
>    hay).

---

**Por qué hace falta un prompt y no solo un script**: identificar "qué
contenido es ahora la Fig. M9" es una comparación semántica (leyendas,
paneles, valores) que ya nos hizo cambiar de opinión una vez en esta
tesis (ver commit donde `composicion_secuencia_library` pasó de M9 a
M10) — no es mecánico. `sync_figure_folders.R` sí es mecánico y seguro
de automatizar del todo: solo renombra carpetas para que coincidan con
lo que ya está en el CSV. Traer los decks desde Drive tampoco requiere
juicio: es una descarga y comparación de bytes.
