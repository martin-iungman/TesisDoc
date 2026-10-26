# Rutina: sincronizar docs/tesis.md desde Google Docs

Corre automáticamente todos los días hábiles a las 8am como tarea
programada (`sync-tesis-md`). Este documento describe la misma rutina
para cuando haga falta pegarla a mano en una sesión de Claude Code con
este repo. Hace dos cosas: actualiza `docs/tesis.md` con el contenido
actual del Google Doc "Tesis 2.0" y, a partir de ese texto, actualiza la
columna `leyenda_lista` de `docs/mapping_figuras.csv`.

---

> **Parte A — sincronizar `docs/tesis.md`**
>
> Actualizá `docs/tesis.md` con el contenido actual del Google Doc
> "Tesis 2.0" (id `1EulYfx64i5t3lmHmP3w03kS3eN9pRH_rDUaYALJxhak`). El
> doc tiene varias pestañas internas - la exportación tiene que traerlas
> todas, no solo la pestaña activa.
>
> 1. Preferido: conector de Google Drive, `download_file_content` con
>    ese `fileId` y `exportMimeType: "text/markdown"`. El resultado es
>    grande y queda guardado en un JSON (`{content (base64), id,
>    mimeType, title}`): decodificá `content` con Python y escribí los
>    bytes en `docs/tesis.md`. Alternativa sin conector: Chrome
>    extension, navegar a
>    `https://docs.google.com/document/d/1EulYfx64i5t3lmHmP3w03kS3eN9pRH_rDUaYALJxhak/export?format=md`
>    logueado con la cuenta que tiene acceso (descarga directa con todas
>    las pestañas); mover el archivo a `docs/tesis.md` y borrar la copia
>    que haya quedado en Descargas.
> 2. Compará con `git diff -- docs/tesis.md`, ignorando diferencias
>    triviales de espacios al final de línea.
> 3. Si el contenido real no cambió, `git checkout -- docs/tesis.md` y
>    seguí igual con la Parte B (el estado de las leyendas puede cambiar
>    porque cambió `numero_actual` en el CSV).
> 4. Si cambió, commiteá solo `docs/tesis.md` con mensaje
>    `docs: sync tesis.md desde Google Docs (YYYY-MM-DD)`.
>
> **Parte B — actualizar `leyenda_lista`**
>
> 5. Si `docs/mapping_figuras.csv` tiene cambios sin commitear antes de
>    empezar, no lo toques: saltá la Parte B y avisalo.
> 6. Buscá las leyendas en `docs/tesis.md`: párrafos que empiezan con el
>    título en negrita, `**Figura R1.2. <título>.** <descripción>`.
>    Asocialas a filas del CSV por número ("Figura R1.2" ↔
>    `numero_actual` "Fig. R1.2"). Las filas con número no resuelto
>    ("?", "/", vacío, o `notas` con "revisar"/"aclarar" sobre el
>    número) quedan como están y se reportan.
>    - A veces el doc conserva debajo un borrador viejo de la sección con
>      numeración vieja (dos "Figura R1.3" distintas). Si un número tiene
>      más de una leyenda, usá la que corresponde por contenido al slug
>      (comparando con slug, `notas`, `archivos`); en general es la
>      primera en orden de lectura. Nunca marques `si` con una leyenda
>      que describe otra figura. Reportá los duplicados.
>    - Una mención ("ver Fig. R2.6") o una línea placeholder del estilo
>      `Fig. R1.8 - comparacion histo-density...` no es una leyenda.
> 7. Clasificá cada fila resuelta:
>    - `si`: leyenda completa - título en negrita más descripción, todos
>      los paneles de la figura descriptos y sin marcadores de pendiente
>      ("[falta la leyenda]", "[completar]", "TODO", "XX", "??").
>    - `desactualizada`: hay leyenda para ese número pero está
>      incompleta (placeholder, paneles sin describir, título sin
>      descripción).
>    - `no`: no hay leyenda de esa figura en el texto.
>    Cambiá solo la columna `leyenda_lista` (con el módulo `csv` de
>    Python, respetando orden de columnas, comillas y fin de línea) y
>    verificá con `git diff` que no se tocó nada más.
> 8. Si el CSV cambió, commiteá solo `docs/mapping_figuras.csv` con
>    mensaje `docs: leyenda_lista desde tesis.md (YYYY-MM-DD)`.
>
> **Push y reporte**
>
> 9. Si hubo algún commit, pusheá directo a `main` - sin pedir
>    confirmación, ya autorizado de antemano para esta rutina puntual.
>    No modifiques ningún otro archivo del repositorio.
> 10. Contame brevemente: si cambió `tesis.md` (qué pestañas/secciones,
>    a alto nivel), qué figuras cambiaron de `leyenda_lista` (slug:
>    viejo → nuevo), qué filas quedaron sin resolver y qué leyendas
>    duplicadas aparecieron.
>
> Si la exportación falla, reportalo - no inventes contenido para
> `docs/tesis.md` ni corras la Parte B sobre un export viejo o parcial.

---

**Notas de diseño**: push directo a main sin revisión, aceptado
explícitamente por el autor, porque la versión que importa de verdad es
la de Google Docs, no la del repo. `leyenda_lista` se decide por el
texto de la tesis (no por el slide del deck): la rutina
`resync-figuras-mapping` no toca esa columna.
