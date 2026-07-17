# Rutina: sincronizar docs/tesis.md desde Google Docs

Prompt para correr en una sesión de Claude Code con este repo (manual o
programada). Actualiza `docs/tesis.md` con el contenido actual del
Google Doc "Tesis 1.0".

---

> Actualizá `docs/tesis.md` en este repo con el contenido actual del
> Google Doc "Tesis 1.0"
> (id `1EQ_qgwtkTYRR5xeG2lch9rZrSURJIqG88PKIwyQRni0`).
>
> El doc tiene varias pestañas internas (a la fecha: To do, Primer
> punteo, Organizacion, Metodos Exp, Metodos bioinfo, Resultados 1,
> Resultados 2) - la exportación tiene que traerlas todas, no solo la
> pestaña activa.
>
> 1. Si hay un conector de Google Drive/Docs disponible y desatendido
>    (no depende de una sesión de Chrome abierta), usalo para exportar
>    el doc completo a Markdown. Si no, usá el navegador (Chrome
>    extension): navegá a
>    `https://docs.google.com/document/d/1EQ_qgwtkTYRR5xeG2lch9rZrSURJIqG88PKIwyQRni0/export?format=md`
>    logueado con la cuenta que tiene acceso - esto dispara una descarga
>    directa a la carpeta de Descargas del sistema con TODAS las
>    pestañas ya en Markdown (validado: 7/7 pestañas, sin necesidad de
>    activar modo lector de pantalla ni leer el canvas a mano).
> 2. Mové el archivo descargado a `docs/tesis.md`, reemplazando el
>    anterior, y borrá cualquier copia residual que haya quedado en
>    Descargas por la descarga de este run.
> 3. Compará contra la versión anterior de `docs/tesis.md` (podés
>    mirarla con `git diff` después de sobrescribir) ignorando
>    diferencias triviales de espacios en blanco al final de línea.
> 4. Si el contenido real no cambió, hacé `git checkout` para descartar
>    la sobrescritura y terminá sin commit ni push.
> 5. Si cambió, commiteá solo `docs/tesis.md` con mensaje
>    `docs: sync tesis.md desde Google Docs (YYYY-MM-DD)` (fecha de
>    hoy) y pusheá directo a `main` - sin pedir confirmación, ya
>    autorizado de antemano para esta rutina puntual.
> 6. No modifiques ningún otro archivo del repositorio.
> 7. Contame en una o dos líneas si hubo cambios y, si los hubo, qué
>    secciones/pestañas se tocaron (a alto nivel, no hace falta el
>    diff completo).

---

**Cambios respecto al primer borrador de esta rutina** (después de
probarla una vez a mano): no asumir que existe un conector de Google
Drive desatendido — el método que sí funcionó fue el navegador con la
URL de exportación nativa (`/export?format=md`), que además es el único
que confirmamos que trae las 7 pestañas del doc. Push directo a main
sin revisión: aceptado explícitamente por el autor, porque la versión
que importa de verdad es la de Google Docs, no la del repo.
