# TesisDoc — contexto del proyecto

Repositorio de análisis (R) para la tesis doctoral "Genome-wide quantitative
profiling of human core promoters..." (Iungman, Uriostegui-Arcos, Fiszbein,
Schor). Gran parte del trabajo está basado en el paper asociado (ver
`docs/paper_referencia.md`). La tesis en español está en un Google
Doc/Word separado (no vive en este repo); acá vive el código, los datos y
las figuras que ese texto referencia.

## Regla central: numeración de figuras vs. identidad del análisis

**El número de figura (Fig. M2, Fig. R4, Fig. 2C del paper...) puede
cambiar en cualquier momento** a medida que se reorganiza la tesis o el
paper. Por eso:

- En **código** (nombres de función, variables, argumentos, contenido de
  scripts) nunca uses el número de figura. Usá siempre el **slug
  estable** de esa figura (ver `docs/mapping_figuras.csv`), por ejemplo
  `gating_facs`, `sesgo_representatividad`, `tf_nfya_sp_contexto_endogeno`.
  Ningún script hardcodea `"figures/M10_..."` — para escribir su salida
  usa `fig_dir("<slug>")` (`R/functions/fig_paths.R`), que resuelve la
  carpeta numerada actual leyendo el CSV. Así un script sigue funcionando
  sin cambios aunque la numeración se reordene.
- La carpeta `figures/<numero>_<slug>/` **sí** lleva el número vigente
  como prefijo (ver "Estructura de carpetas" abajo) para poder navegarla
  ordenada — es la única excepción a la regla de arriba, y es un dato
  derivado/sincronizado, no algo que se edite a mano: al cambiar
  `numero_actual` en el CSV, correr `Rscript
  R/maintenance/sync_figure_folders.R` para renombrar las carpetas
  afectadas. `R/<numero>_<slug>/` (los scripts) NO se renombran nunca —
  ese número es de orden de pipeline, no de figura (ver más abajo).
- La numeración actual, la sección de la tesis/paper a la que corresponde,
  y el estado de cada figura viven **únicamente** en
  `docs/mapping_figuras.csv`. Es la fuente de verdad para "qué número es
  hoy esta figura". Si el usuario menciona "Fig. R6" o similar, buscar
  primero en ese CSV a qué slug corresponde antes de asumir cuál es.
- Si el CSV no tiene una correspondencia clara para un número mencionado
  (columna `notas` con "revisar" o "aclarar", o `numero_actual` con "?"
  o "/"), señalarlo en vez de adivinar. `fig_dir()` trata esos casos como
  no resueltos y guarda en `figures/sinnum_<slug>/` hasta que se aclare.
- Si hace falta un export con los PDFs nombrados por número (por ejemplo
  para mandar a la revista), generarlo con un script que lea el CSV y
  copie a `figures/_export_numerado/` — ese export es un producto
  derivado y descartable, nunca la fuente de verdad.
- Cuando cambie `docs/Fig MyM.pptx` (o el correspondiente de Resultados),
  releer su contenido slide por slide (el número de slide = número de
  figura) y comparar contra `docs/mapping_figuras.csv` antes de asumir
  que la numeración sigue igual — ver `docs/rutina_resync_figuras.md`.

## Tipos de figura — no todo se genera por código

La columna `tipo` en `docs/mapping_figuras.csv` distingue tres orígenes:

- **codigo**: se genera desde datos con un script R. Vive en `R/<carpeta
  del slug>/` y su salida en `figures/<numero>_<slug>/` (resuelto vía
  `fig_dir("<slug>")`, nunca hardcodeado).
- **foto**: fotografía real (microscopía, gel de electroforesis, etc.).
  No hay script que la genere. El archivo de imagen final va directo en
  `figures/<numero>_<slug>/`, sin carpeta correspondiente en `R/`.
- **esquema**: diagrama hecho a mano (PowerPoint, Illustrator, BioRender),
  no proviene de datos ni de código. La imagen final va directo en
  `figures/<numero>_<slug>/`, igual que una foto. El archivo editable (si
  existe) no se versiona en este repo.
- **mixto**: figura con paneles de distinto origen (ej. `vector_construccion`
  tiene un esquema + dos fotos). Anotarlo panel por panel en la columna
  `notas`.

No asumir que porque una figura aparece en el CSV existe un script para
regenerarla — revisar la columna `tipo` antes de buscar o crear código
para una figura de tipo `foto` o `esquema`.

## Estructura de carpetas

```
R/<numero>_<slug>/     scripts que generan cada figura de tipo "codigo",
                       numerados solo por orden de pipeline (no por
                       número de figura) - este número NUNCA cambia
R/functions/           helpers reutilizables (temas ggplot, fig_dir()/
                       fig_number_prefix() en fig_paths.R, funciones de
                       filtrado/normalización comunes)
R/maintenance/         scripts de mantenimiento del repo (no generan
                       figuras), ej. sync_figure_folders.R
data/raw/              datos originales — nunca se editan a mano
data/external/         referencias externas (CAGE, EPD, ChIP-seq, PUFFIN...)
data/processed/        intermedios generados por scripts
figures/<numero>_<slug>/  salida final de cada figura, código, foto o
                       esquema. <numero> viene de mapping_figuras.csv y
                       se sincroniza con R/maintenance/sync_figure_folders.R
                       - nunca a mano. Si el número no está resuelto en
                       el CSV, la carpeta usa el prefijo "sinnum".
tables/                tablas suplementarias y resultados de stats
docs/mapping_figuras.csv            tabla puente figura ↔ script ↔ sección
docs/paper_referencia.md            copia de trabajo del paper
docs/notas_pendientes.md            lista de pendientes migrada de la tesis
docs/glosario_convenciones.md       nombres de variables/criterios comunes
docs/rutina_resync_figuras.md       prompt para resincronizar numeración
                       cuando cambia docs/Fig MyM.pptx (u otro pptx de figuras)
```

## Convenciones generales

- Cada script en `R/<slug>/` debe poder correr de forma independiente
  siempre que sus insumos ya estén en `data/processed/`.
- Guardar resultados intermedios costosos (ej. procesamiento de
  secuenciación) en `data/processed/` con nombre versionado, para no
  tener que recorrer todo el pipeline al retocar una figura.
- Al agregar o modificar una figura, actualizar la fila correspondiente
  en `docs/mapping_figuras.csv` (estado, script, notas) en el mismo
  cambio.
