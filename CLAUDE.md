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

- Nunca uses el número de figura como nombre de archivo, carpeta, función
  o variable. Usá siempre el **slug estable** de esa figura (ver
  `docs/mapping_figuras.csv`), por ejemplo `gating_facs`,
  `sesgo_representatividad`, `tf_nfya_sp_contexto_endogeno`.
- La numeración actual, la sección de la tesis/paper a la que corresponde,
  y el estado de cada figura viven **únicamente** en
  `docs/mapping_figuras.csv`. Es la fuente de verdad para "qué número es
  hoy esta figura". Si el usuario menciona "Fig. R6" o similar, buscar
  primero en ese CSV a qué slug corresponde antes de asumir cuál es.
- Si el CSV no tiene una correspondencia clara para un número mencionado
  (columna `notas` con "revisar" o "aclarar"), señalarlo en vez de
  adivinar.
- Si hace falta un export con los PDFs nombrados por número (por ejemplo
  para mandar a la revista), generarlo con un script que lea el CSV y
  copie a `figures/_export_numerado/` — ese export es un producto
  derivado y descartable, nunca la fuente de verdad.

## Tipos de figura — no todo se genera por código

La columna `tipo` en `docs/mapping_figuras.csv` distingue tres orígenes:

- **codigo**: se genera desde datos con un script R. Vive en `R/<carpeta
  del slug>/` y su salida en `figures/<slug>/`.
- **foto**: fotografía real (microscopía, gel de electroforesis, etc.).
  No hay script que la genere. El archivo de imagen final va directo en
  `figures/<slug>/`, sin carpeta correspondiente en `R/`.
- **esquema**: diagrama hecho a mano (PowerPoint, Illustrator, BioRender),
  no proviene de datos ni de código. Vive en
  `assets/diagramas_manuales/<slug>/` con el archivo editable si existe.
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
                       número de figura)
R/functions/           helpers reutilizables (temas ggplot, funciones
                       de filtrado/normalización comunes)
data/raw/              datos originales — nunca se editan a mano
data/external/         referencias externas (CAGE, EPD, ChIP-seq, PUFFIN...)
data/processed/        intermedios generados por scripts
figures/<slug>/        salida final de cada figura, código, foto o esquema
assets/diagramas_manuales/<slug>/   esquemas hechos a mano
tables/                tablas suplementarias y resultados de stats
docs/mapping_figuras.csv            tabla puente figura ↔ script ↔ sección
docs/paper_referencia.md            copia de trabajo del paper
docs/notas_pendientes.md            lista de pendientes migrada de la tesis
docs/glosario_convenciones.md       nombres de variables/criterios comunes
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
