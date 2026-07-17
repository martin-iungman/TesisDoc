# TesisDoc

Repositorio de análisis (R) para la tesis doctoral **"Genome-wide
quantitative profiling of human core promoters..."**
(Iungman, Uriostegui-Arcos, Fiszbein, Schor). Acá vive el código, los
datos procesados y las figuras que referencia el texto de la tesis
(el texto en sí está en un Google Doc separado, ver `docs/tesis.md`
para una copia sincronizada en Markdown).

## Cómo está organizado

```
R/<numero>_<slug>/     scripts que generan cada figura de tipo "código",
                        numerados por orden de pipeline (no por número de figura)
R/functions/            helpers reutilizables (tema de figuras, fig_dir()...)
R/maintenance/          scripts de mantenimiento del repo (no generan figuras)
data/raw/                datos originales — nunca se editan a mano
data/external/           referencias externas (CAGE, EPD, ChIP-seq, PUFFIN...)
data/processed/          intermedios generados por scripts
figures/<numero>_<slug>/ salida final de cada figura (código, foto o esquema)
tables/                  tablas suplementarias y resultados de stats
docs/mapping_figuras.csv tabla puente: slug ↔ número de figura vigente ↔ script ↔ estado
docs/tesis.md            copia en Markdown del Google Doc de la tesis
docs/paper_referencia.md copia de trabajo del paper asociado
```

`data/` no se versiona (es pesado y regenerable) salvo los `.gitkeep`
que marcan la estructura de carpetas.

**Regla clave**: los números de figura (Fig. M9, Fig. R4...) cambian a
medida que se reorganiza la tesis. El código nunca usa el número, solo
el **slug estable** de cada figura (`gating_facs`,
`composicion_secuencia_library`, etc.). `docs/mapping_figuras.csv` es la
única fuente de verdad de qué número le corresponde hoy a cada slug —
ver `CLAUDE.md` para el detalle completo de esta y otras convenciones
del repo (pensado para trabajar con Claude Code, pero aplica igual
trabajando a mano).

## Requisitos

- R 4.4+ con: `tidyverse`, `BSgenome.Hsapiens.UCSC.hg38`, `rtracklayer`,
  `plyranges`, `GenomicRanges`, `Biostrings`, `biomaRt`,
  `EnsDb.Hsapiens.v86`, `patchwork`, `ggpubr`.
- Este repo asume que `../transcriptional_library` (el repo de datos del
  labo) existe como carpeta hermana — varios scripts leen archivos
  externos pesados directamente de ahí en vez de duplicarlos acá (ver
  comentarios "heavy_data_paths" en `R/00_prom_features/`).

## Correr un script de figura

Todos los scripts asumen que se corren desde la raíz del repo:

```r
# 1. (una sola vez, o cuando cambien los datos fuente) generar la tabla
#    de features compartida:
source("R/00_prom_features/build_prom_features.R")

# 2. correr el script de la figura que quieras regenerar, ej:
source("R/06_composicion_library/composicion_library.R")
```

La salida queda en `figures/<numero>_<slug>/`, resuelto automáticamente
según `docs/mapping_figuras.csv` — no hace falta (ni conviene) hardcodear
la carpeta numerada en el script.

## Mantenimiento

- `R/maintenance/sync_figure_folders.R`: renombra las carpetas de
  `figures/` para que coincidan con la numeración vigente en
  `docs/mapping_figuras.csv`. Correr después de editar `numero_actual`
  a mano.
- `docs/rutina_sync_tesis.md` y `docs/rutina_resync_figuras.md`:
  rutinas (documentadas como prompts) para resincronizar `tesis.md` y la
  numeración de figuras cuando cambian sus fuentes en Google
  Docs/PowerPoint. Ambas corren programadas de lunes a viernes a las 8am.
