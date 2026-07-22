# To do

Por hacer:

- Pipeline gral a español  
- Potencia: todo  
- Esquema de density (creo que ya tengo y solo haria falta pasar a español)  
- M2C como figura aparte?  
- M5 (qpcr): agrandar textos  
- M6: editar para indicar calles  
- Anexo Primers y secuencias: falta index de secuenciacion usado y primer qPCR EGFP   
- M9C (heatmap cage): modificar por el del review. modificar leyenda, texto y discutir sobre el upstream  
- M9B: agrandar texto  
- Caracterización de los promotores: reorganizar en subsecciones  
- M10 agrandar texto  
- M11A agrandar texto  
- Nueva figura sobre shape y tissue specificity??  
- Seccion coocurrencia  
- Tabla de features con criterio de dicotomia  
- M12 pasar a español  
- MyM: PUFFIN. Pensar si poner figura, armar texto  
- Revisar texto ROC-AUC. Hacer alguna figura explicativa?

- Los datos de validacion de ruido y replicabiliad que dan feos… los incluyo???  
- R2C: cambiar el theme  
- R6: tienen distinto x label. hay problema con el x text. texto de violin muy chico. Reemplazar el TRUE \- FALSE  
- Leyenda R6  
- Incluir CCAAT en R6  
- Plot para NFY y SP? Pre y post review?   
- Leyenda R7 y R8

Para ir tachando Fig y textos que ya estan (draft)  
Ref:  
Cerrado 100%  
Listo figura pero falta leyenda  
Tengo los elementos de la figura pero falta trabajo (por ej, pasar a español)  
Revisar mas adelante (por ahora ok)  
Falta referenciar y/o explicar la figura  
Iniciado

**Materiales y metodos experimentales:**  
Descripción general

- Pipeline

Construcción de los plasmidos conteniendo la library de promotores

- Esquema vector  
- micro cmv-wt y cmv-strong   
- PCR\_GA\_all\_prom\_cloning\_control 

	Cultivo celular y transfección  
	Citometría de flujo y sorting

- Potencia??  
- Esquema de density con los gates

	Control por spike-in de celulas

	Co-extraccion de ADN/ARN y RT-qPCR de EGFP

- qPCR   
- IMAGEN PCR\_loxP\_gates\_ctl

Secuenciacion de amplicones  
	Anexo Primers y secuencias

**Metodología bioinformática:**  
	Procesamiento de los datos de secuenciacion

- Pipeline preprocessing  
- S1A (sample effort) \+ S6A (bimodal)

	Caracterización de la library

- General: Esquema promotores y enhancer 252pb con posicion del TSS y primers \+ acompañado por el heatmap de CAGE \+ barplot seq types   
- Motivos: cont GC \+ EPD barplot  \+ CGI+ INR  
- Evo: phyloP \+ young;  
- tissue specificity \+ shape??  
- coocurrence?  
- PUFFIN??


Asociacion de las características del promotor con actividad y ruido transcripcional  
	Clasificación de promotores alternativos

**Resultados**  
Estimación masiva de las propiedades transcripcionales de promotores basales humanos

- N prom por replica   
-  comparacion con splicing22 \+ library bias  
- spike-in   
-  Mean replicate   
- 1B-C-D (validacion)

Efectos de la secuencia promotor sobre la fuerza transcripcional

- figura explicativa clara sobre metodologia ??  
-  CpG \+ TATA \+CCAAT  
-  NFYA \+ SP1/2  
- summary seq  
- tissue sp \+ ?? \+ summary endo  
- integrador?

Efectos de la secuencia promotor sobre el ruido transcripcional

- 3A-B \+TATA2(S6B)  
- summary (3C)  
- S6C \+ 3D-E (GSEA MLL1)

La relevancia de la secuencia del promotor basal en el contexto endogeno

- 2E \+ S4A (hela \- hek etc)  
- S2(A-B, gsea pol)  
- 2F S4B (housekeeping y predictibilidad)  
- 2G \+ S5 (puffin)  
- sure?  
- xpresso?  
- podria sumar (en alguna) tissue sp y ruido

Influencia de la multiplicidad de promotores por gen sobre las caracteristica transcripcionales

- 4A-B  
- estratificacion por tissue specificity  
- caracterizacion (?)  
- S7A \+ promoter similarity (o esto va a otra seccion? veremos que tan largo queda)  
- clasificacion: esquemas fig4C \+ S7B  
- caracterizacion (?)  
- 4C \+ ruido?

# Primer punteo

Además de las figuras del paper, que podria sumar?

- análisis de potencia del experimento?  
- similitud de promotores  
- ~~downstream y upstream~~  
- caracterización de promotores alternativos  
- análisis integrador (que da feo)  
- ~~transient validation~~  
- pcr control loxp  
- alguna foto de microscopía  
- comparación con splicing 2022  
- library bias  
- overlapped promoters  
- **analizar un poco casos cancer y enhancer?**  
- sumar algunas caracteristicas extra de genes?  
- xpresso  
- sure?

# Organizacion

**Introducción**

**Materiales y metodos**  
Descripción general

- Pipeline

**Experimentales:**  
Construcción de los plasmidos conteniendo la library de promotores

- Esquema vector \+ PCR \+ microscopía

	Cultivo celular y transfección  
	Citometría de flujo y clasificación de células

- Potencia??  
- Esquema de density con los gates?

	Control por spike-in de celulas  
	Co-extraccion de ADN/ARN y secuenciacion de amplicones  
	RT-qPCR de EGFP

- qPCR   
- PCR loxp (cel no transfectadas)

	Anexo Primers y secuencias

**Metodología bioinformática:**  
	Procesamiento de los datos de secuenciacion

- Pipeline preprocessing  
- S1A (sample effort) \+ S6A (bimodal)

	Caracterización de la library

- Descripción library  
- caracteristicas de seq de la lib  
- evo library (young 2015 \+ phylop)

	Asociacion de las caracteristicas del promotor con actividad y ruido transcripcional  
	Clasificación de promotores alternativos

**Resultados**  
Estimación masiva de las propiedades transcripcionales de promotores basales humanos

- Barplot N prom por replica \+ comparacion con splicing22 y otras seq \+ library bias  
- Mean replicate \+ spike-in  
- 1B-C-D (validacion)

Efectos de la secuencia promotor sobre la fuerza transcripcional

- figura explicativa clara sobre metodologia \+ CpG \+ TATA \+ NFYA   
- \+ SP1/2 \+ NFYA   
-   
- summary seq  
- tissue sp \+ ?? \+summary endo  
- integrador?

Efectos de la secuencia promotor sobre el ruido transcripcional

- 3A-B \+TATA2(S6B)  
- summary (3C)  
- S6C \+ 3D-E (GSEA MLL1)

La relevancia de la secuencia del promotor basal en el contexto endogeno

- 2E \+ S4A (hela \- hek etc)  
- S2(A-B, gsea pol)  
- 2F S4B (housekeeping y predictibilidad)  
- 2G \+ S5 (puffin)  
- sure  
- xpresso  
- podria sumar (en alguna) tissue sp y ruido

Influencia de la multiplicidad de promotores por gen sobre las caracteristica transcripcionales

- 4A-B  
- caracterizacion (?)  
- S7A \+promoter similarity (o esto va a otra seccion? veremos que tan largo queda)  
- clasificacion: esquemas fig4C \+ S7B  
- caracterizacion (?)  
- 4C \+ ruido?

**Discusión**

# Metodos Exp

**Materiales y metodos experimentales:**  
Descripción general del enfoque

Con el fin de cuantificar las propiedades transcripcionales que se desprenden intrínsecamente de la secuencia del promotor basal, diseñamos un ensayo reportero masivo en paralelo (MPRA, por sus siglas en inglés) que permite evaluar simultáneamente la media y el ruido transcripcional, similar a las técnicas de sort-seq utilizadas principalmente en levaduras y bacterias. Utilizamos una *library* con 23908  secuencias de ADN de 252pb que, en su gran mayoría, comprenden un TSS anotado y su región flanqueante (235pb río arriba, 16pb río abajo), abarcando los elementos esenciales del promotor basal (7,8). La misma fue producida previamente por el grupo de la Dra. Fiszbein (Boston University), en cuyo laboratorio, y en colaboración con la Dra. Uriostegui-Arcos, realicé la totalidad de los experimentos de la presente tesis. Estas secuencias se clonaron en un plásmido reportero, río arriba de la secuencia codificante para EGFP. Mediante el uso de la técnica de intercambio de casete mediado por recombinasa (RMCE), generamos líneas estables de HEK293T, donde cada célula posee una copia única del constructo conteniendo al gen EGFP bajo control de un promotor de la *library*, todas en la misma posición del genoma. 

Para obtener la distribución de niveles de expresión de cada promotor, las células se fraccionaron en siete poblaciones según la señal de fluorescencia de EGFP, se extrajo su ADN genómico (ADNg), y se analizó la abundancia de cada promotor en las distintas poblaciones mediante Amplicon-seq. A partir de estos datos, reconstruimos la distribución de la expresión de cada promotor para estimar su nivel de expresión media y la varianza asociada. El procedimiento se esquematiza en la Figura M1.

El procedimiento se realizó enteramente a partir de dos poblaciones distintas de células, generando dos réplicas independientes. Asimismo, para preservar la representatividad de la *library* y evitar cuellos de botella, los pasos de biología molecular y celular se ejecutaron con múltiples réplicas técnicas y volúmenes de trabajo superiores a los estándares.

Fig. M1 (pipeline exp)

Construcción de los plásmidos reporteros conteniendo la *library* de promotores basales

Los oligonucleótidos de cadena simple que conforman la *library* (cuya composición y detalle se especifican en la Sección Métodos Computacionales), se amplificaron inicialmente mediante PCR de extensión por solapamiento (*overlap extension PCR*) usando Platinum SuperFi II Master Mix (ThermoFisher, 12368010\) y los *primers* PromLib Forward y Reverse (Anexo) durante 14 ciclos, seguidos de una purificación por columna (QIAGEN, 28104). 

El vector *backbone* (Fig. M2A) —término utilizado para este plásmido a lo largo de este trabajo— en el que se clonó la *library*, es una variante del reportero fluorescente bicromático (EGFP \+ DsRed) de *splicing* alternativo desarrollado por Orengo et al. (74). En esta variante, mutaciones específicas en el sitio de *splicing* 3' del exón alternativo aseguran la expresión exclusiva de EGFP. A su vez, el plásmido empleado comprende sitios LoxP no idénticos flanqueando la construcción reportera, que permiten su integración sitio-específica en células mediante un sistema de intercambio de casetes mediado por recombinasa (RMCE).

La elección de un vector derivado de un reportero de *splicing* alternativo obedeció a la disponibilidad y compatibilidad con líneas celulares desarrolladas previamente por nuestras colaboradoras (Fiszbein Lab, Boston University). En la Figura M2B se observan células HEK293T donde se evidencia la diferencia entre el reportero de *splicing* original y la versión optimizada para este ensayo, ambos bajo el promotor viral CMV. Además, la ausencia de expresión de dsRed fue validada previamente mediante RT-qPCR.

Para la integración de la *library* al plásmido *backbone,* se digirió a este último con las enzimas de restricción NdeI y NheI (ThermoFisher, FD0583 y FD0974), y se purificó el fragmento del vector lineal por electroforesis preparativa en gel de agarosa, utilizando un gel 1% y el kit QIAQuick Gel Extraction (QIAGEN, 28704). Con el objetivo de mejorar la eficiencia del clonado posterior, se utilizó el kit de purificación de productos de PCR QiaQuick (QIAGEN, 28104\) seguido por el kit de purificación MinElute (QIAGEN, 28004), que a su vez permite la concentración del fragmento en menor volúmen. La *library* se clonó en el plásmido digerido mediante ensamblaje Gibson (*Gibson assembly mix*, NEB, E2611L) en una relación molar inserto-vector de 3:1, con 0.250pmoles de ADN total. La reacción de ensamblaje Gibson se precipitó con isopropanol para remover las sales que puedan interferir con la electroporación posterior. Se introdujo el precipitado por electroporación en células MegaX DH10B T1R Electrocomp™ (ThermoFisher, C640003) para su transformación. Las bacterias fueron plaqueadas en placas cuadradas (245mm de largo) de LB ágar con ampicilina 100mg/L. Finalmente, los plásmidos se aislaron usando el kit de purificación de plásmidos PureLink™ Expi Endotoxin-Free Maxi (Invitrogen, A31217). 

Para comprobar la correcta inserción del promotor, se realizó un ensayo de *colony* PCR (directo de bacterias, sin purificación del ADN) utilizando los *primers* PromLib. No se detectó la presencia del promotor CMV original (\~700 pb), mientras que las bandas obtenidas corresponden al tamaño esperado de la *library* (\~300 pb) (Fig. M2C)

Fig M2. ESQUEMA VECTOR \+ microscopia \+ PCR

Cultivo celular y transfección

Se utilizaron las células HEK293T-A2 (34), que ya han sido caracterizadas y poseen un único locus para RMCE en su genoma. Las mismas se cultivaron en medio DMEM (Dulbecco’s Modified Eagle Medium) con alto contenido de glucosa y piruvato (Gibco, 11965118), suplementado con suero fetal bovino al 10% (Gibco, A31406-02). Las células se mantuvieron en incubadora humidificada a 37°C con 5% de CO2.

Para generar líneas celulares estables con integración genómica del reportero, la *library* se co-transfectó con un 10% (m/m) de un plásmido que codifica la recombinasa Cre. Se transfectaron un total de 15 µg de ADN plasmídico en ocho placas de 10 cm utilizando Lipofectamine 3000 (ThermoFisher, L3000-015) y Opti-MEM (Gibco, 31985-070), siguiendo las instrucciones del fabricante. Después de la transfección, las células se seleccionaron con 2 µg/mL de puromicina (Gibco, A1113802) durante dos semanas. Posteriormente, las células se distribuyeron en dos grupos distintos, los cuales se trataron como réplicas biológicas independientes. 

Para las transfecciones individuales de los plásmidos reporteros, se empleó el mismo procedimiento escalando proporcionalmente las cantidades: 2 µg de ADN plasmídico en un único well de una placa de seis wells. Esto incluye tanto los controles con el plásmido *backbone* con CMV, tanto en su versión bicromática como mutada, como aquellos con promotores específicos de la *library* utilizados para el control de consistencia. La confección de estos últimos implicó en primera instancia la amplificación de los promotores desde el *pool* de la *library* (en Anexo se presentan los *primers* correspondientes a cada uno), y su clonado por ensamblaje de Gibson en el plásmido *backbone* digerido.

Las células se visualizaron en un microscopio ECHO Revolve (RVL-100M) con un objetivo de 20X.

Citometría de flujo y *sorting*

Para el *sorting* de la *library* de promotores, se utilizó un equipo Beckman Coulter MoFlo Astrios de seis vías. Los datos de citometría de flujo se analizaron con FlowJo v10.10.0.

Para tener seguridad sobre el nivel de discretización suficiente para la distinción de los patrones de ruido y media que pudiéramos observar, se realizó un análisis de potencia considerando múltiples alternativas. Es importante comentar que los costos del ensayo crecen en gran medida con el aumento de las fracciones recolectadas, tanto por el tiempo de *sorting* creciente como los costos de secuenciación. Por lo tanto, la elección final resulta de un compromiso entre nuestra capacidad de medir las variables de interés (expresión media y dispersión) y el costo del experimento.

Fig M3. Potencia(hacer)

Se estableció que todas las fracciones tengan idéntico número de células para evitar sesgos diferenciales entre fracciones en el tratamiento posterior de las mismas. Para la determinación de los umbrales, se usó como base una citometría de flujo de las mismas células. La fracción \#6 sería aquella con menos células y por lo tanto es la que definiría el tiempo total (y con ello el costo) del *sorting*, y no poseería un límite superior. Se decidió por lo tanto que comprenda el 1% de la población total de células, estableciendo así la posición del umbral más alto. Respecto al menor, se lo ubicó a partir en del pico del control negativo de células HEK293T-A2 no transfectadas y de las celulas transfectadas(Fig. M4A). Una vez establecidos el límite inferior y el umbral más alto, se dividió el rango que queda en seis intervalos de tamaño constante en escala logarítmica (Fig. M4B). El diseño con intervalos equivalentes en tamaño es esencial para la correcta comparación entre fracciones. Como control positivo se emplearon células que portan el reportero bajo el promotor CMV. 

SPara la preparación del experimento en sí mismas, se analizaron las células HEK293T-A2 transfectadas con la *library* hasta recolectar un total de 250.000 células por gate, las cuales se preservaron en DNA/RNA Shield (Zymo Research, R1100-250) para los experimentos posteriores. Este procedimiento se realizó tres veces a partir de cultivos independientes de células transfectadas con la *library*. La decisión de recolectar idéntica cantidad de células para todas las fracciones, y que no sea proporcional al porcentaje de células en la población total, busca evitar sesgos en los procesos posteriores de extracción de ADN y PCR, principalmente. Eventualmente, sin embargo, esto llevará a distorsiones en las distribuciones reconstruidas de los promotores, que serán corregidas computacionalmente (ver Cuantificación de la actividad media y el ruido de los promotores de la *library)*

Fig M4. Figura gates

Controles de *spike-in* basados en células

Para evaluar posibles sesgos técnicos entre las distintas fracciones introducidas durante la manipulación de las muestras *post-sorting* (desde la extracción de ácidos nucleicos hasta la secuenciación), decidimos introducir moléculas identificables en cantidades conocidas (*spike-in*) a las mismas. La incorporación de *spike-in* es una práctica frecuente en análisis genómicos/transcriptómicos que implican la comparación cuantitativa precisa entre muestras procesadas por separado. Aunque el procedimiento estándar es el agregado de moléculas de ADN, ARN o células de otras especies, dadas las características del presente ensayo decidimos generar un control *spike-in* similar pero discernible al material con el que se iba a trabajar: células HEK293T-A2 con el reportero incorporado en el genoma. 

Brevemente, se clonaron tres secuencias de aproximadamente el mismo tamaño y contenido de GC que las de la *library* (ver Anexo) en el vector *backbone* y se co-transfectaron con el plásmido de expresión de la recombinasa Cre en células HEK293T-A2, para obtener líneas estables, en la forma descrita anteriormente. Se agregaron 30, 300 o 3000 células de cada línea de *spike-in* respectivamente en cada tubo de recolección de FACS, juntándolas con las células separadas por el *sorter*.

Co-extraccion de ADN/ARN  y RT-qPCR de EGFP

### 

El ADNg y el ARN total se co-extrajeron de las fracciones obtenidas utilizando el *Quick-DNA/RNA Miniprep Kit* (Zymo Research, D7001). Mientras que la extracción del ADNg es central para el ensayo, el ARN se utilizó para evaluar la correlación entre el *sorting* asociado a la señal de EGFP y los niveles de transcripto de dicho gen. Dado que el objetivo es vincular la secuencia promotora con su actividad transcripcional a través de la intensidad de fluorescencia, resulta crucial verificar la correspondencia entre la señal proteica y la abundancia del ARNm de EGFP.

El ADNc se sintetizó a partir de las muestras de ARN de cada *bin* de células utilizando el *cDNA Synthesis Kit* (Thermo Scientific, K1622) siguiendo las instrucciones del fabricante. La transcripción reversa se realizó con 500ng de ARN total. Posteriormente, se determinaron los niveles de ARN de EGFP mediante PCR cuantitativa (qPCR) con primers específicos (ver Anexo), utilizando la mezcla de reacción  *Maxima SYBR Green/ROX qPCR Master Mix (2X)* (*Thermo Scientific,* K0222) en un sistema *ABI 7900HT Fast Real-Time PCR* (Applied Biosystems). Los niveles de expresión génica relativa se calcularon mediante el método 

2−ΔCt 

(75), donde 

ΔCt \= (CtEGFP \- CtGAPDH)

, usando GAPDH como control interno (Fig. M5). Se realizaron tres réplicas técnicas para cada determinación, utilizando el promedio de las mismas para el cálculo de la abundancia.

Fig. M5 QPCR 

En primera instancia, los resultados (Fig. M5) indican que la fracción 0 (correspondiente al nivel mínimo de fluorescencia, apenas por sobre el de las células sin transfectar) no presentaba resultados consistentes. Este fenómeno es atribuible a la dificultad de discernir la señal de EGFP de la autofluorescencia propia de las células para células con niveles muy bajos de expresión. Por consiguiente, se decidió excluir esta población de la secuenciación. Asimismo, se descartó la réplica 3 debido a que no presentó un patrón de qPCR consistente (Fig. M5, derecha) y mostró un rendimiento insuficiente en la extracción de ADNg en varias de sus poblaciones (no mostrado).

Adicionalmente, se evaluó la posible persistencia de células no transfectadas tras la selección con puromicina. Si bien estas células no generarían amplicones de la *library*, su acumulación en los gates de menor expresión podría afectar la representatividad del número de células positivas. Para detectar la correcta presencia de los promotores de la *library*, se utilizaron *primers* flanqueantes a los sitios LoxP de las células HEK293T-A2 (XXXX, ver Anexo) para realizar una PCR sobre el ADNg de las distintas fracciones, permitiendo discernir por tamaño la presencia del constructo de la *library* (Fig. M6). Aunque se detectaron células sin el reportero, su abundancia relativa fue uniforme en todas las muestras, por lo que no se consideró un factor de sesgo para los análisis posteriores.

Fig. M6 IMAGEN PCR\_loxP\_gates\_ctl

Secuenciación de amplicones (Amp-seq)

Para la secuenciación paralela masiva del ADNg de las distintas fracciones, se empleó un enfoque de PCR de dos pasos. Primero, se realizó una PCR primaria usando los *primers* PromLib Forward y Reverse (Anexo) con *Platinum™ SuperFi II PCR Master Mix* (Thermo Fisher Scientific, Cat. No. 12-368-010), siguiendo las instrucciones del fabricante. Luego, se llevó a cabo una segunda PCR para incorporar las secuencias de los adaptadores de Illumina, manteniendo constante el *primer* *forward* con el índice i5 y usando un *primer* *reverse* con el índice i7 específico para cada muestra (Tabla S2), empleando nuevamente la *Platinum™ SuperFi II PCR Master Mix*. Los productos de la segunda PCR se purificaron mediante extracción de bandas del gel (QIAGEN, 28704\) y posteriormente se enviaron para la secuenciación de amplicones a MedGenome[^1]. Las secuenciacioónes se realizaronrealizó con una cobertura de 1000× utilizando lecturas *paired-end* de 150 pb en la plataforma *NovaSeq*. La incorporación de un índice i7 específico por muestra permite mezclar las muestras en una única corrida de secuenciación y evitar un posible *batch effect*. Sin embargo, cada réplica debió mantenerse en corridas en paralelo debido a la ausencia de más índices. 

Secuencias de interes y primers utilizados

| Nombre | Secuencia (5'→3') |
| :---- | :---- |
| PromLib\_Forward | TACGAAGTTATATGGATCCATATG |
| PromLib\_Reverse | TGGAAGCTTAAGTTTAAACGCTAG |
| TMEM87A\_1\_Forward | TACGAAGTTATATGGATCCATATGGGCCAGGCTGGCATGTAG |
| TMEM87A\_1\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGTTCACAGCCGTGGAGTG |
| PPP1R14B\_3\_Forward | TACGAAGTTATATGGATCCATATGCCCCACCCCCAGGGCCC |
| PPP1R14B\_3\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGGCCACGGGCCTGGAAGAC |
| ZKSCAN2\_1\_Forward | TACGAAGTTATATGGATCCATATGGGCAGGTTCCCTAGAATC |
| ZKSCAN2\_1\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGTCGGCCCGCGGAGAGCG |
| METAP2\_1\_Forward | TACGAAGTTATATGGATCCATATGTTGCTTCGGGAATGC |
| METAP2\_1\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGGAGAGCGCGAGGGAA |
| KIAA0753\_Forward | TACGAAGTTATATGGATCCATATGGCCACAACACGATGA |
| KIAA0753\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGCTGACAGAGCAAAAG |
| BTG1-Forward | TACGAAGTTATATGGATCCATATGGTCTCCAGCCGCCAC |
| BTG1-Reverse | TGGAAGCTTAAGTTTAAACGCTAGCCAGCTCCGAGAGGC |
| ETS1\_1\_Forward | TACGAAGTTATATGGATCCATATGCGCCAGCCCTTCCTTTCG |
| ETS1\_1\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGGGCGGCTGCCTCGTTCG |
| LSM1\_1\_Forward | TACGAAGTTATATGGATCCATATGAGGTGGGTGTACCGG |
| LSM1\_1\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGGGTTCGGCAGCAGAAGG |

Tabla1: Secuencia de los primers utilizados para los clonados de la *library* en el plásmido *backbone*, así como aquellos utilizados para amplificar los promotores específicos de la *library*. 

| Nombre | Secuencia (5'→3') |
| :---- | :---- |
| **Control de integración** |  |
| LoxP\_EF | CCAGCTTGGCACTTGATGT |
| LoxP\_WR | GGGCCACAACTCCTCATAAA |
| **Secuencias de Spike-in** |  |
| spikeIn\_SV40 | CATACACGGTGCCTGACTGGCGTTAGCAATTTAACTGTGTGATAAACTACCGCATTAAAGCTTTTTGCAAAAGCCTAGGCCTCCAAAAAAGCCTCCTCACTACTTCTGGAATAGCTCAGAGGCCGAGGCGGCCTCGGCCTCTGCATAAATAAAAAAATTAGTCAGCCATGGGGCGGAGAATGGGCGGAACTGGGCGGAGTT |
| spikeIn\_CMVe | CTAGGACAATTGATTATTGACTAGTTTATTAATAGTATCAATTACGGGGTCATTAGTTCATAGCCCATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGACTTTCCATTGACGTCAAT |
| spikeIn\_CMVeMut | CTAGGACAATTGATTATTGACTAGTTTATTAATAGTATCAATTACGGGGTCATTAGTTCATAGCCCATATATGGAGTTCCGCGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGACTTTCCATTGACGTCAAT |
| **qPCR primers** |  |
| GAPDH\_qPCR\_Forward | TACGAAGTTATATGGATCCATATGGAAAGAGTGACACCCCG |
| GAPDH\_qPCR\_Reverse | TGGAAGCTTAAGTTTAAACGCTAGCAGGGCTGTGGGTCCTGG |
| EGFP\_qPCR\_Forward | AAGTTCAGCGTGTCCGGC |
| EGFP\_qPCR\_Reverse | TCAGGGTGGTCACGAGGG |
| **Amplicon-seq** |  |
| TruSeq\_Universal\_primer\_Fwd (i7)  | AATGATACGGCGACCACCGAGATCTACAC |
| TruSeq10\_Rev\_Gate\_1 (i5+P5) | CAAGCAGAAGACGGCATACGAGATTAGCTT |
| TruSeq1\_Rev\_Gate\_2 (i5+P5) | CAAGCAGAAGACGGCATACGAGATATCACG |
| TruSeq23\_Rev\_Gate\_3 (i5+P5) | CAAGCAGAAGACGGCATACGAGATATGAGTGG |
| TruSeq13\_Rev\_Gate\_4 (i5+P5) | CAAGCAGAAGACGGCATACGAGATTGAGTCAA |
| TruSeq6\_Rev\_Gate\_5 (i5+P5) | CAAGCAGAAGACGGCATACGAGATGCCAAT |
| TruSeq14\_Rev\_Gate\_6 (i5+P5) | CAAGCAGAAGACGGCATACGAGATACAGTTCC |

Tabla 2: Secuencia completa de los Spike-in generados e insertados, de los *primers* utilizados para el control de células con el plásmido integrado, para la qPCR y los adaptadores utilizados durante la secuenciación.

# Metodos bioinfo

**Materiales y métodos computacionales:**

**Procesamiento de datos y cuantificación de la actividad promotora**  
   
Procesamiento de los datos de secuenciación

Las lecturas *paired-end* de *Illumina* se procesaron inicialmente con *cutadapt* (56) para el filtrado por calidad y el recorte de *primers* y de calidad, y posteriormente con *FASTP* para el recorte de colas de poli-G y de baja calidad (77) (Fig. M7). Se descartaron aquellas lecturas que carecían de superposición de *primers*. El alineamiento contra las secuencias *FASTA* de la *library* se realizó mediante *HISAT2* (78). Únicamente se consideraron fragmentos con longitudes entre 230 y 270 pb, omitiendo los alineamientos secundarios. La abundancia de cada promotor en cada muestra se cuantificó con *SAMTOOLS* (1.20) (79). El control de calidad de las lecturas crudas se llevó a cabo con *FASTQC[^2]* y el desempeño del *pipeline* se reportó mediante *multiQC* (80).   
Fig. M7 pipeline bioinfo

Cuantificación de la actividad media y el ruido de los promotores de la *library*

Salvo que se indique lo contrario, los análisis estadísticos y el análisis exploratorio de datos se realizaron utilizando R 4.2.  
Los conteos de lecturas por promotor se corrigieron según la desviación del tamaño de cada *library* respecto al tamaño medio. Además, dado que el esfuerzo de muestreo para obtener 250.000 células varía para cada rango de intensidades de EGFP, los conteos también se corrigieron mediante un término de relativización determinado por:

Fk=pkmin(p)										(E1)

donde [![][image1]](https://www.codecogs.com/eqnedit.php?latex=p#0) es la proporción de células en la fracción [![][image2]](https://www.codecogs.com/eqnedit.php?latex=k#0) respecto a la población total, medida por citometría de flujo (Fig M8A). Es decir que las fracciones que, al medir por citometría de flujo, presentan mayor número de células, se les amplificará, en términos relativos y proporcionales, el número de *counts* por promotor observados. Esta corrección, más típica de la ecología que del campo de la biología molecular, es fundamental para comparar correctamente entre las fracciones de células y evitar un sesgo hacia los gates más altos, donde hay menos células en la población total. 

Se calcularon estadísticas descriptivas, como la media y la varianza, para cada secuencia que superara un umbral de 1000 conteos corregidos en cada réplica, asignando a cada fracción un valor de expresión coincidente su índice numérico (1 al 6). Las secuencias con una diferencia en la mediana superior a 2 entre ambas réplicas se descartaron por inconsistencia. Asimismo, se descartaron las secuencias con patrones altamente bimodales al considerarlas potenciales artefactos (Fig. M8B). Para detectar estos patrones, se ideó la estrategia de considerar aquellas con una varianza un 20% superior a la esperada para una distribución uniforme, dado un número [![][image3]](https://www.codecogs.com/eqnedit.php?latex=j#0) de *gates* con conteos. La varianza esperada de la distribución uniforme se calculó como:

Si [![][image4]](https://www.codecogs.com/eqnedit.php?latex=j#0) es impar:

[![][image5]](https://www.codecogs.com/eqnedit.php?latex=%5Ctext%7Bvar%7D\(j\)%20%3D%202%20%5Csum_%7Bi%3D1%7D%5E%7Bj%2F2%7D%20\(i%20-%200.5\)%5E2#0)

Si [![][image6]](https://www.codecogs.com/eqnedit.php?latex=j#0) es par:

[![][image7]](https://www.codecogs.com/eqnedit.php?latex=%5Ctext%7Bvar%7D\(j\)%20%3D%202%20%5Csum_%7Bi%3D1%7D%5E%7B\(j-1\)%2F2%7D%20\(i%20-%200.5\)%5E2#0)								(E2)

Fig M8. sample effort \+ bimodal

**Caracterización de la *library***

Composición y diseño de lLa *library*

La *library* utilizada consiste en 23908 secuencias de 300 pb de longitud, de las cuales 20851 incluyen una región de promotor basal (core promoter) de \-251 a \+16 con respecto a la posición anotada como TSS principal en la Eukaryotic Promoter Database[^3] (EPD) (73), con adaptadores de 24pb en cada extremo para facilitar la clonación (Fig. M9 A-B). EPD define sus TSS (29512 sitios en total) a partir de datos del repositorio del consorcio FANTOM que utilizan la técnica *Cap Analysis of Gene Expression* (CAGE) para la determinación del extremo 5’ de los transcriptos con la precisión de un nucleótido y con ello inferir la actividad promotora a lo largo de todo el genoma en diversas muestras biológicas. Agrupando únicamente las muestras de dicho repositorio correspondientes a cultivos primarios de células, se observa como el TSS indicado por EPD para los promotores seleccionados de la *library* coincide con aquel más usado, si bien son mínimos los casos donde no se observa actividad proveniente de las bases vecinas (Fig. M9C). Es preciso aclarar que en ningún momento se evaluó el TSS en el contexto del reportero utilizado, por lo que cada vez que me refiera a dicho sitio a lo largo de esta tesis, será teniendo en cuenta al anotado.

De las restantes secuencias presentes en la *library*, 2910 son regiones de enhancer de 152pb de longitud inmediatamente río arriba de 100pb del promotor FN1(Fig. M9 A-B). También hay un subgrupo de 147 promotores asociados al cáncer, tanto en sus versiones wild-type como mutantes, y que no responden a la estructura tal cual fue definida para los promotores EPD. Estos últimos grupos, aunque presentes en la *library* y secuenciados, no fueron considerados durante el análisis bioinformático posterior al recuento de actividad.

Fig. M9   
   
Elementos de secuencia del promotor basalCaracterización de los promotores

La base de datos EPD provee información respecto a la presencia de motivos típicos en sus promotores anotados, a partir de la búsqueda de patrones y su localización respecto al TSS: TATA-box, CCAAT-box, GC-box e INR. Según sus datos, del 23851 secuencias, el 47.7% contiene un GC-box, el 32% un INR fuerte, el 16.3% un CCAAT-box y solo el 8% TATA-box (Fig. M10A).

Fig. M10

La presencia de islas CpG fue determinada a partir de las anotaciones provenientes de UCSC, considerando aquellos casos donde se superponen más del 100pb con la secuencia promotora en cuestión. Las islas CpG son bastante prevalentes en promotores, observándose en el 64% de los mismos (Fig. M10B). También se buscó el solapamiento con elementos repetitivos y/o transponibles, anotados en la base de datos RepeatMasker[^4]. Además de las clases de transposones presentes en dicho repositorio, se incorporó una clase denominada “Repeticiones de Baja Complejidad”, agrupando las categorías *Simple Repeats* y *Low Complexity*, asi como una de “Elementos transponibles”. 

La frecuencia nucleotídica y su identidad en sitios específicos fue evaluada con el paquete de R-Bioconductor Biostrings[^5]. Se evaluó el contenido de G y C (Fig. M10C), así como la identidad del dinucleótido del TSS (Fig. XD). Los patrones buscados de forma estricta fueron YCASW para el INR “fuerte”, TCT para el clásico motivo de proteínas ribosomales y el dinucleótido YR (PyPu)[^6]. Respecto a este último patrón se discernir también en sus cuatro posibilidades (CA, CG, TA, TG). Aquellos promotores que no cuadran en su TSS con alguno de los patrones mencionados, fueron catalogados como “No canónicos”. Se incluye el dinucleótido GC para evidenciar que la prevalencia del CG no es simplemente producto de alto contenido de dichos nucleótidos.  

Patrones de conservación de los promotores basales

A su vez, se incorporaron datos que reflejan la historia evolutiva de los promotores, tanto a nivel de secuencia como funcional. Por un lado, se extrajeron datos de PhyloP score (81) provenientes de la comparación entre el genoma hg38 y 100 especies de mamíferos. Evaluando la mediana de dicho valor del “metapromotor” a cada base (Fig. M11A), se observa claramente una mayor conservación en la región proximal al TSS (+16 a \-50), con claros picos alrededor del \+1 y del \-30, asociado al TATA-box. En una región de intermedia cercanía (-50 a \-150) hay un progresivo decaimiento de la conservación, mientras que se acerca mucho a valores de evolución neutra para la región más distal (-150 a \-235). A su vez, para cada una de estas regiones, en cada promotor, se evaluó el PhyloP score medio. En términos de conservación funcional, nos basamos en datos de Young et al. (19), quienes utilizan datos de actividad promotora en tejidos de humano y ratón para considerar si, las secuencias que se pueden considerar homólogas, están activas en ambas especies, si perdieron actividad promotora en humanos o en ratón o si, por el contrario, la adquirieron en alguna de estas especies (Fig M11B). 

Fig. M11. evo library

Patrones de actividad endógena de los promotores

Se obtuvieron datos de accesibilidad de cromatina en HEK293 provenientes de ENCODE (accesible como ENCSR956YZJ). La presencia de Módulos Regulatorios en *cis* (CRM) se determino con datos de Remap basados en cientos de muestras de ChIP-seq [(44, ver](https://www.zotero.org/google-docs/?sWJeqi) Análisis de datos de ChIP-seq[)](https://www.zotero.org/google-docs/?sWJeqi). La posición de los enhancers anotados se obtuvo de UCSC. 

La actividad promotora endógena fue cuantificada a partir de datos de CAGE, provistos por el consorcio FANTOM5. Estos datos se utilizaron tanto para la determinación de la actividad en muestras específicas, así como para obtener valores del comportamiento de los promotores a lo largo de todas las muestras, como la especificidad de tejido del promotor o su forma (la distribución de sus TSS). Para las muestras específicas (HEK293, HeLa, músculo esquelético), se utilizó el paquete de R CAGEr v2.12.0 [(83)](https://www.zotero.org/google-docs/?OtTt7L). Específicamente, aplicamos una normalización de *power-law* para calcular los niveles de expresión en *Tags* Por Millón (TPM), filtrando aquellas señales que estuvieran por debajo de un umbral de 1 TPM. Los clusters de TSS se definieron mediante el algoritmo *paraclu* (84). Para refinar la señal, se conservó el 80% central del rango de expresión y se intersectó con las coordenadas genómicas de la *library* para evaluar la actividad endógena específica de cada promotor.

Si bien esta metodología ofrece una mayor precisión, su alta demanda computacional (en términos de tiempo y memoria) la hizo inviable para evaluar la actividad de los promotores en cada una de las muestras de células primarias y tejidos del dataset del consorcio FANTOM5. Dado que la representación de tipos celulares y tejidos en el dataset es muy heterogénea, decidimos agruparlos por ontología utilizando las categorías definidas por Andersson et al. 2014 (85). Las muestras se fusionaron, se normalizaron por el tamaño de la *library* y se restringieron a los rangos genómicos de la *library*; finalmente, la actividad de los promotores en cada grupo se evaluó como la suma de los *counts* que solapaban con el rango genómico de cada promotor. Cabe aclarar que, si bien la base de datos de FANTOM5 incluye otras categorías de muestras humanas, como por ejemplo líneas celulares,  Para distinguir entre promotores *housekeeping* y aquellos con alta especificidad de tejido, utilizamos el índice de Gini de la actividad del promotor a lo largo de todos los grupos de ontología y dividimos a los promotores en terciles. Los promotores inactivos se excluyeron de este análisis (menos de 1 TPM en cualquier grupo) y se etiquetaron como "Sin actividad detectable". La fórmula para el índice de Gini, luego de ordenar cada uno de los *n* grupos de ontología por actividad decreciente, es:  
G=1-i=1n(pi-pi-1)(qi+qi-1)							(E3)  
donde, en este caso, *pi* y *qi* son las proporciones acumuladas de los grupos de muestras y de la actividad promotora, respectivamente.

La forma de los promotores se calculó a partir de la combinación de todas las muestras de FANTOM5 analizadas (células primarias y tejidos). En este caso, también se utilizó el paquete *CAGEr* y se aplicó una normalización de *power-law*. Los *clusters* de TSS se obtuvieron mediante el método *distclu*, con una distancia máxima entre TSSs de 5 y un mínimo de 10 *counts* por TSS. Se utilizó el ancho intercuantílico (0.05-0.95) como *proxy* de la forma. Los valores se dividieron en terciles de ancho, conservando el primer y tercer tercil como promotores focalizados y anchos, respectivamente. En caso de que se obtuvieran múltiples *clusters* de TSS para una misma región promotora, solo se utilizó aquel que solapaba con la posición establecida por la EPD.

Coocurrencia de las características de los promotores basales

Las características de los promotores, tanto aquellas basadas en la secuencia como aquellas que surgen de estudiarlos en sus contextos endógenos, no son completamente independientes entre sí (Fig. M12). Esto implica que frecuentemente, sea complejo poder asignar a un efecto observado, una característica particular, con confianza de que no se trate de un efecto confusor de otra característica con alto grado de coocurrencia. Si bien esto se puede resolver en ciertos casos con una estratificación por a potencial característica confusora, en los casos mas extremos de co-presencia, esta tarea resulta prácticamente imposible y es un limitante en este tipo de enfoques experimentales, basados en secuencias naturales.  

Fig. M12 coocurrencia

Algoritmos de predicción de la actividad promotora basados en la secuencia

Asociación de las características del promotor con la actividad y ruido transcripcional

Dado que las diferencias en la representación de las secuencias entre los conjuntos de datos de ambas réplicas pueden complicar el análisis si se juntasen sus resultados, decidimos realizar el análisis de asociación de características para cada réplica por separado y conservar únicamente aquellas que estuvieran significativamente asociadas en ambas. A su vez, se dejaron de lado aquellas secuencias que no correspondieran a promotores anotados en la base EPD.  
Para evaluar la asociación entre cada característica y la actividad de los promotores evaluados anotados en EPDl promotor, se realizó un test de Wilcoxon utilizando el paquete de R Coin v1.4.3 (86). Para establecer la significancia, la réplica se consideró como una variable de efectos aleatorios. Sin embargo, para obtener el tamaño del efecto, realizamos el test en cada réplica por separado e informamos ambos valores. A su vez, para hacer foco en la replicabilidad, se limitaron los resultados a aquellos casos donde los efectos fueran consistentes en sentido entre ambas réplicas.   
En todos los casos donde se realizaron múltiples comparaciones, se aplicó la corrección de Benjamini-Hochberg sobre los *p-values* (87). 

Asociación de las características del promotor con el ruido transcripcional

La estimación del ruido y su asociación con las características de los promotores, tiene una complejidad extra. Existe una asociación natural, y también experimental, entre la media y la varianza, que deseamos desacoplar. Para ello, La magnitud del efecto de una característica sobre el ruido del promotor se estimó mediante un enfoque de rendimiento de una variable binaria, utilizando curvas del tipo *Receiver Operating Characteristics* (ROC) (32). En ese sentido, se buscó desacoplar la asociación entre la media y la varianza y para ello los promotores se agruparon por su media en diferentes bins y, dentro de cada uno, se asignó un *rank* a la varianza. Este último valor se utilizó como la variable continua para separar dos grupos de promotores: aquellos con y sin la presencia de una característica específica. La calidad de esta clasificación, y en definitiva la métrica de relevancia de la característica  para el ruido,se midió mediante el área bajo la curva de la curva ROC (*Receiver Operating Characteristics (32))* (AUC-ROC), con la ayuda del paquete de R pROC v1.18.5 (88). Ello implica evaluar para cada valor del ranking de varianza (1-100) cual es la proporción de promotores a cada lado de dicho umbral y según la presencia de la característica del promotor que se desea evaluar. Así se construye la curva ROC (sensibilidad en función de 1-especificidad) y finalmente la estimación del área bajo su curva.  A su vezAdemás, se realizó un bootstrapping de 1000 remuestreos para obtener intervalos de confianza del 95% sobre el AUC estimado. Las características cuyo intervalo de confianza incluye el valor 0.5, aquel propio de una distribución al azar, en alguna réplica se descartaron por no presentar un efecto significativo.

En sistemas biológicos, la media y la varianza de la expresión génica están intrínsecamente correlacionadas: los promotores con mayor expresión media tienden también a exhibir mayor varianza absoluta. Esta relación (conocida como ruido proporcional o efecto de Fano) impide comparar directamente el nivel de ruido entre promotores sin antes controlar por su nivel de expresión. Si no se desacopla esta asociación, cualquier característica que influya sobre la media aparecerá artificialmente como moduladora del ruido, generando asociaciones espurias.

#### **Estrategia de desacoplamiento: ranking intra-bin**

Para eliminar este confundidor, se adoptó el siguiente procedimiento:

1. **Agrupamiento por media (*binning*):** los promotores se ordenaron según su expresión media y se distribuyeron en *bins* de 100 promotores cada uno. Dentro de cada bin, los promotores comparten un rango similar de expresión media, de modo que las diferencias de varianza observadas no pueden atribuirse a diferencias en la media.  
2. **Asignación de *rank* de varianza intra-bin:** dentro de cada bin, los promotores se ordenaron por su varianza y se les asignó un valor de *rank* (rango percentil). Este ranking relativo —no el valor absoluto de varianza— es la variable continua resultante del desacoplamiento. Valores altos de rank indican promotores con mayor ruido *de lo esperado para su nivel de expresión*; valores bajos indican promotores más silenciosos de lo esperado.

#### **Evaluación de la relevancia de características: AUC-ROC**

Una vez obtenido el rank de varianza desacoplado, se evaluó si las características de los promotores (secuencia, estructura, factores de transcripción, etc.) se asocian con este ruido residual. Para ello se empleó el área bajo la curva ROC (AUC-ROC) (32), una métrica que mide la capacidad discriminatoria de una variable continua para separar dos grupos, de forma independiente al umbral de clasificación elegido.

El procedimiento fue el siguiente:

* Para cada característica de interés, los promotores se dividieron en dos grupos: aquellos que poseen la característica (*presencia*) y aquellos que no la poseen (*ausencia*).  
* El rank de varianza intra-bin actuó como variable de puntuación (*score*) para clasificar los promotores entre los dos grupos.  
* Se construyó la curva ROC evaluando, para cada posible umbral del rank (valores 1 a 100), la sensibilidad (proporción de promotores con la característica que superan el umbral) y la especificidad (proporción de promotores sin la característica que no lo superan). La curva ROC representa la sensibilidad en función de 1 − especificidad al barrer todos los umbrales posibles.  
* El AUC-ROC resume esta curva en un único valor: 0.5 indica que la característica no discrimina mejor que el azar; valores superiores a 0.5 indican que la presencia de la característica se asocia con mayor ruido relativo; valores inferiores a 0.5 indican asociación con menor ruido relativo. Los análisis se realizaron con el paquete `pROC` v1.18.5 en R.

#### **Estimación de la incertidumbre: *bootstrapping***

Para cuantificar la incertidumbre del AUC estimado y descartar asociaciones no significativas, se realizaron 1.000 remuestreos con reposición (*bootstrapping*). En cada remuestreo se recalculó el AUC, generando una distribución empírica que permitió construir intervalos de confianza del 95%.

Las características cuyo intervalo de confianza incluye el valor 0.5 en alguna de las réplicas fueron descartadas, ya que no pueden distinguirse de una clasificación aleatoria y, por tanto, no presentan un efecto significativo sobre el ruido transcripcional.

Análisis de datos de ChIP-seq

Se obtuvieron picos de ChIP-seq no redundantes de la base de datos ReMap y se intersectaron con la *library* de promotores para determinar el estado de unión de cada proteína. La intersección positiva de al menos un pico, fue evidencia suficiente para considerar la presencia de un CRM. El análisis se restringió a los factores de transcripción (TFs) que presentaron picos en al menos 100 promotores en ambas réplicas biológicas.  
Las asociaciones entre la ocupación proteica y tanto el ruido como la actividad transcripcional se evaluaron utilizando el mismo marco de trabajo aplicado a todas las características binarias. De manera paralela, se llevó a cabo un análisis para modificaciones de histonas y marcas epigenéticas obtenidas del conjunto de datos de ChIP-Atlas (58).

Análisis de enriquecimiento funcional

Se realizó un análisis de enriquecimiento de conjuntos de genes (*GSEA*) (89) sobre el subconjunto de factores de transcripción identificados utilizando términos de *Gene Ontology* (*GO*), ordenados ya sea por el tamaño del efecto del test de suma de rangos de *Wilcoxon* (*Wilcoxon rank-sum test*) o por el *AUC-ROC*.Estos análisis se implementaron mediante *scripts* propios utilizando los paquetes de *R/Bioconductor* *AnnotationDbi*, *clusterProfiler* (90) y *enrichplot*.

Clasificación de promotores alternativos

Para clasificar los promotores de acuerdo a su relación con otros regulando el mismo gen, unimos todos los promotores anotados en la base EPD junto con los datos de CAGE de FANTOM5 (tejidos y cultivos primarios).  Como ya se mencionó, las muestras de FANTOM5 fueron unificadas en base a su ontología para evitar redundancias y normalizadas por el tamaño de la *library*. En primer lugar, se clasificaron como “Sin actividad detectada” a aquellos promotores que no alcanzaran 1 TPM en ninguna muestra. Si bien podría llamar la atención que dichos promotores estuvieran incluidos en la base de datos de EPD, que también está basada en datos de FANTOM5, esto asumo que se debe a que son promotores que se registraron activos en líneas celulares, que estoy excluyendo del presente análisis. Se agruparon los promotores restantes por gen, y se separaron aquellos ”Promotores únicos”, para los genes sin promotores alternativos con actividad detectable.   
Para aquellos genes con múltiples promotores, se buscaron aquellos con mayor actividad para clasificar como “Promotores principales”. Para ello se consideraron aquellos promotores que:  
1- Concentran el máximo porcentaje de counts normalizadas sumando todas las muestras.   
2- Contienen el máximo número de muestras con la actividad más alta.   
Un pequeño grupo de promotores cumplió únicamente una de aquellas condiciones, que fueron considerados “No clasificables”. El resto se clasificaron como “Promotores secundarios”.  
En segunda instancia cada par Principal-Secundario fue clasificado en base a la correlación en sus actividades endógenas. Para ello, se inició por la identificación de casos de “Alternancia” (o *switch*), término con el que nos referimos a muestras donde se observa una clara alternancia en el rol Principal/Secundario de los promotores. Se definió como tal cuando el Secundario supera los 5 TPM y, o bien la relación de actividad en escala logarítmica fuera 50% mayor para el Secundario, o bien el Principal no tuviese actividad detectable en la muestra (\<1 TPM). Los pares de promotores sin un caso de alternancia, se clasificaron como “Correlacionados” o “Independientes” según si la correlación (o bien de Spearman o de Pearson) fuese mayor a 0,5. 

# Resultados 1

Estimación masiva de las propiedades transcripcionales de promotores basales humanos

Una vez obtenidos los datos de secuenciación, en sus dos réplicas, el primer paso fue el alineamiento de las lecturas obtenidas con las secuencias de la *library*. Como se ha mencionado previamente, una preocupación a lo largo del proceso experimental fue evitar un cuello de botella que implique que, al analizar los datos, las lecturas observadas correspondieran todas a un subgrupo pequeño de secuencias de la *library*. Afortunadamente, este no fue el caso, y contamos con lecturas, en al menos una réplica, del 80,9% de las secuencias, y el 67,3% en ambas réplicas (Fig. R1 venn) . 

Fig. R1.1 venn  rep1-rep2-tot.library

Si bien los análisis a realizar posteriormente serán todos a nivel comparativo entre las secuencias de las que podemos extraer datos confiablemente, nos planteamos la posibilidad de que haya, más allá de un cuello de botella aleatorio, un sesgo en la representatividad de las secuencias en las células. Comprender este aspecto nos serviría tanto para poder dimensionar las limitaciones en la extrapolación de los datos a la totalidad de los promotores anotados en humanos, así como para identificar puntos clave a mejorar en la metodología experimental a futuro. Con este objetivo, tomamos datos de un experimento previo realizado en el laboratorio de la Dra, Fiszbein, utilizando la misma *library,* con los adaptadores ya ligados, y en el que el proceso fue idéntico al ya expuesto, generando células HEK293T-A2 con una variante del reportero integrado. En dicho ensayo, se realizó una secuenciación de los promotores en dichas células previo a cualquier *sorting* celular o filtro por expresión. Encontramos una fuerte asociación entre los promotores seleccionados en ambos ensayos (OR \= 13.7, p \< 2.2e-16, *Fisher’s exact test*) (Fig. R2A venn\_presort). Esto indicaría que la incorporación de una secuencia a las células en el ensayo previo aumenta drásticamente la probabilidad de ser detectado en el ensayo actual, sugiriendo un sesgo común. Ante esto, surgen dos hipótesis: o bien hay un diferencia de partida en la representatividad de las secuencias, o alguno/s de los pasos experimentales presentan un sesgo sistemático.

Una posibilidad en este último sentido es el contenido de G y C en la secuencia, que podría haber generado pequeños cambios de eficiencia en la amplificación de los promotores por PCR, por ejemplo. Evidentemente, las secuencias no observadas en nuestros resultados, tienen valores de contenido G+C más extremos, en ambos sentidos (Fig. R2B gc\_bias\_violin). A su vez, este patrón se repite si se tiene en cuenta la distribución en el contenido de G+C en todas las lecturas: al comparar con la distribución hipotética e ideal en la que todos los promotores tuvieran igual representatividad en los resultados, se evidencia una depleción en promotores de alto y bajo contenido de G+C (Fig. R2C gc\_bias\_violin). Si bien en las amplificaciones por PCR se utilizó un *kit* que ha sido probado como eficiente para proporciones extremas de AT y GC en los amplicones, las eficiencias para dichos casos podrían diferir, evidenciando, luego de muchos procesos amplificadores, el patrón que observamos. Una posible solución para ensayos futuros podría ser la combinación de productos de PCR con protocolos optimizados para distintas proporciones de AT/GC. De cualquier manera, este sesgo, si bien limita el universo que podrá ser abordado en el presente trabajo, no afecta los resultados internos del mismo.

Fig. R2 Sesgo de representatividad

Dado que una parte esencial del análisis implica la comparación precisa entre las muestras correspondientes a los distintos *bins* de células por expresión, realizamos un control por *spike-in* celulares (ver Métodos). Evaluamos la abundancia de los *counts* (normalizando únicamente por *library size*) de las tres secuencias incorporadas post-*sorting* (Fig. R3). Afortunadamente, encontramos una relación relativamente constante para cada una de ellas, sin grandes desvíos en ninguna muestra que nos hubieran hecho plantearnos la posibilidad de reescalar los datos o tener que descartar alguna muestra.

Fig. R3. Spike-in

A partir de los datos de secuenciación obtenidos, en sus dos réplicas, y luego de los ya mencionados controles y procesamiento de los datos, se logró la reconstrucción de las distribuciones subyacentes de expresión asociadas a cada uno de las secuencias regulatorias. Esto se traduce visualmente en histogramas indicando la frecuencia de lecturas normalizadas con que se observó un cierto promotor en cada fracción de células (Fig. R4A \- aka 1B histo). Fue posible con estos datos obtener la media y la dispersión de los conteos, indicadores de actividad promotora y ruido transcripcional respectivamente. Para evaluar el grado de concordancia entre las distribuciones discretas reconstruidas y las distribuciones subyacentes, se seleccionaron una serie de secuencias de la *library* para replicar el procedimiento experimental pero con promotores individuales aislados. En estos casos, se generaron líneas celulares estables con el reportero integrado y regulado por cada una de dichas secuencias, a las que se le midió la fluorescencia de EGFP en cada célula por citometría de flujo (Fig. R4B \- aka 1C densities), con tres réplicas técnicas. Esto permite comparar la correlación entre la media calculada a partir de los datos discretizados y los continuos, observando la esperada asociación positiva (Fig. R4C \- aka 1D). A nivel ruido, las métricas no son lo suficientemente robustas como para permitir una buena determinación a nivel individual de los promotores y replicar esta prueba de consistencia.

Fig. R4 \- Distribuciones

La correlación entre las medias obtenidas para cada promotor en ambas réplicas tiene un coeficiente de Pearson de 0.72 (Fig. R5A mean\_replicates\_prefilter.jpg), contabilizando solamente aquellas secuencias con un umbral mínimo de *counts* y presentes en ambas réplicas, reflejando la robustez técnica del experimento realizado. Una vez incorporados los filtros de consistencia a nivel mediana y de exclusión de los patrones bimodales, que serán tenidos en cuenta para los análisis subsiguientes, se llega a un coeficiente de 0.85  (Fig. R5B mean\_replicates\_postfilter.jpg). En este último subconjunto de secuencias, de alta confianza, se encuentran los 12214 promotores que serán utilizados en los siguientes análisis, presentándose 10401 de ellos en ambas réplicas (Fig. R5B fig venn o bar post filter). De las distribuciones observadas en la Fig. 5A-B, y en mas detalle en Fig.5D, se desprende que las poblaciones de ambas réplicas presentan distribuciones distintas a lo largo de las subpoblaciones predefinidas. La distribución de un promotor es en última instancia una relación respecto a la distribución de todos los promotores de esa réplica, y sus valores de media y ruido obtenidos no tienen mucho sentido sin dicho contexto. Es por ello que se consideró, para la evaluación de efectos estadísticos en los análisis subsiguientes, a cada réplica por separado y de forma independiente. 

Fig. R5 Filtros

En resumen, el MPRA realizado permitió una cuantificación robusta de las propiedades transcripcionales intrínsecas de gran parte de los promotores basales humanos conocidos, ofreciendo la posibilidad de asociarlas con las características arquitectónicas de los promotores.

# Resultados 2

Efectos del promotor sobre la fuerza transcripcional

En primer lugar, analizamos la influencia de elementos conocidos del promotor basal sobre la fuerza transcripcional. Mientras que la importancia de estos motivos en el funcionamiento del promotor están bien establecidos, el campo todavía escasea de análisis masivos dirigidos y focalizados en los efectos de la secuencia sobre los niveles transcripcionales. Para evaluar el impacto de una dada característica sobre la actividad promotora, las secuencias evaluadas fueron divididas en grupos de acuerdo a la media de su distribución de expresión, calculando en cada uno la proporción de promotores conteniendo dicha característica. A su vez, para evaluar el tamaño del efecto y la significancia estadística de cada característica sobre la actividad media, se utilizó el estadístico asociado al test no paramétrico de Wilcoxon y su valor *p* (ver Métodos). Las islas CpG, por ejemplo, han sido asociadas con alta actividad promotora y autonomía (36-38), probablemente por la presencia de motivos ricos en GC para factores de transcripción (como las GC-box, (39-40)) y efectos sobre el posicionamiento nucleosomal (41). En concordancia con ello, encontramos un efecto positivo de este carácter, que se evidencia por el mayor número de promotores que superponen con islas CpG a crecientes niveles de actividad media (Fig. R2.1A)

Otro elemento del promotor muy bien estudiado aunque menos abundante es la TATA-box, un motivo posicionado con precisión que se conoce como un determinante positivo de una alta actividad promotora (42). De hecho, nuestro análisis muestra una fuerte asociación entre la presencia de una TATA-box canónica y una elevada actividad del promotor, explicada por la casi exclusiva presencia del motivo en los grupos de promotores más activos (Fig. R2.1B).

Fig. R2.1- TATA y CGI

Utilizamos este enfoque general para estudiar la influencia de otras características de la secuencia sobre la actividad promotora (Fig. R2. 2). Observamos asociaciones positivas con otros motivos promotores bien estudiados, como las GC-boxes y las CCAAT-boxes (Fig. R2.2, barras naranjas), en línea con sus roles conocidos en la regulación transcripcional (43, 44).

Fig R2.2- summary

Para evaluar si el efecto de estos motivos estaba relacionado con el reclutamiento de sus factores de transcripción asociados, utilizamos datos de ChIP-seq de la base de datos ReMap para analizar la asociación entre la fuerza del promotor y la presencia de picos de ChIP-seq en los promotores endógenos (Fig. R2.3). La asociación positiva de las CCAAT-boxes y GC-boxes fue reforzada por asociaciones igualmente positivas de los picos de ChIP-seq de NFYA (Fig. R2.3A) y *SP1/2* (Fig. R2.3B), respectivamente. Por el contrario, los promotores que solapan con elementos móviles o que exhiben TSS no canónicos tienden a tener una menor actividad dirigida por la secuencia (Fig. R2.2, barras verdes). Esto podría reflejar un estado no óptimo de los promotores más recientes o de promotores que no han estado bajo una fuerte presión selectiva.

Fig R2.3-TFs nfya y sp

Para profundizar en la relación mecanística entre la secuencia del promotor *core*, la actividad y la unión de TF, evaluamos cómo la actividad promotora se relaciona con el reclutamiento de más de 900 factores de transcripción (TF) en el *locus* del promotor endógeno, medido a través de datos de ChIP-seq de la base de datos ReMap (45) (Fig. S3).

La mayoría de los TF presentan una asociación positiva con la actividad del promotor (Fig. S3A, barras naranjas) y, al realizar un GSEA (*Gene Set Enrichment Analysis*) con estos TF *rankeados* según su efecto en la actividad (Tabla S1), hallamos que los términos enriquecidos más significativamente en promotores de alta actividad se vinculaban con la actividad de los factores de iniciación general de la RNA Polimerasa II (Fig. S3B). Este hallazgo brinda evidencia adicional de que la capacidad de las regiones fuertes del promotor *core* para ensamblar el complejo de preiniciación y luego iniciar la transcripción está *hard-wired* en la secuencia de DNA. Sorprendentemente, solo la 5-metilcitosina y tres TF, que incluyeron a los factores asociados a Polycomb CBX7 y JARID2, se asociaron negativamente con la actividad (Fig. S3A, barras verdes, recuadro).

[^1]:  https://diagnostics.medgenome.com/research-service/

[^2]:  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[^3]:  https://epd.expasy.org/epd/

[^4]:  https://www.repeatmasker.org/

[^5]:  https://bioconductor.org/packages/release/bioc/html/Biostrings.html

[^6]:  Las bases subrayadas refieren a la posición del TSS anotado

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAgAAAAJBAMAAAD9fXAdAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVKvN74lEInYy3WYQmbv8EmWgAAAAEUlEQVR4XmP8z8DAwMRAHAEALZYBEce9Dw0AAAAASUVORK5CYII=>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAALBAMAAABBvoqbAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAECJEVHaJq7vd75nNZjIrqulnAAAAEklEQVR4XmP8z8DwkYkBCEgiAGhhAgZfdY9rAAAAAElFTkSuQmCC>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAANBAMAAACX52mGAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVLvvmRCrIt0yRGbNdol/YymBAAAAEklEQVR4XmP8z8DwkYkBCMgkAHy7AgpX6bqBAAAAAElFTkSuQmCC>

[image4]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAANBAMAAACX52mGAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVLvvmRCrIt0yRGbNdol/YymBAAAAEklEQVR4XmP8z8DwkYkBCMgkAHy7AgpX6bqBAAAAAElFTkSuQmCC>

[image5]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAKEAAAAwBAMAAACLRQdjAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMARO/dzburiWYQVHaZMiKgQ028AAAAOklEQVR4Xu3MoRHAIADAQGD/ZZkAfC1vepeXEZlnWHt9y7OORkejo9HR6Gh0NDoaHY2ORkejo/GH4wVwAAJQRw8HJwAAAABJRU5ErkJggg==>

[image6]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAANBAMAAACX52mGAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVLvvmRCrIt0yRGbNdol/YymBAAAAEklEQVR4XmP8z8DwkYkBCMgkAHy7AgpX6bqBAAAAAElFTkSuQmCC>

[image7]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAALMAAAAwBAMAAACoHla2AAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMARO/dzburiWYQVHaZMiKgQ028AAAAOUlEQVR4Xu3MoQEAIACAMPX/Z71A+7I2FgnMMz7Zy/JOa7RGa7RGa7RGa7RGa7RGa7RGa7RGa3xcX/9kAlCB4WOCAAAAAElFTkSuQmCC>