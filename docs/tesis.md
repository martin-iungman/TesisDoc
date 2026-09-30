# To do

Por hacer:

- Potencia: todo  
- M2C como figura aparte?  
- M6: editar para indicar calles  
-   
- M9C (heatmap cage): modificar por el del review. modificar leyenda, texto y discutir sobre el upstream  
-   
-   
- Nueva figura sobre shape y tissue specificity??  
- Seccion coocurrencia  
- Tabla de features con criterio de dicotomia  
- M12 pasar a español  
- MyM: PUFFIN. Pensar si poner figura, armar texto  
- Revisar texto ROC-AUC. Hacer alguna figura explicativa?

- Los datos de validacion de ruido y replicabiliad que dan feos… los incluyo???  
- R2C: cambiar el theme  
-   
- Leyenda R6  
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

**Experimentales:**  
	Descripción general del enfoque

- Pipeline

Construcción de los plásmidos reporteros conteniendo la library de promotores basales

- Esquema vector \+ microscopía  
- PCR control integración

	Cultivo celular y transfección  
	Citometría de flujo y clasificación de células

- Potencia??  
- Esquema de density con los gates

	Control por spike-in de celulas  
	Co-extraccion de ADN/ARN y secuenciación de amplicones  
	RT-qPCR de EGFP

- qPCR   
- PCR loxp (cel no transfectadas)

	Anexo Primers y secuencias

**Metodología bioinformática:**  
**Procesamiento de datos y cuantificación de la actividad promotora**  
	Procesamiento de los datos de secuenciacion

- Pipeline preprocessing  
- sample effort \+ bimodal 

	**Caracterización de la library**  
	Composicion y diseño de la library

- type \+ esquema lib \+ heatmap cage

Elementos de secuencia del promotor basal

- epd+cgi+tss+gc

	Patrones de conservación de los promotores basales

- young 2015 \+ phylop

Patrones de actividad endógena de los promotores

- tissue sp \+ shape

Coocurrencia de las características de los promotores basales

- coocurrencia

	Algoritmos de predicción de la actividad promotora basados en la secuencia

- esquema puffin

	Análisis de datos masivos de ChIP-seq  
**Evaluación del efecto de las características del promotor sobre los patrones transcripcionales medidos**  
	Asociacion de las características del promotor con la actividad transcripcional  
	Asociacion de las características del promotor con el ruido transcripcional  
Análisis de enriquecimiento funcional de la unión de TFs  
**Promotores alternativos**  
	Clasificación de promotores alternativos

**Resultados**  
Estimación masiva de las propiedades transcripcionales de promotores basales humanos

- representatividad por replica venn  
- sesgo representatividad: venn presort 2022 \+ gc bias (bars y violin)  
- spike-in  
- densidad individual \+ densidad ejemplo \+ histograma ejemplo \+ correlacion medias validacion  
- correlacion medias (pre y post filtro) \+ venn final replicas \+ var corr replicas

Efectos de la secuencia promotor sobre la fuerza transcripcional

- CpG \+ TATA  
- summary seq  
- summary endo \+ tissue sp?  
- SP1/2 \+ NFYA (subdividido en 2 figuras)  
- remap activity \+ gsea rna pol \+ remap act low  
- integrador?

Efectos de la secuencia promotor sobre el ruido transcripcional

- scatter ruido \+ auc: cgi \+ tata (separado por expresion)  
- summary ruido  
- remap ruido \+ mll1 gsea  
- chipatlas \+ h3k4me3 \+h3k4me3 noCGI

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

# Esqueleto intro

**0- Orden y azar en la materia viva**   
2da ley de la termodinamica: entropia. Universo tiende al desorden  
Schordinger 1944: la vida parece resistir esa tendencia, y su respuesta es que lo hace generando orden de dos maneras distintas — la estadística (orden-desde-desorden) y la "nueva", basada en un código estable y heredable (orden-desde-orden), que él predijo sin conocer aún la estructura del ADN.   
Ese código estable no actúa por sí solo — se ejecuta mediante química de colisiones al azar, difusión, encuentros probabilísticos entre moléculas. El programa es determinista; su ejecución, no.   
Traer a Monod para reformular esto en términos evolutivos: el azar molecular no es ruido a eliminar sino materia prima sobre la que actúa la necesidad   
Traer a Waddington para el otro extremo: cuando el azar amenaza la fidelidad del desarrollo, la evolución también puede construir robustez/canalización para amortiguarlo.   
Wagner — esa misma robustez, lejos de ser solo protección, es lo que permite explorar y acumular variación silenciosa, y eventualmente innovar.   
La arquitectura del promotor basal y su ruido transcripcional intrínseco son, en este marco, un caso concreto y medible de esa tensión general: orden codificado, ejecución estocástica, robustez que no es pasiva sino que habilita innovación regulatoria (nuevos promotores, nuevos usos tisulares). 

**1- Del genotipo al fenotipo regulatorio**   
Objetivo central de la genómica: cómo la secuencia de ADN determina la expresión génica y el fenotipo  
La transcripción como paso regulatorio clave en eucariotas  
Cuánto del comportamiento regulatorio de un gen está "escrito" en su secuencia, independientemente del contexto? (introducción a la pregunta, no abordaje completo)  
Adelantar que esta tesis aborda dos fenotipos regulatorios distintos (actividad y ruido) y que ambos son parte de esta misma pregunta de "lectura" de la secuencia.

**2- Más allá del enhancer: el rol regulatorio del promotor basal**   
Core promoter: región que recluta RNAPII para iniciar la transcripción  
Definición operacional y elementos de secuencia (TATA-box, Inr, DPE, BRE), PIC y factores generales  
El paradigma histórico enhancer-céntrico: los enhancers como responsables de controlar el output transcripcional  
Evidencia reciente que lo cuestiona: la secuencia del promotor sola puede ser predictiva de actividad  
Autonomía transcripcional como propiedad general de los promotores humanos, medida de forma no sesgada a nivel genómico   
El promotor también selecciona qué enhancers vecinos pueden regularlo — no es solo receptor pasivo de señales 

**3- Diversidad arquitectónica de promotores: robustez versus especificidad**   
Promotores focused vs. dispersed/broad, y su implicancia evolutiva en robustez a variación genética  
Asociación housekeeping vs TATA-less/GC-rich/CpG islands/TSS disperso, y tejido-específico vs TATA-box/TSS focused   
*R*obustez y estabilidad vs. capacidad de respuesta precisa, como dos soluciones evolutivas a distintas demandas regulatorias  
Islas CpG y su asociación con alta actividad/autonomía (motivos GC-box, efecto sobre posicionamiento nucleosomal)  
Nucleosomas y marcas de histonas (H3K4me3) como contexto de cromatina asociado a estas arquitecturas  
*Plantear:* ¿la cromatina establece el régimen transcripcional o es un reflejo de un estado ya definido por la secuencia?

**4\. Un gen, varios promotores** 

* La mayoría de los genes humanos tienen múltiples promotores alternativos  
* Origen de promotores alternativos: enhancers, elementos retrovirales endógenos, splicing-mediated activation  
* Plantear la pregunta de economía evolutiva: si ya existe un promotor funcional, ¿qué ventaja da mantener promotores "de repuesto"?  
* Necesidad de un abordaje promotor-céntrico en vez de gen-céntrico  
* *I*dea de promotor principal vs. secundario y el fenómeno de tissue-switching, como pistas de que la multiplicidad podría ser plasticidad y no simple redundancia

**5- Ruido transcripcional. El componente estocástico de la transcripción** 

Definición: componente de variabilidad célula a célula específico de cada gen   
Bursting transcripcional como mecanismo subyacente, ligado a switches de estado del promotor   
Explicar el modelo de bursting (on/off, frecuencia vs. tamaño de burst) para dar base mecanística a la pregunta de desacople media-varianza  
Base genética heredable del ruido, sujeta a selección natural  
¿puede un promotor ser muy activo pero poco ruidoso, o alta actividad y alto ruido van necesariamente ligados? (adelantar el caso TATA como ejemplo de acoplamiento)  
Relevancia en robustez biológica y desarrollo (determinación de destino celular) Antecedentes en levadura (nucleosomas, sitios de unión a TF) y en mamíferos vía scRNA-seq  
Limitación central de scRNA-seq: no separa el efecto del promotor de otras señales genéticas/epigenéticas, ni discrimina promotores alternativos

**6- Abordajes experimentales para el estudio del promotor basal** 

Qué es un MPRA y por qué permite desacoplar secuencia de contexto genómico/epigenético (vuelve a la pregunta de la seccion 3\)  
Variantes existentes (episomales, STARR-seq, SuRE, lentiMPRA) y sus limitaciones frente a un sistema de integración en locus único  
Conectar con la meta-pregunta de predictibilidad: ¿se puede predecir el comportamiento regulatorio de un promotor solo mirando su secuencia? (y reconectar con seccion 1\)  
Modelos como Xpresso o PUFFIN intentan predecir comportamiento regulatorio desde una secuencia.   
Y para que nos sirve todo esto?

**7-**   
**hipotesis, objetivos generales y especificos**  
**descripcion breve de metodologia**  
**hilo de secciones de resultados**

# Prologo

**Orden, azar e innovación en la materia viva**

Un cristal de sal y un huracán son, en un sentido físico preciso, dos maneras distintas de producir orden. El cristal alcanza su estructura minimizando energía libre, cerca del equilibrio termodinámico: una vez formado, no necesita nada del exterior para persistir. El huracán, en cambio, es una estructura disipativa, en el sentido que le dio Ilya Prigogine: existe únicamente mientras un flujo continuo de energía lo atraviesa, y se disuelve en el aire quieto apenas ese flujo se interrumpe. La vida, con toda su complejidad, pertenece inequívocamente a la segunda familia. Ningún organismo alcanza jamás el equilibrio con su entorno mientras vive, y sostiene su organización interna a fuerza de un metabolismo que nunca se detiene. Al final de cuentas, para casi cualquier célula, el equilibrio termodinámico equivale a la muerte. Esto no contradice la ley más inapelable de la física, según la cual la entropía de un sistema cerrado nunca disminuye: la vida no es una excepción a esa ley, sino una de sus consecuencias más elaboradas, orden que se paga constantemente con entropía exportada al entorno.[^1]

Erwin Schrödinger se hizo una pregunta contigua a esta en *¿Qué es la vida?* (1944), un pequeño libro que terminaría influyendo decisivamente en el nacimiento de la biología molecular: si el orden del metabolismo es disipativo, ¿de qué tipo es el orden, mucho más preciso y heredable, de la información genética? Para describirlo, eligió una imagen tomada justamente de la otra familia, la del equilibrio: la de un cristal. Pero le agregó un adjetivo que la vuelve extraña: aperiódico. Un cristal ordinario es orden repetitivo y, por eso mismo, pobre en información, ya que toda su estructura se deduce de una celda mínima que se repite. Un cristal *aperiódico*, en cambio, conservaría la estabilidad física de un cristal (de equilibrio, no disipativa; no necesita energía continua para no degradarse) sin sacrificar la posibilidad de portar información variable, casi ilimitada, en su secuencia. Esa doble condición, estable como un cristal pero informativo como un mensaje, es lo que Schrödinger proponía como sustrato de la herencia, casi una década antes de que Watson y Crick describieran la estructura del ADN.

Pero identificar el ADN como este "cristal aperiódico" no resuelve la paradoja, apenas la desplaza a una escala menor. Porque ese código estable, de equilibrio, guardado sin costo energético en cada célula, no actúa en el vacío. Leerlo, ejecutarlo, es un proceso disipativo como cualquier otro metabolismo: requiere energía continua, y ocurre a través de una química que, a escala molecular, es fundamentalmente estocástica. Las moléculas se encuentran por difusión y colisión al azar; los factores de transcripción exploran el núcleo celular en una caminata browniana hasta encontrar, o no, su sitio de unión; la polimerasa inicia la transcripción en eventos discretos e impredecibles a nivel individual. El programa genético es determinista y estable; su lectura, molécula por molécula, célula por célula, es disipativa y azarosa.

Jacques Monod, en *El azar y la necesidad* (1970), no vio en esta convivencia un defecto a corregir sino una condición productiva. Las variaciones que alimentan la evolución (mutaciones, pero también, podríamos agregar hoy, la variabilidad estocástica en la expresión de un gen) surgen del puro azar molecular, y es únicamente la selección natural, la "necesidad", la que opera después sobre esa variación, reteniendo lo que funciona. Monod resumió esta idea en una imagen que se volvió célebre: el azar, escribió, queda *"pris sur l'aile"*, atrapado al vuelo, preservado y reproducido por la maquinaria de la herencia, y así convertido en orden, regla, necesidad.

Pero si ese azar molecular fuera completamente libre, el desarrollo de un organismo, que depende de que miles de decisiones celulares ocurran de manera coordinada y reproducible, sería inviable. Conrad Waddington capturó el límite de esa libertad con la metáfora del paisaje epigenético: el desarrollo transcurre como una bola rodando por un paisaje de valles y colinas, y aunque pequeñas perturbaciones puedan desviar levemente su trayectoria, los valles (las trayectorias de desarrollo canalizadas) tienden a devolverla al mismo destino final. A esta capacidad de amortiguar la variación, ya sea genética o ambiental, sin alterar el fenotipo resultante, Waddington la llamó canalización, y hoy la entendemos como una forma particular de robustez biológica: el freno necesario para que el azar de Monod no desborde el desarrollo.

Podría pensarse que ese freno es puramente conservador, un mecanismo que protege al organismo del cambio y que por lo tanto se opone a la innovación evolutiva. Andreas Wagner (*The Origins of Evolutionary Innovations*, 2011\) mostró que ocurre casi exactamente lo contrario. En sistemas biológicos tan diversos como el plegado de moléculas de RNA o las redes metabólicas, la robustez frente a la variación, genética o estocástica, no impide la exploración de nuevas soluciones fenotípicas: la habilita. Un sistema robusto puede acumular, de manera silenciosa y sin consecuencias fenotípicas inmediatas, una cantidad considerable de variación "críptica" en sus componentes subyacentes, disponible para producir innovación el día en que el contexto cambia. La canalización de Waddington, entonces, no clausura la posibilidad de cambio que abrió Monod, la acumula en reserva. Conviene aclarar, sin embargo, que no toda la robustez observada en los genomas necesita leerse en estos términos adaptativos: el genetista de poblaciones Michael Lynch ha argumentado que buena parte de la arquitectura genómica de los eucariotas es, más bien, un subproducto de procesos neutrales, como la deriva genética en poblaciones de tamaño efectivo reducido, antes que el resultado de una selección positiva por sus beneficios de evolvabilidad (Lynch, 2007). La robustez, en otras palabras, no siempre se elige: a veces, simplemente, se tolera.

Esa reserva de variación silenciosa, sin embargo, rara vez se traduce en piezas genuinamente nuevas cuando la innovación finalmente ocurre. François Jacob lo señaló con una imagen tan célebre como la de Waddington: la evolución no procede como un ingeniero que diseña una solución óptima desde cero a partir de un plano, sino como un *bricoleur*, el que arregla con lo que tiene a mano, el que ata todo con alambres, y que resuelve cada problema nuevo reutilizando y recombinando piezas con las que ya contaba, aunque hubieran sido "diseñadas" originalmente para otra cosa. Un enhancer que, en un gen ya establecido, comienza a funcionar como promotor alternativo, o un elemento retroviral endógeno que aporta, sin haber sido nunca seleccionado para eso, un sitio de inicio de transcripción funcional, son, en el sentido más literal, actos de bricolage molecular. La innovación regulatoria, vista así, rara vez inventa una pieza nueva; casi siempre reutiliza una vieja con un propósito distinto.

Antes de bajar a ese promotor concreto, falta todavía un último ingrediente: un lenguaje capaz de medir cuánto de lo discutido hasta acá, orden codificado, indeterminación, está efectivamente contenido en una secuencia particular. Curiosamente, la palabra "entropía" reaparece un siglo después de Boltzmann en un contexto completamente distinto para ofrecer exactamente eso. Cuando Claude Shannon formalizó en 1948 una teoría matemática de la comunicación, necesitaba una medida de cuánta incertidumbre contiene una fuente de señales, o, en sentido inverso, cuánta información transmite un mensaje al resolver esa incertidumbre. La fórmula a la que llegó es, salvo por una constante, formalmente idéntica a la entropía de Boltzmann, y no por casualidad adoptó el mismo nombre. Bajo esta luz, una secuencia de ADN puede leerse también como una fuente de información: un arreglo no aleatorio de bases cuyo grado de estructura, cuán lejos está de la equiprobabilidad, es en principio cuantificable, y del cual depende cuánto puede "saberse" de antemano sobre el proceso que esa secuencia dirige.

Todas estas ideas convergen en un mismo escenario concreto y medible: el promotor basal de un gen. Un promotor es, en esencia, un fragmento acotado de secuencia (estable, de equilibrio, en el sentido en que Schrödinger hablaba de un cristal aperiódico) que codifica un programa: cuándo, cuánto y con qué precisión debe iniciarse la transcripción de un gen. En el sentido en que Shannon hablaba de una fuente informativa, cuanto más estructurada y menos arbitraria es esa secuencia, más predecible debería ser, en principio, el proceso que gobierna. Pero ese programa se ejecuta como todo proceso disipativo, mediante encuentros moleculares estocásticos, célula por célula, en el sentido en que Monod hablaba de azar productivo, y solo permanece dentro de límites viables gracias a mecanismos de robustez análogos a la canalización que describió Waddington. La *actividad* transcripcional que un promotor dirige mide cuánto se ejecuta ese programa: cuánto ARN produce, en promedio, cada vez que se pone en marcha. El *ruido* transcripcional que ese mismo promotor tolera o amplifica mide, en cambio, cuán azaroso es cada evento individual de esa ejecución, cuán impredecible es el instante exacto en que la transcripción se enciende o se apaga, algo que, al cambiar de escala, se vuelve visible en la variabilidad generada entre células. Y la arquitectura particular de cada promotor (focalizada o dispersa, robusta o precisa, única o multiplicada en variantes alternativas dentro de un mismo gen, muchas veces reclutadas por bricolage evolutivo a partir de enhancers o elementos móviles preexistentes, como describió Jacob) es, posiblemente, el resultado visible de una historia evolutiva que, como sugiere Wagner, no solo ha tenido que contener esta tensión entre orden y azar, sino que se ha valido de ella para innovar.

Esta tesis parte de esa tensión general para hacerla, en las páginas que siguen, una pregunta molecular concreta y experimentalmente abordable: ¿qué elementos de la secuencia de un promotor basal humano determinan cuánto se ejecuta ese programa, y cuán al azar ocurre cada uno de esos eventos de ejecución, y qué nos dice esto sobre la lógica evolutiva de la regulación génica? ¿Cómo tensiona la evolución, en cada promotor, la información con el caos, la robustez con la posibilidad de cambiar?

# Intro

**1- Del genotipo al fenotipo regulatorio** 

Explicar cómo la información contenida en el genotipo se traduce en las características observables de un organismo —su fenotipo— es una de las preguntas fundacionales de la genética, presente desde que Wilhelm Johannsen acuñó ambos términos a comienzos del siglo XX para distinguir la constitución hereditaria de un organismo de su manifestación observable (Johannsen, 1911\). Durante buena parte del siglo XX, esta pregunta se abordó casi exclusivamente a través de la secuencia codificante: qué proteína produce un gen, y cómo las variantes en esa secuencia alteran su función. Sin embargo, el advenimiento de la genómica funcional y de los estudios de asociación de genoma completo (GWAS) reveló un panorama considerablemente más complejo. La gran mayoría de las variantes genéticas asociadas a enfermedades y rasgos complejos en humanos no se localiza en regiones codificantes, sino que se concentra en el ADN regulatorio no codificante (Maurano et al., 2012\). Este hallazgo desplazó buena parte del interés de la genética funcional desde la pregunta de qué hace una proteína hacia una pregunta distinta, y en muchos sentidos más difícil: qué determina cuánto, cuándo y en qué célula se expresa un gen.

Un objetivo central de la genómica moderna es, entonces, comprender cómo la secuencia de ADN determina la expresión génica y, a través de ella, el fenotipo. En la enorme mayoría de los organismos eucariotas, la transcripción constituye el paso regulatorio primario de este proceso: es en la decisión de transcribir o no un gen, y en qué magnitud hacerlo, donde se define buena parte del destino de la información genética (Levine y Tjian, 2003; Ptashne y Gann, 1997). Esta decisión no ocurre de manera espontánea. Requiere el reclutamiento y ensamblado de la ARN Polimerasa II (RNAPII) sobre el ADN, en un proceso que depende de una región de secuencia acotada alrededor del sitio de inicio de la transcripción (*transcription start site*, TSS): el promotor basal o *core promoter* (Sainsbury et al., 2015).

Durante buena parte del desarrollo de la biología molecular, los promotores basales fueron entendidos principalmente como plataformas de ensamblaje relativamente pasivas, cuya actividad dependía en gran medida de señales provistas por elementos regulatorios distales. Este paradigma tiene su origen en el descubrimiento de los primeros enhancers a comienzos de la década de 1980 (Banerji et al., 1981), y se consolidó en las décadas siguientes a partir de estudios genómicos y del desarrollo en distintos organismos modelo, que establecieron a los enhancers como los principales elementos organizadores de los programas de expresión espaciotemporal durante el desarrollo animal (Levine, 2010; Spitz y Furlong, 2012; Long et al., 2016). En organismos metazoos en particular, se asumió durante décadas que eran los enhancers, y no los promotores, los principales responsables de determinar cuánto, cuándo y en qué tipo celular se transcribe un gen (Bulger y Groudine, 2011). 

Trabajos más recientes han puesto en cuestión esta jerarquía. Distintos estudios, basados en el diseño de bibliotecas masivas de secuencias promotoras sintéticas o aleatorias, mostraron que la secuencia del promotor basal, evaluada de manera aislada, puede ser sorprendentemente predictiva de su actividad transcripcional (de Boer et al., 2020; Kwasnieski et al., 2012; Agarwal y Shendure, 2020). En la misma línea, un mapeo genómico no sesgado de la actividad autónoma de promotores humanos —es decir, la capacidad de un fragmento de ADN de iniciar transcripción en ausencia de otros elementos regulatorios— encontró que esta autonomía es, de hecho, una propiedad común a la mayoría de los promotores humanos, y no una excepción restringida a un subconjunto particular (van Arensbergen et al., 2017). A esto se suma evidencia de que los propios promotores basales participan activamente en la selección de qué enhancers vecinos podrán regularlos, en lugar de responder de manera indiscriminada a cualquier señal distal disponible (Zabidi et al., 2015; Haberle y Stark, 2018). En conjunto, estos hallazgos sugieren que una porción sustancial del comportamiento regulatorio de un gen podría estar, en efecto, escrita directamente en la secuencia de su promotor basal.

Esta reevaluación del rol del promotor no resuelve, sin embargo, una pregunta más amplia y todavía abierta: si tanta información regulatoria está efectivamente codificada en esta región relativamente corta de secuencia, ¿cuánto del comportamiento regulatorio de un gen está, en efecto, "escrito" en ella, independientemente del contexto genómico en el que se encuentre? Esta tesis aborda esa pregunta a través de dos fenotipos regulatorios distintos, aunque estrechamente relacionados entre sí: la *actividad* transcripcional que un promotor es capaz de dirigir, y el *ruido*, o variabilidad célula a célula, que esa misma actividad conlleva. Ambos pueden entenderse como dos maneras diferentes pero complementarias de "leer" cuánta información regulatoria contiene la secuencia de un promotor basal. Responder esta pregunta requiere, en primer lugar, precisar qué constituye un promotor basal y qué elementos de secuencia lo definen. 

**2- Más allá del enhancer: el rol regulatorio del promotor basal**

El promotor basal, o *core promoter*, es la región de ADN que rodea el sitio de inicio de la transcripción y que recluta directamente a la ARN Polimerasa II para dar comienzo a este proceso. Conviene aclarar, antes de avanzar, que cualquier definición de esta región es necesariamente operacional y no una categoría natural con límites fijos. Distintos autores la delimitan de formas ligeramente distintas según el criterio metodológico empleado (footprint de unión de TFIID, mapeo funcional por deleción, ventana de distancia al TSS), y no existe un elemento de secuencia universal que permita trazar un borde inequívoco (**Smale y Kadonaga, 2003**; **Kadonaga, 2012**). En la práctica, suele considerarse el tramo comprendido entre aproximadamente 35 y 40 pares de bases hacia arriba y hacia abajo del TSS, pero esa cifra es una convención útil, no una frontera biológica estricta.

Dentro de esa región ocurre el ensamblado del complejo de preiniciación (*preinitiation complex*, PIC). El modelo mejor descripto, obtenido en gran medida a partir de reconstituciones in vitro, propone un orden aproximadamente jerárquico: el factor general TFIID, a través de su subunidad TBP y de un conjunto de proteínas asociadas (TAFs), reconoce primero los elementos de secuencia del promotor; TFIIA estabiliza esta unión; TFIIB actúa como puente hacia el complejo formado por la RNAPII y TFIIF; y finalmente se suman TFIIE y TFIIH, completando el PIC. TFIIH cumple un rol particular, ya que su actividad helicasa desenrolla localmente el ADN alrededor del TSS, y su actividad quinasa fosforila el dominio carboxilo-terminal de la RNAPII, un paso necesario para que la polimerasa escape del promotor e inicie la elongación (Sainsbury et al., 2015). Vale aclarar que este esquema describe un caso relativamente canónico, y no un mecanismo único ni fijo: existen complejos TFIID de composición alternativa, factores específicos de ciertos promotores o tipos celulares que pueden reemplazar a componentes considerados generales (como TBP), y evidencia de que el orden real de ensamblado en la célula no siempre reproduce con exactitud la secuencia observada en sistemas reconstituidos (**Cramer, 2019**). La maquinaria basal de transcripción, en otras palabras, admite variantes.

Buena parte de este ensamblado depende, además, de qué elementos de secuencia estén presentes en el promotor, y ningún elemento por sí solo alcanza para explicarlo. La caja TATA, reconocida directamente por TBP, se ubica típicamente entre 25 y 30 pares de bases antes del TSS y es uno de los elementos mejor caracterizados, aunque está presente en apenas una fracción minoritaria de los promotores humanos (Smale y Kadonaga, 2003). El iniciador (*Inr*) abarca el propio sitio de inicio y admite dos niveles de definición: una versión mínima, que solo requiere una pirimidina en la posición −1 y una purina en la posición \+1 (Py-Pu), y una versión más estricta, reconstruida a partir de datos genómicos de alta resolución, que identifica un consenso más completo (BBCABW) en promotores con un único sitio de inicio bien definido (**Vo ngoc et al., 2017**). Ese mismo trabajo mostró que los promotores con un Inr consenso tienden a carecer de caja TATA, y que su presencia no está particularmente asociada a la de islas CpG. A estos se suman el elemento de reconocimiento de TFIIB (*BRE*), adyacente a un subconjunto de cajas TATA, y el elemento promotor *downstream* (*DPE*), que requiere de un Inr para funcionar y es especialmente frecuente en promotores sin caja TATA (Smale y Kadonaga, 2003; Kadonaga, 2012).

Otros dos elementos frecuentes son la caja CCAAT, reconocida por el factor trimérico NF-Y, y la caja GC, reconocida principalmente por proteínas de la familia Sp (Sp1, Sp2, Sp3). Ambos suelen coexistir en un mismo promotor con mayor frecuencia de la esperada por azar, y su presencia está fuertemente asociada a la ausencia de caja TATA. La caja GC, en particular, suele encontrarse en contextos ricos en dinucleótidos CG, es decir, en o cerca de islas CpG.

Como se planteó antes, el promotor basal fue entendido durante décadas como una plataforma más bien pasiva, cuya activación dependía casi por completo de señales aportadas por enhancers distales. Bajo este modelo, conocido como el *enhanceosoma*, la especificidad de un programa transcripcional residía en la combinación particular de factores de transcripción que se ensamblaban sobre un enhancer, mientras que el promotor cumplía una función mayormente indiferenciada de recepción de esa señal (**Spitz y Furlong, 2012**). La evidencia reunida en los últimos años, de que la secuencia del promotor por sí sola puede predecir buena parte de su actividad transcripcional, obliga a matizar esta idea.

Quizás el hallazgo más contundente en esta dirección es que los propios promotores basales no responden de manera indiferenciada a cualquier enhancer disponible, sino que existe una compatibilidad selectiva entre ciertos tipos de promotores y ciertos tipos de enhancers. Un estudio en *Drosophila* mostró que los promotores asociados a genes *housekeeping* y los promotores asociados a genes de desarrollo responden preferentemente a enhancers de su misma categoría, y que esta preferencia depende de elementos de secuencia específicos dentro del propio *core promoter*, y no solo de la identidad del enhancer (Zabidi et al., 2015). En la misma línea, un análisis más amplio de secuencias de promotores eucariotas propuso que la diversidad de elementos del *core promoter* codifica, en la práctica, distintas gramáticas regulatorias, cada una compatible con un conjunto particular de mecanismos de activación (Haberle y Stark, 2018). El promotor, entonces, no solo puede generar actividad transcripcional de manera autónoma: también parece participar activamente en decidir qué señales regulatorias externas está dispuesto a aceptar.

Esta evidencia reposiciona al promotor basal como un actor con información propia, y no como un simple ejecutor de decisiones tomadas en otro lugar del genoma. Pero también deja planteada una pregunta que todavía no fue respondida: si la arquitectura del promotor determina tanto su nivel de actividad como su compatibilidad con distintos enhancers, ¿determina también su comportamiento estocástico, es decir, cuán variable es esa actividad de una célula a otra? Para abordar esta pregunta hace falta, primero, entender por qué no todos los promotores tienen la misma arquitectura, y qué lógica evolutiva podría explicar esa diversidad.

**3- Diversidad arquitectónica de promotores: robustez y especificidad**

No todos los promotores inician la transcripción de la misma manera. Un mapeo genómico de sitios de inicio de transcripción en mamíferos, realizado mediante la técnica de CAGE (*cap analysis of gene expression*), mostró que los promotores pueden agruparse en dos grandes clases según su arquitectura: promotores angostos o focalizados, donde la transcripción comienza en una posición única o en un grupo estrecho de posiciones vecinas, y promotores anchos o dispersos, donde el inicio se distribuye entre múltiples posiciones alternativas a lo largo de una ventana más amplia (Carninci et al., 2006). Esta distinción, conocida como la "forma" (*shape*) del promotor, suele cuantificarse con una medida de entropía análoga a la que Shannon propuso para cuantificar la incertidumbre de una fuente de información: cuanto más repartida está la probabilidad de inicio entre distintas posiciones, mayor es la entropía, y más disperso es el promotor.

La forma de un promotor no es independiente del tipo de gen al que pertenece. Los promotores focalizados están fuertemente asociados a la presencia de una caja TATA y suelen encontrarse en genes de expresión regulada, específica de tejido o de condición. Los promotores dispersos, en cambio, están asociados a la presencia de islas CpG, carecen mayormente de caja TATA, y son característicos de genes housekeeping, expresados de manera más constante en distintos tipos celulares (Carninci et al., 2006). Esta asociación no parece casual: se ha propuesto que los promotores dispersos son más robustos frente a variación genética que los focalizados, ya que contar con múltiples posiciones de inicio alternativas amortigua el efecto de una mutación puntual sobre cualquiera de ellas, mientras que un promotor focalizados, dependiente de un único sitio bien definido, es más vulnerable a que esa misma mutación altere su comportamiento (Carninci et al., 2006).

Esta asociación sugiere una lógica de diseño evolutivo relativamente intuitiva. Un gen housekeeping necesita mantener un nivel de expresión estable en prácticamente cualquier contexto celular, y para eso puede beneficiarse de una arquitectura robusta, tolerante a variación, aunque menos precisa. Un gen regulado de forma específica, en cambio, necesita poder activarse (o no) con precisión en el momento y el tejido correctos, algo que se logra mejor con un sitio de inicio único y bien definido, aun a costa de una mayor fragilidad frente a mutaciones. Robustez y especificidad, en este sentido, parecen operar como dos soluciones distintas a dos demandas regulatorias distintas, más que como una jerarquía de arquitecturas "buenas" y "malas".

Esta arquitectura no ocurre en el vacío, sino embebida en un contexto de cromatina particular. Las islas CpG (regiones cortas, ricas en GC y en dinucleótidos CG, que escapan a la metilación característica del resto del genoma) tienden a desestabilizar la ocupación de nucleosomas y a favorecer un estado de cromatina permisivo para la transcripción, lo cual explica en parte por qué los promotores que las contienen suelen ser transcripcionalmente más activos y autónomos (Deaton y Bird, 2011\). Este mismo contexto de cromatina, evaluado de manera más general, distingue a los genes según su capacidad de responder a cambios de condición: en levadura, los promotores con un nucleosoma bien posicionado inmediatamente antes del sitio de inicio muestran mayor plasticidad transcripcional frente a distintos estímulos, mientras que aquellos con una región libre de nucleosomas amplia y estable en esa misma posición tienden a mostrar una expresión más constitutiva. Ambas estrategias, además, se asocian a niveles distintos de ruido transcripcional (Tirosh y Barkai, 2008).

Entre las marcas de cromatina asociadas a este contexto, la trimetilación de la lisina 4 de la histona H3 (H3K4me3) ocupa un lugar particular. Se trata de una de las marcas más consistentemente enriquecidas en los sitios de inicio de transcripción activos, presente en la gran mayoría de los promotores humanos con actividad detectable, independientemente del nivel de expresión del gen (Guenther et al., 2007). Su relación con la maquinaria basal no es solo correlativa: TAF3, una de las subunidades de TFIID, reconoce directamente a H3K4me3 y contribuye a anclar el complejo de preiniciación sobre el promotor, lo que sugiere que esta marca no solo acompaña a los promotores activos, sino que participa activamente en su reconocimiento por parte de la maquinaria de transcripción (Vermeulen et al., 2007). H3K4me3 tiende además a coincidir con la presencia de islas CpG, reforzando el mismo entorno de cromatina permisivo ya descripto.

Esta última pregunta no es un detalle menor. Un estudio en *Drosophila* mostró que la forma de un promotor puede evolucionar de manera independiente de su nivel de actividad, y que las variantes genéticas naturales que modifican esa forma suelen aumentar el ruido transcripcional, con evidencia de que la minimización del ruido podría ser, en sí misma, una presión selectiva relevante en la evolución de los promotores (Schor et al., 2017). Es decir, la arquitectura de un promotor no solo determina cuánto se transcribe un gen, sino también cuán variable es esa transcripción de una célula a otra, lo cual empieza a insinuar que actividad y ruido, lejos de ser fenotipos independientes, podrían estar ambos escritos en la misma arquitectura de secuencia.

Queda, sin embargo, una pregunta sin resolver: si la arquitectura de secuencia predice tanto la actividad como (al menos en parte) el ruido de un promotor, ¿lo hace de manera directa, o actúa a través del contexto de cromatina que esa misma secuencia genera a su alrededor? Distinguir estas dos posibilidades requiere herramientas capaces de aislar la secuencia del promotor de su contexto genómico habitual, algo que retomamos más adelante. Antes de eso, sin embargo, conviene atender a otra fuente de complejidad que la descripción hecha hasta acá pasó por alto: la mayoría de los genes no tiene un único promotor, sino varios.

**4- Un gen, varios promotores**

Hasta acá describimos al promotor basal como si cada gen tuviera uno solo. En la práctica, esa es la excepción y no la regla: la mayoría de los genes humanos poseen múltiples promotores alternativos, capaces de iniciar la transcripción de distintas isoformas de un mismo gen desde posiciones genómicas diferentes (Carninci et al., 2006). Esta multiplicidad no es un fenómeno marginal ni exclusivo de un puñado de genes particulares, sino una característica extendida del genoma humano, lo cual complica cualquier intento de describir la regulación transcripcional gen por gen, como si cada uno tuviera un único punto de entrada a la transcripción.

¿De dónde surgen estos promotores adicionales? La evidencia disponible señala más de un origen posible. En algunos casos, elementos que originalmente cumplían función de enhancer en un gen ya establecido adquieren, a lo largo de la evolución, la capacidad de iniciar transcripción por sí mismos, funcionando de hecho como un promotor alternativo intragénico (Kowalczyk et al., 2012; Carelli et al., 2016, 2018). En otros casos, son elementos retrovirales endógenos, remanentes de infecciones virales ancestrales incorporadas al genoma, los que aportan sitios de inicio de transcripción funcionales, muchas veces sin haber sido nunca seleccionados originalmente para cumplir ese rol (Li et al., 2018; Young et al., 2015). Existe además un mecanismo distinto, ligado no a la aparición de un nuevo elemento de secuencia sino a la propia arquitectura del gen: la activación mediada por splicing, en la que el propio proceso de corte y empalme del pre-ARNm favorece el uso de un sitio de inicio alternativo (Fiszbein et al., 2019). Los promotores alternativos, en otras palabras, no tienen un único origen evolutivo, sino que emergen repetidamente por vías distintas.

Esta multiplicidad plantea una pregunta que no es solo mecanística, sino también de economía evolutiva: si un gen ya cuenta con un promotor funcional, capaz de sostener su expresión, ¿qué ventaja hay en mantener uno o más promotores adicionales, aparentemente redundantes? Una respuesta posible es que no lo sean: que en lugar de una redundancia sin costo, la multiplicidad de promotores represente una forma de plasticidad regulatoria, que le permite a un mismo gen sostener un programa de expresión estable en su contexto principal mientras conserva la capacidad de responder de manera distinta en contextos secundarios.

Esta pregunta tiene, además, una consecuencia práctica ineludible para cualquier estudio de regulación transcripcional: si un gen puede expresarse desde más de un promotor, y esos promotores pueden tener propiedades regulatorias distintas entre sí, entonces medir "la" actividad o "el" ruido de un gen, como si fuera un valor único, oculta información relevante. Se vuelve necesario, en cambio, adoptar una perspectiva centrada en el promotor individual, y no en el gen como unidad indivisible.

Algo de esta lógica de plasticidad ya empieza a insinuarse en la relación entre promotores dentro de un mismo gen. Cuando un gen tiene varios promotores, suele haber uno que domina la expresión en la mayoría de los tejidos, el promotor principal, y otro u otros que contribuyen de manera minoritaria, los promotores secundarios. Pero esta jerarquía no es fija: en ciertos tejidos, un promotor normalmente secundario puede convertirse en el promotor dominante, desplazando al principal, en lo que suele describirse como un patrón de *tissue-switching*. La existencia de este fenómeno sugiere que, al menos para un subconjunto de genes, la multiplicidad de promotores no es simplemente un excedente regulatorio tolerado por la evolución, sino un mecanismo activo que permite a un gen adaptar su programa de expresión según el contexto celular, sin necesidad de rediseñar su regulación desde cero.

**5- Ruido transcripcional. El componente estocástico de la transcripción**

Incluso dentro de una población de células genéticamente idénticas, del mismo tipo y en el mismo ambiente, el nivel de expresión de un gen dado puede variar considerablemente de una célula a otra. Esta variabilidad no es un artefacto técnico ni un error de medición: es, en gran medida, una propiedad genuina del proceso biológico. Elowitz y colaboradores, en uno de los trabajos fundacionales del campo, distinguieron dos componentes de esta variabilidad usando un sistema de doble reportero fluorescente en *Escherichia coli*: el **ruido extrínseco**, que afecta simultáneamente a muchos o todos los genes de una célula (por diferencias en el microambiente, el estado del ciclo celular o la capacidad general de expresión de esa célula en particular), y el **ruido intrínseco**, específico de cada gen, que depende de sus propios mecanismos moleculares de expresión y ocurre incluso cuando esos factores más generales se mantienen constantes (Elowitz, 2002). Esta tesis se ocupa exclusivamente del segundo: el ruido transcripcional intrínseco, atribuible a la identidad y el comportamiento de un promotor particular.

El origen mecanístico de este ruido intrínseco está bien establecido: la transcripción no ocurre de manera continua, sino en pulsos discretos, conocidos como ráfagas o *bursts*, separados por períodos de inactividad. Este comportamiento se ha observado de manera consistente en organismos tan distintos como bacterias, levaduras y células de mamífero (Raser y O'Shea, 2005). El modelo más simple para describirlo asume que el promotor alterna entre un estado activo, capaz de iniciar transcripción, y uno inactivo, y que esta dinámica puede caracterizarse mediante dos parámetros: el tamaño medio de cada ráfaga (cuánto ARN se produce durante el período activo) y su frecuencia (cuán seguido ocurre ese período). Ambos parámetros afectan de manera distinta al fenotipo observable: aumentar el tamaño de las ráfagas incrementa el nivel medio de expresión sin modificar demasiado el ruido, mientras que aumentar su frecuencia incrementa el nivel medio y, al mismo tiempo, reduce el ruido, ya que un mayor número de eventos independientes por unidad de tiempo tiende a promediarse con mayor precisión (Dar et al., 2012; Hornung et al., 2012). La actividad y el ruido de un promotor, entonces, no son magnitudes independientes entre sí: comparten un origen mecanístico común en la dinámica de estas ráfagas, aunque no varíen necesariamente juntas ni en la misma dirección.

Buena parte de lo que se sabe sobre qué características de un promotor determinan esta dinámica proviene de estudios sistemáticos en levadura. Se ha mostrado que la ocupación de nucleosomas en la región del promotor, así como el número y la afinidad de los sitios de unión a factores de transcripción, afectan tanto la fuerza transcripcional como el ruido, aunque a través de mecanismos distintos y con efectos que no siempre van en la misma dirección (Tirosh y Barkai, 2008; Dadiani et al., 2013). Estudios con bibliotecas de miles de promotores sintéticos, diseñados específicamente para variar de manera controlada estas características, confirmaron que modificaciones incluso mínimas en la secuencia de un promotor (a veces un cambio de una sola base) pueden alterar sustancialmente su nivel de ruido, independientemente de su efecto sobre la actividad media (Sharon et al., 2014).

Esta sensibilidad del ruido a variantes de secuencia relativamente sutiles tiene una consecuencia evolutiva importante: si el ruido depende de la secuencia del promotor, entonces también puede estar sujeto a selección natural. La evidencia en este sentido proviene de organismos distintos. En levadura, un análisis del promotor del gen *TDH3* mostró que los polimorfismos naturales asociados a un aumento del ruido transcripcional eran considerablemente menos frecuentes en la población de lo esperado si las mutaciones ocurrieran al azar, lo cual sugiere la acción de selección purificadora en contra de variantes ruidosas (Metzger et al., 2015). En *Drosophila*, la comparación entre los haplotipos naturales más frecuentes de un promotor y una variante recombinante artificial mostró que esta última generaba, en la mayoría de los genes evaluados, más ruido que cualquiera de los haplotipos naturales, sugiriendo que las combinaciones de variantes favorecidas en la población son precisamente las que minimizan el ruido (Schor et al., 2017). El ruido transcripcional, lejos de ser un subproducto tolerado pasivamente, parece ser en muchos casos un blanco activo de la selección.

Esto no implica, sin embargo, que el ruido bajo sea siempre la opción evolutivamente favorecida. Distintos genes parecen requerir distintos regímenes de variabilidad: genes esenciales o con fenotipos deletéreos ante niveles bajos de expresión tienden a mostrar menor ruido que el resto, mientras que genes asociados a respuesta a estrés muestran niveles de ruido superiores al promedio (Lehner, 2008; Newman et al., 2006). El ruido transcripcional adquiere una relevancia particular durante el desarrollo, donde la variabilidad célula a célula puede amplificarse en procesos de decisión de destino celular, y donde su correcta regulación (ya sea para amortiguarla o para aprovecharla) parece ser importante para el desarrollo normal (Arias y Hayward, 2006; Urban y Johnston, 2018).

Durante mucho tiempo, extender estos hallazgos a mamíferos a escala genómica resultó difícil, dado que estudiar el ruido de un gen endógeno, célula por célula, requería métodos gen-específicos que no escalaban bien a todo el genoma. El desarrollo de la secuenciación de transcriptoma en células únicas (*single-cell RNA sequencing*, scRNA-seq) cambió este panorama, permitiendo estimar simultáneamente los niveles de expresión de miles de genes en miles de células individuales. Usando estos datos, un estudio en células madre embrionarias de ratón encontró que la arquitectura del promotor está asociada al ruido a escala genómica: los promotores con islas CpG se asocian a bajo ruido, mientras que la presencia de una caja TATA se asocia a alto ruido, un patrón confirmado posteriormente en otros estudios con células de mamífero (Faure et al., 2017; Morgan y Marioni, 2018; Sun y Zhang, 2020).

Este enfoque, sin embargo, tiene limitaciones importantes. La medición de scRNA-seq está sujeta a un ruido técnico considerable, particularmente para genes de baja expresión, lo cual introduce incertidumbre adicional sobre el valor de ruido biológico que se intenta estimar. Además, estos estudios miden niveles de ARN en estado estacionario, que dependen no solo de la transcripción sino también de procesos posteriores, como la tasa de degradación específica de cada transcripto, por lo que no es posible atribuir toda la variabilidad observada exclusivamente al proceso transcripcional. Y quizás la limitación más relevante para esta tesis: el promotor endógeno de un gen actúa siempre dentro de un contexto genómico particular, moldeado por elementos regulatorios cercanos y por la conformación de la cromatina circundante, lo cual dificulta aislar el efecto específico de la secuencia del promotor de todos esos factores adicionales. A esto se suma que, en genes con múltiples promotores alternativos (como vimos en la sección anterior), los datos de scRNA-seq no permiten distinguir qué promotor está efectivamente dirigiendo la transcripción detectada.

Superar estas limitaciones requiere un abordaje distinto: uno que permita medir la actividad y el ruido de un promotor aislándolo de su contexto genómico habitual, y que además preserve la posibilidad de estudiar cada promotor alternativo de un gen por separado.

**6- Abordajes experimentales para el estudio del promotor basal**

Buena parte de la evidencia discutida hasta acá (la capacidad predictiva de la secuencia del promotor, su asociación con arquitecturas particulares, su relación con el ruido) proviene de datos observacionales: genomas humanos tal como son, con sus promotores insertos en su contexto genómico habitual. Esta clase de evidencia es poderosa para detectar asociaciones, pero limitada para establecer causalidad, ya que la secuencia de un promotor y su contexto genómico (cromatina, elementos regulatorios vecinos, el propio gen al que pertenece) no varían de manera independiente en la naturaleza. Para preguntar si una característica de secuencia *causa* un determinado nivel de actividad o de ruido, y no simplemente lo acompaña, hace falta poder variar la secuencia mientras se mantiene fijo todo lo demás.

Antes de describir cómo se mide la actividad de un promotor aislado de su contexto, conviene precisar cómo se mide su actividad *en* ese contexto, es decir, su uso endógeno. La técnica más extendida con este fin es CAGE (*Cap Analysis of Gene Expression*), que aprovecha la estructura *cap* presente en el extremo 5' de todo ARN mensajero maduro como ancla para secuenciar específicamente los primeros 25 a 27 pares de bases de cada transcripto. Al mapear estas lecturas cortas contra el genoma, es posible identificar la posición exacta del sitio de inicio de la transcripción con resolución de un nucleótido, y cuantificar la frecuencia de uso de cada uno a partir del número de lecturas que se originan en él (**Kanamori-Katayama et al., 2011**; **Takahashi et al., 2012**). Esta técnica fue la base del consorcio FANTOM, cuyo mapeo sistemático de promotores en decenas de tipos celulares y tejidos humanos constituye, hasta la fecha, el catálogo de referencia más utilizado en el campo (The FANTOM Consortium and the RIKEN PMI and CLST (DGT), 2014). RAMPAGE, una variante posterior, incorpora secuenciación de a pares para exigir que cada lectura provenga de un ADNc verdaderamente completo en su extremo 5', mejorando la especificidad de la identificación de TSS a costa de un protocolo algo más laborioso (**Batut y Gingeras, 2013**). Un enfoque distinto, PRO-seq, no secuencia el ARN maduro sino que marca y captura directamente el extremo 3' del ARN naciente asociado a moléculas de RNAPII activamente comprometidas con la transcripción, lo cual permite detectar eventos de inicio incluso cuando el transcripto resultante es inestable y nunca llega a acumularse lo suficiente para ser detectado por CAGE (**Kwak et al., 2013**). Estas tres técnicas, y otras derivadas de ellas, miden en definitiva lo mismo desde ángulos distintos: cuánto se usa un promotor en su locus natural, rodeado de todo su contexto genómico y regulatorio habitual. Esto las distingue conceptualmente de los MPRA, que miden la actividad de una secuencia de promotor aislada de ese contexto.

Los ensayos reporteros masivamente paralelos (*Massively Parallel Reporter Assays*, MPRA) fueron desarrollados con este objetivo. La idea general es sintetizar o clonar miles de secuencias candidatas, colocarlas en un contexto reportero común (típicamente, dirigiendo la expresión de un gen fácilmente cuantificable) e inferir la actividad de cada una a partir de la abundancia relativa de sus transcriptos, identificados mediante secuenciación masiva (**Melnikov et al., 2012**; **Patwardhan et al., 2012**). Este diseño permite evaluar, en un único experimento, el efecto de miles de variantes de secuencia sobre la actividad transcripcional, algo impensado con ensayos reporteros tradicionales de a uno por vez.

Existen distintas variantes de este enfoque, cada una con sus propios compromisos. STARR-seq, una de las más utilizadas para el estudio de enhancers, invierte el diseño clásico: en lugar de colocar la secuencia candidata corriente arriba del gen reportero, la coloca dentro de la región no traducida 3', aprovechando que muchos enhancers pueden actuar sobre su propio ARN mensajero y auto-transcribirse; esto permite ensayar fragmentos de ADN genómico fragmentado al azar, sin necesidad de síntesis dirigida, escalando el enfoque a genomas completos (**Arnold et al., 2013**). Tanto los MPRA clásicos como STARR-seq comparten, sin embargo, una limitación importante en su forma original: son ensayos episomales, es decir, las secuencias candidatas se introducen en la célula como parte de un plásmido que no se integra al genoma, y por lo tanto no adopta el empaquetamiento de cromatina característico del ADN genómico, y su número de copias por célula puede variar considerablemente y de forma no controlada.

Distintas estrategias intentaron resolver esta limitación integrando las construcciones al genoma. El desarrollo de lentiMPRA, que emplea vectores lentivirales para insertar cada secuencia candidata en el genoma de manera estable, mostró una correlación considerablemente mayor con anotaciones regulatorias basadas en cromatina endógena que su contraparte episomal, evaluada con la misma biblioteca de secuencias (**Inoue et al., 2017**; **Gordon et al., 2020**). Sin embargo, incluso en este caso, cada célula de la población recibe la construcción en un sitio de integración distinto y aleatorio del genoma, lo que introduce variabilidad adicional asociada al contexto cromatínico particular de cada sitio de inserción, un factor que resulta especialmente problemático si lo que se busca medir es, precisamente, cuánta variabilidad genera la propia secuencia del promotor.

Un enfoque distinto, más cercano al que se sigue en esta tesis, consiste en integrar cada secuencia candidata en un único locus genómico predefinido, idéntico para toda la biblioteca, eliminando por completo la variabilidad de número de copias y de contexto de inserción entre secuencias. Bajo este diseño, cualquier diferencia observada entre promotores puede atribuirse, con mayor confianza, a la secuencia misma y no al lugar del genoma donde terminó integrada. Un antecedente influyente en esta línea utilizó justamente un ensayo de integración basado en recombinasas para mapear la actividad autónoma de fragmentos aleatorios de ADN a lo largo de todo el genoma humano, encontrando que la capacidad de iniciar transcripción de manera independiente del contexto es una propiedad extendida, y no excepcional (van Arensbergen et al., 2017).

En paralelo al desarrollo de estas herramientas experimentales, ha crecido un campo dedicado a predecir la actividad transcripcional directamente a partir de la secuencia, mediante modelos computacionales entrenados con datos genómicos a gran escala (Agarwal y Shendure, 2020; Dudnyk et al., 2024). Es importante notar que estos modelos, en su gran mayoría, se entrenan justamente con datos de actividad endógena obtenidos mediante las técnicas recién descriptas: PUFFIN, por ejemplo, fue entrenado usando datos de CAGE (tanto de FANTOM como de ENCODE), RAMPAGE y PRO-cap combinados, aprendiendo a predecir cuánto se usa un promotor en su contexto genómico natural, y no cuánta actividad autónoma es capaz de generar aislado de él (Dudnyk et al., 2024). Estos modelos, en cierto sentido, formalizan la misma pregunta que venimos planteando desde el prólogo: cuánta información regulatoria puede extraerse de la secuencia por sí sola, aunque, a diferencia de los MPRA de locus único, lo hacen a partir de actividad medida en contexto, no en aislamiento.

Con estas herramientas es posible, entonces, retomar de manera experimental las preguntas que quedaron planteadas en las secciones anteriores: qué elementos de secuencia determinan la actividad y el ruido de un promotor, si esa determinación ocurre de manera directa o mediada por el contexto de cromatina que la propia secuencia genera, y si la multiplicidad de promotores de un gen refleja plasticidad regulatoria o simple redundancia. Los objetivos concretos de esta tesis, que se desprenden de estas preguntas, se detallan a continuación.

# Metodos Exp

**Materiales y metodos experimentales:**  
Descripción general del enfoque

Con el fin de cuantificar las propiedades transcripcionales que se desprenden intrínsecamente de la secuencia del promotor basal, diseñamos un ensayo reportero masivo en paralelo (MPRA, por sus siglas en inglés) que permite evaluar simultáneamente la media y el ruido transcripcional, similar a las técnicas de sort-seq utilizadas principalmente en levaduras y bacterias. Utilizamos una *library* con 23908  secuencias de ADN de 252pb que, en su gran mayoría, comprenden un TSS anotado y su región flanqueante (235pb río arriba, 16pb río abajo), abarcando los elementos esenciales del promotor basal (7,8). La misma fue producida previamente por el grupo de la Dra. Fiszbein (Boston University), en cuyo laboratorio, y en colaboración con la Dra. Uriostegui-Arcos, realicé la totalidad de los experimentos de la presente tesis. Estas secuencias se clonaron en un plásmido reportero, río arriba de la secuencia codificante para EGFP. Mediante el uso de la técnica de intercambio de casete mediado por recombinasa (RMCE), generamos líneas estables de HEK293T, donde cada célula posee una copia única del constructo conteniendo al gen EGFP bajo control de un promotor de la *library*, todas en la misma posición del genoma. 

Para obtener la distribución de niveles de expresión de cada promotor, las células se fraccionaron en siete poblaciones según la señal de fluorescencia de EGFP, se extrajo su ADN genómico (ADNg), y se analizó la abundancia de cada promotor en las distintas poblaciones mediante Amplicon-seq. A partir de estos datos, reconstruimos la distribución de la expresión de cada promotor para estimar su nivel de expresión media y la varianza asociada. El procedimiento se esquematiza en la Figura M1.

El procedimiento se realizó enteramente a partir de dos poblaciones distintas de células, generando dos réplicas independientes. Asimismo, para preservar la representatividad de la *library* y evitar cuellos de botella, los pasos de biología molecular y celular se ejecutaron con múltiples réplicas técnicas y volúmenes de trabajo superiores a los estándares.

Fig. M1 (pipeline exp)

Construcción de los plásmidos reporteros conteniendo la *library* de promotores basales

Los oligonucleótidos de cadena simple que conforman la *library* (cuya composición y detalle se especifican en la Sección Métodos Computacionales), se amplificaron inicialmente mediante PCR de extensión por solapamiento (*overlap extension PCR*) usando Platinum SuperFi II Master Mix (ThermoFisher, 12368010\) y los *primers* PromLib Forward y Reverse (ver Anexo Tabla 2\) durante 14 ciclos, seguidos de una purificación por columna (QIAGEN, 28104). 

El vector *backbone* (Fig. M2A) —término utilizado para este plásmido a lo largo de este trabajo— en el que se clonó la *library*, es una variante del reportero fluorescente bicromático (EGFP \+ DsRed) de *splicing* alternativo desarrollado por Orengo et al. (74). En esta variante, mutaciones específicas en el sitio de *splicing* 3' del exón alternativo aseguran la expresión exclusiva de EGFP. A su vez, el plásmido empleado comprende sitios LoxP no idénticos flanqueando la construcción reportera, que permiten su integración sitio-específica en células mediante un sistema de intercambio de casetes mediado por recombinasa (RMCE).

La elección de un vector derivado de un reportero de *splicing* alternativo obedeció a la disponibilidad y compatibilidad con líneas celulares desarrolladas previamente por nuestras colaboradoras (Fiszbein Lab, Boston University). En la Figura M2B se observan células HEK293T donde se evidencia la diferencia entre el reportero de *splicing* original y la versión optimizada para este ensayo, ambos bajo el promotor viral CMV. Además, la ausencia de expresión de dsRed fue validada previamente mediante RT-qPCR.

Para la integración de la *library* al plásmido *backbone,* se digirió a este último con las enzimas de restricción NdeI y NheI (ThermoFisher, FD0583 y FD0974), y se purificó el fragmento del vector lineal por electroforesis preparativa en gel de agarosa, utilizando un gel 1% y el kit QIAQuick Gel Extraction (QIAGEN, 28704). Con el objetivo de mejorar la eficiencia del clonado posterior, se utilizó el kit de purificación de productos de PCR QiaQuick (QIAGEN, 28104\) seguido por el kit de purificación MinElute (QIAGEN, 28004), que a su vez permite la concentración del fragmento en menor volúmen. La *library* se clonó en el plásmido digerido mediante ensamblaje Gibson (*Gibson assembly mix*, NEB, E2611L) en una relación molar inserto-vector de 3:1, con 0.250pmoles de ADN total. La reacción de ensamblaje Gibson se precipitó con isopropanol para remover las sales que puedan interferir con la electroporación posterior. Se introdujo el precipitado por electroporación en células MegaX DH10B T1R Electrocomp™ (ThermoFisher, C640003) para su transformación. Las bacterias fueron plaqueadas en placas cuadradas (245mm de largo) de LB ágar con ampicilina 100mg/L. Finalmente, los plásmidos se aislaron usando el kit de purificación de plásmidos PureLink™ Expi Endotoxin-Free Maxi (Invitrogen, A31217). 

Fig M2. ESQUEMA VECTOR \+ microscopia 

Para comprobar la correcta inserción del promotor, se realizó un ensayo de *colony* PCR (directo de bacterias, sin purificación del ADN) utilizando los *primers* PromLib. No se detectó la presencia del promotor CMV original (\~700 pb), mientras que las bandas obtenidas corresponden al tamaño esperado de la *library* (\~300 pb) (Fig. M3)

Fig M3. PCR integracion

Cultivo celular y transfección

Se utilizaron las células HEK293T-A2 (34), que ya han sido caracterizadas y poseen un único locus para RMCE en su genoma. Las mismas se cultivaron en medio DMEM (Dulbecco’s Modified Eagle Medium) con alto contenido de glucosa y piruvato (Gibco, 11965118), suplementado con suero fetal bovino al 10% (Gibco, A31406-02). Las células se mantuvieron en incubadora humidificada a 37°C con 5% de CO2.

Para generar líneas celulares estables con integración genómica del reportero, la *library* se co-transfectó con un 10% (m/m) de un plásmido que codifica la recombinasa Cre. Se transfectaron un total de 15 µg de ADN plasmídico en ocho placas de 10 cm utilizando Lipofectamine 3000 (ThermoFisher, L3000-015) y Opti-MEM (Gibco, 31985-070), siguiendo las instrucciones del fabricante. Después de la transfección, las células se seleccionaron con 2 µg/mL de puromicina (Gibco, A1113802) durante dos semanas. Posteriormente, las células se distribuyeron en dos grupos distintos, los cuales se trataron como réplicas biológicas independientes. 

Para las transfecciones individuales de los plásmidos reporteros, se empleó el mismo procedimiento escalando proporcionalmente las cantidades: 2 µg de ADN plasmídico en un único well de una placa de seis wells. Esto incluye tanto los controles con el plásmido *backbone* con CMV, tanto en su versión bicromática como mutada, como aquellos con promotores específicos de la *library* utilizados para el control de consistencia. La confección de estos últimos implicó en primera instancia la amplificación de los promotores desde el *pool* de la *library* (en Anexo Tabla 1 se presentan los *primers* correspondientes a cada uno), y su clonado por ensamblaje de Gibson en el plásmido *backbone* digerido.

Las células se visualizaron en un microscopio ECHO Revolve (RVL-100M) con un objetivo de 20X.

Citometría de flujo y *sorting*

Para el *sorting* de la *library* de promotores, se utilizó un equipo Beckman Coulter MoFlo Astrios de seis vías. Los datos de citometría de flujo se analizaron con FlowJo v10.10.0.

Para tener seguridad sobre el nivel de discretización suficiente para la distinción de los patrones de ruido y media que pudiéramos observar, se realizó un análisis de potencia considerando múltiples alternativas. Es importante comentar que los costos del ensayo crecen en gran medida con el aumento de las fracciones recolectadas, tanto por el tiempo de *sorting* creciente como los costos de secuenciación. Por lo tanto, la elección final resulta de un compromiso entre nuestra capacidad de medir las variables de interés (expresión media y dispersión) y el costo del experimento.

Fig M4. Potencia(hacer)

Se estableció que todas las fracciones tengan idéntico número de células para evitar sesgos diferenciales entre fracciones en el tratamiento posterior de las mismas. Para la determinación de los umbrales, se usó como base una citometría de flujo de las mismas células. La fracción \#6 sería aquella con menos células y por lo tanto es la que definiría el tiempo total (y con ello el costo) del *sorting*, y no poseería un límite superior. Se decidió por lo tanto que comprenda el 1% de la población total de células, estableciendo así la posición del umbral más alto. Respecto al menor, se lo ubicó a partir en del pico del control negativo de células HEK293T-A2 no transfectadas y de las celulas transfectadas(Fig. M5A). Una vez establecidos el límite inferior y el umbral más alto, se dividió el rango que queda en seis intervalos de tamaño constante en escala logarítmica (Fig. M5B). El diseño con intervalos equivalentes en tamaño es esencial para la correcta comparación entre fracciones. Como control positivo se emplearon células que portan el reportero bajo el promotor CMV. 

SPara la preparación del experimento en sí mismas, se analizaron las células HEK293T-A2 transfectadas con la *library* hasta recolectar un total de 250.000 células por gate, las cuales se preservaron en DNA/RNA Shield (Zymo Research, R1100-250) para los experimentos posteriores. Este procedimiento se realizó tres veces a partir de cultivos independientes de células transfectadas con la *library*. La decisión de recolectar idéntica cantidad de células para todas las fracciones, y que no sea proporcional al porcentaje de células en la población total, busca evitar sesgos en los procesos posteriores de extracción de ADN y PCR, principalmente. Eventualmente, sin embargo, esto llevará a distorsiones en las distribuciones reconstruidas de los promotores, que serán corregidas computacionalmente mediante un reescalado por la proporción real obtenida en el perfil de citometría previo  (ver Cuantificación de la actividad media y el ruido de los promotores de la *library)*

Fig M5. Figura gates

Controles de *spike-in* basados en células

Para evaluar posibles sesgos técnicos entre las distintas fracciones introducidas durante la manipulación de las muestras *post-sorting* (desde la extracción de ácidos nucleicos hasta la secuenciación), decidimos introducir moléculas identificables en cantidades conocidas (*spike-in*) a las mismas. La incorporación de *spike-in* es una práctica frecuente en análisis genómicos/transcriptómicos que implican la comparación cuantitativa precisa entre muestras procesadas por separado. Aunque el procedimiento estándar es el agregado de moléculas de ADN, ARN o células de otras especies, dadas las características del presente ensayo decidimos generar un control *spike-in* similar pero discernible al material con el que se iba a trabajar: células HEK293T-A2 con el reportero incorporado en el genoma. 

Brevemente, se clonaron tres secuencias de aproximadamente el mismo tamaño y contenido de GC que las de la *library* (ver Anexo Tabla 2\) en el vector *backbone* y se co-transfectaron con el plásmido de expresión de la recombinasa Cre en células HEK293T-A2, para obtener líneas estables, en la forma descrita anteriormente. Se agregaron 30, 300 o 3000 células de cada línea de *spike-in* respectivamente en cada tubo de recolección de FACS, juntándolas con las células separadas por el *sorter*.

Co-extraccion de ADN/ARN  y RT-qPCR de EGFP

### 

El ADNg y el ARN total se co-extrajeron de las fracciones obtenidas utilizando el *Quick-DNA/RNA Miniprep Kit* (Zymo Research, D7001). Mientras que la extracción del ADNg es central para el ensayo, el ARN se utilizó para evaluar la correlación entre el *sorting* asociado a la señal de EGFP y los niveles de transcripto de dicho gen. Dado que el objetivo es vincular la secuencia promotora con su actividad transcripcional a través de la intensidad de fluorescencia, resulta crucial verificar la correspondencia entre la señal proteica y la abundancia del ARNm de EGFP.

El ADNc se sintetizó a partir de las muestras de ARN de cada *bin* de células utilizando el *cDNA Synthesis Kit* (Thermo Scientific, K1622) siguiendo las instrucciones del fabricante. La transcripción reversa se realizó con 500ng de ARN total. Posteriormente, se determinaron los niveles de ARN de EGFP mediante PCR cuantitativa (qPCR) con primers específicos (ver Anexo Tabla 2), utilizando la mezcla de reacción  *Maxima SYBR Green/ROX qPCR Master Mix (2X)* (*Thermo Scientific,* K0222) en un sistema *ABI 7900HT Fast Real-Time PCR* (Applied Biosystems). Los niveles de expresión génica relativa se calcularon mediante el método 

2−ΔCt 

(75), donde 

ΔCt \= (CtEGFP \- CtGAPDH)

, usando GAPDH como control interno (Fig. M6). Se realizaron tres réplicas técnicas para cada determinación, utilizando el promedio de las mismas para el cálculo de la abundancia.

Fig. M6 QPCR 

En primera instancia, los resultados (Fig. M6) indican que la fracción 0 (correspondiente al nivel mínimo de fluorescencia, apenas por sobre el de las células sin transfectar) no presentaba resultados consistentes. Este fenómeno es atribuible a la dificultad de discernir la señal de EGFP de la autofluorescencia propia de las células para células con niveles muy bajos de expresión. Por consiguiente, se decidió excluir esta población de la secuenciación. Asimismo, se descartó la réplica 3 debido a que no presentó un patrón de qPCR consistente (Fig. M6, derecha) y mostró un rendimiento insuficiente en la extracción de ADNg en varias de sus poblaciones (no mostrado).

Adicionalmente, se evaluó la posible persistencia de células no transfectadas tras la selección con puromicina. Si bien estas células no generarían amplicones de la *library*, su acumulación en los gates de menor expresión podría afectar la representatividad del número de células positivas. Para detectar la correcta presencia de los promotores de la *library*, se utilizaron *primers* flanqueantes a los sitios LoxP de las células HEK293T-A2 (ver Anexo Tabla 2\) para realizar una PCR sobre el ADNg de las distintas fracciones, permitiendo discernir por tamaño la presencia del constructo de la *library* (Fig. M7). Aunque se detectaron células sin el reportero, su abundancia relativa fue uniforme en todas las muestras, por lo que no se consideró un factor de sesgo para los análisis posteriores.

Fig. M7 IMAGEN PCR\_loxP\_gates\_ctl

Secuenciación de amplicones (Amp-seq)

Para la secuenciación paralela masiva del ADNg de las distintas fracciones, se empleó un enfoque de PCR de dos pasos. Primero, se realizó una PCR primaria usando los *primers* PromLib Forward y Reverse (Tabla 1\) con *Platinum™ SuperFi II PCR Master Mix* (Thermo Fisher Scientific, Cat. No. 12-368-010), siguiendo las instrucciones del fabricante. Luego, se llevó a cabo una segunda PCR para incorporar las secuencias de los adaptadores de Illumina, manteniendo constante el *primer* *forward* con el índice i5 y usando un *primer* *reverse* con el índice i7 específico para cada muestra (ver Anexo Tabla S2), empleando nuevamente la *Platinum™ SuperFi II PCR Master Mix*. Los productos de la segunda PCR se purificaron mediante extracción de bandas del gel (QIAGEN, 28704\) y posteriormente se enviaron para la secuenciación de amplicones a MedGenome[^2]. Las secuenciacioónes se realizaronrealizó con una cobertura de 1000× utilizando lecturas *paired-end* de 150 pb en la plataforma *NovaSeq*. La incorporación de un índice i7 específico por muestra permite mezclar las muestras en una única corrida de secuenciación y evitar un posible *batch effect*. Sin embargo, cada réplica debió mantenerse en corridas en paralelo debido a la ausencia de más índices. 

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

Las lecturas *paired-end* de *Illumina* se procesaron inicialmente con *cutadapt* (56) para el filtrado por calidad y el recorte de *primers* y de calidad, y posteriormente con *FASTP* para el recorte de colas de poli-G y de baja calidad (77) (Fig. M8). Se descartaron aquellas lecturas que carecían de superposición de *primers*. El alineamiento contra las secuencias *FASTA* de la *library* se realizó mediante *HISAT2* (78). Únicamente se consideraron fragmentos con longitudes entre 230 y 270 pb, omitiendo los alineamientos secundarios. La abundancia de cada promotor en cada muestra se cuantificó con *SAMTOOLS* (1.20) (79). El control de calidad de las lecturas crudas se llevó a cabo con *FASTQC[^3]* y el desempeño del *pipeline* se reportó mediante *multiQC* (80). 

Fig. M8 pipeline bioinfo

Cuantificación de la actividad media y el ruido de los promotores de la *library*

Salvo que se indique lo contrario, los análisis estadísticos y el análisis exploratorio de datos se realizaron utilizando R 4.2.  
Los conteos de lecturas por promotor se corrigieron según la desviación del tamaño de cada *library* respecto al tamaño medio. Además, dado que el esfuerzo de muestreo para obtener 250.000 células varía para cada rango de intensidades de EGFP, los conteos también se corrigieron mediante un término de relativización determinado por:

${F}_{k}=\frac{{p}_{k}}{min(p)}$										(E1)

donde [![][image1]](https://www.codecogs.com/eqnedit.php?latex=p#0) es la proporción de células en la fracción [![][image2]](https://www.codecogs.com/eqnedit.php?latex=k#0) respecto a la población total, medida por citometría de flujo (Fig M8A). Es decir que las fracciones que, al medir por citometría de flujo, presentan mayor número de células, se les amplificará, en términos relativos y proporcionales, el número de *counts* por promotor observados. Esta corrección, más típica de la ecología que del campo de la biología molecular, es fundamental para comparar correctamente entre las fracciones de células y evitar un sesgo hacia los gates más altos, donde hay menos células en la población total. 

Se calcularon estadísticas descriptivas, como la media y la varianza, para cada secuencia que superara un umbral de 1000 conteos corregidos en cada réplica, asignando a cada fracción un valor de expresión coincidente su índice numérico (1 al 6). Las secuencias con una diferencia en la mediana superior a 2 entre ambas réplicas se descartaron por inconsistencia. Asimismo, se descartaron las secuencias con patrones altamente bimodales al considerarlas potenciales artefactos (Fig. M8B). Para detectar estos patrones, se ideó la estrategia de considerar aquellas con una varianza un 20% superior a la esperada para una distribución uniforme, dado un número [![][image3]](https://www.codecogs.com/eqnedit.php?latex=j#0) de *gates* con conteos. La varianza esperada de la distribución uniforme se calculó como:

Si [![][image4]](https://www.codecogs.com/eqnedit.php?latex=j#0) es impar:

[![][image5]](https://www.codecogs.com/eqnedit.php?latex=%5Ctext%7Bvar%7D\(j\)%20%3D%202%20%5Csum_%7Bi%3D1%7D%5E%7Bj%2F2%7D%20\(i%20-%200.5\)%5E2#0)

Si [![][image6]](https://www.codecogs.com/eqnedit.php?latex=j#0) es par:

[![][image7]](https://www.codecogs.com/eqnedit.php?latex=%5Ctext%7Bvar%7D\(j\)%20%3D%202%20%5Csum_%7Bi%3D1%7D%5E%7B\(j-1\)%2F2%7D%20\(i%20-%200.5\)%5E2#0)								(E2)

Fig M9. sample effort \+ bimodal

**Caracterización de la *library***

Composición y diseño de lLa *library*

La *library* utilizada consiste en 23908 secuencias de 300 pb de longitud, de las cuales 20851 incluyen una región de promotor basal (core promoter) de \-251 a \+16 con respecto a la posición anotada como TSS principal en la Eukaryotic Promoter Database[^4] (EPD) (73), con adaptadores de 24pb en cada extremo para facilitar la clonación (Fig. M10 A-B). EPD define sus TSS (29512 sitios en total) a partir de datos del repositorio del consorcio FANTOM que utilizan la técnica *Cap Analysis of Gene Expression* (CAGE) para la determinación del extremo 5’ de los transcriptos con la precisión de un nucleótido y con ello inferir la actividad promotora a lo largo de todo el genoma en diversas muestras biológicas. Agrupando únicamente las muestras de dicho repositorio correspondientes a cultivos primarios de células, se observa como el TSS indicado por EPD para los promotores seleccionados de la *library* coincide con aquel más usado, si bien son mínimos los casos donde no se observa actividad proveniente de las bases vecinas (Fig. M10C). Es preciso aclarar que en ningún momento se evaluó el TSS en el contexto del reportero utilizado, por lo que cada vez que me refiera a dicho sitio a lo largo de esta tesis, será teniendo en cuenta al anotado.

De las restantes secuencias presentes en la *library*, 2910 son regiones de enhancer de 152pb de longitud inmediatamente río arriba de 100pb del promotor FN1(Fig. M10 A-B). También hay un subgrupo de 147 promotores asociados al cáncer, tanto en sus versiones wild-type como mutantes, y que no responden a la estructura tal cual fue definida para los promotores EPD. Estos últimos grupos, aunque presentes en la *library* y secuenciados, no fueron considerados durante el análisis bioinformático posterior al recuento de actividad.

Fig. M10   
   
Elementos de secuencia del promotor basalCaracterización de los promotores

La base de datos EPD provee información respecto a la presencia de motivos típicos en sus promotores anotados, a partir de la búsqueda de patrones y su localización respecto al TSS: TATA-box, CCAAT-box, GC-box e INR. Según sus datos, del 23851 secuencias, el 47.7% contiene un GC-box, el 32% un INR fuerte, el 16.3% un CCAAT-box y solo el 8% TATA-box (Fig. M11A).

Fig. M11

La presencia de islas CpG fue determinada a partir de las anotaciones provenientes de UCSC, considerando aquellos casos donde se superponen más del 100pb con la secuencia promotora en cuestión. Las islas CpG son bastante prevalentes en promotores, observándose en el 64% de los mismos (Fig. M11B). También se buscó el solapamiento con elementos repetitivos y/o transponibles, anotados en la base de datos RepeatMasker[^5]. Además de las clases de transposones presentes en dicho repositorio, se incorporó una clase denominada “Repeticiones de Baja Complejidad”, agrupando las categorías *Simple Repeats* y *Low Complexity*, asi como una de “Elementos transponibles”. 

La frecuencia nucleotídica y su identidad en sitios específicos fue evaluada con el paquete de R-Bioconductor Biostrings[^6]. Se evaluó el contenido de G y C (Fig. M11C), así como la identidad del dinucleótido del TSS (Fig. 11D). Los patrones buscados de forma estricta fueron YCASW para el INR “fuerte”, TCT para el clásico motivo de proteínas ribosomales y el dinucleótido YR (PyPu)[^7]. Respecto a este último patrón se discernir también en sus cuatro posibilidades (CA, CG, TA, TG). Aquellos promotores que no cuadran en su TSS con alguno de los patrones mencionados, fueron catalogados como “No canónicos”. Se incluye el dinucleótido GC para evidenciar que la prevalencia del CG no es simplemente producto de alto contenido de dichos nucleótidos.  

Patrones de conservación de los promotores basales

A su vez, se incorporaron datos que reflejan la historia evolutiva de los promotores, tanto a nivel de secuencia como funcional. Por un lado, se extrajeron datos de PhyloP score (81) provenientes de la comparación entre el genoma hg38 y 100 especies de mamíferos. Evaluando la mediana de dicho valor del “metapromotor” a cada base (Fig. M12A), se observa claramente una mayor conservación en la región proximal al TSS (+16 a \-50), con claros picos alrededor del \+1 y del \-30, asociado al TATA-box. En una región de intermedia cercanía (-50 a \-150) hay un progresivo decaimiento de la conservación, mientras que se acerca mucho a valores de evolución neutra para la región más distal (-150 a \-235). A su vez, para cada una de estas regiones, en cada promotor, se evaluó el PhyloP score medio. En términos de conservación funcional, nos basamos en datos de Young et al. (19), quienes utilizan datos de actividad promotora en tejidos de humano y ratón para considerar si, las secuencias que se pueden considerar homólogas, están activas en ambas especies, si perdieron actividad promotora en humanos o en ratón o si, por el contrario, la adquirieron en alguna de estas especies (Fig M12B). 

Fig. M12. evo library

Patrones de actividad endógena de los promotores

Se obtuvieron datos de accesibilidad de cromatina en HEK293 provenientes de ENCODE (accesible como ENCSR956YZJ). La presencia de Módulos Regulatorios en *cis* (CRM) se determino con datos de Remap basados en cientos de muestras de ChIP-seq [(44, ver](https://www.zotero.org/google-docs/?sWJeqi) Análisis de datos de ChIP-seq[)](https://www.zotero.org/google-docs/?sWJeqi). La posición de los enhancers anotados se obtuvo de UCSC. 

La actividad promotora endógena fue cuantificada a partir de datos de CAGE, provistos por el consorcio FANTOM5. Estos datos se utilizaron tanto para la determinación de la actividad en muestras específicas, así como para obtener valores del comportamiento de los promotores a lo largo de todas las muestras, como la especificidad de tejido del promotor o su forma (la distribución de sus TSS). Para las muestras específicas (HEK293, HeLa, músculo esquelético), se utilizó el paquete de R CAGEr v2.12.0 [(83)](https://www.zotero.org/google-docs/?OtTt7L). Específicamente, aplicamos una normalización de *power-law* para calcular los niveles de expresión en *Tags* Por Millón (TPM), filtrando aquellas señales que estuvieran por debajo de un umbral de 1 TPM. Los clusters de TSS se definieron mediante el algoritmo *paraclu* (84). Para refinar la señal, se conservó el 80% central del rango de expresión y se intersectó con las coordenadas genómicas de la *library* para evaluar la actividad endógena específica de cada promotor.

Si bien esta metodología ofrece una mayor precisión, su alta demanda computacional (en términos de tiempo y memoria) la hizo inviable para evaluar la actividad de los promotores en cada una de las muestras de células primarias y tejidos del dataset del consorcio FANTOM5. Dado que la representación de tipos celulares y tejidos en el dataset es muy heterogénea, decidimos agruparlos por ontología utilizando las categorías definidas por Andersson et al. 2014 (85). Las muestras se fusionaron, se normalizaron por el tamaño de la *library* y se restringieron a los rangos genómicos de la *library*; finalmente, la actividad de los promotores en cada grupo se evaluó como la suma de los *counts* que solapaban con el rango genómico de cada promotor. Cabe aclarar que, si bien la base de datos de FANTOM5 incluye otras categorías de muestras humanas, como por ejemplo líneas celulares,  Para distinguir entre promotores *housekeeping* y aquellos con alta especificidad de tejido, utilizamos el índice de Gini de la actividad del promotor a lo largo de todos los grupos de ontología y dividimos a los promotores en terciles. Los promotores inactivos se excluyeron de este análisis (menos de 1 TPM en cualquier grupo) y se etiquetaron como "Sin actividad detectable". La fórmula para el índice de Gini, luego de ordenar cada uno de los *n* grupos de ontología por actividad decreciente, es:  
$G=1-\sum\limits_{i=1}^{n}({p}_{i}-{p}_{i-1})({q}_{i}+{q}_{i-1})$							(E3)  
donde, en este caso, *pi* y *qi* son las proporciones acumuladas de los grupos de muestras y de la actividad promotora, respectivamente.

La forma de los promotores se calculó a partir de la combinación de todas las muestras de FANTOM5 analizadas (células primarias y tejidos). En este caso, también se utilizó el paquete *CAGEr* y se aplicó una normalización de *power-law*. Los *clusters* de TSS se obtuvieron mediante el método *distclu*, con una distancia máxima entre TSSs de 5 y un mínimo de 10 *counts* por TSS. Se utilizó el ancho intercuantílico (0.05-0.95) como *proxy* de la forma. Los valores se dividieron en terciles de ancho, conservando el primer y tercer tercil como promotores focalizados y anchos, respectivamente. En caso de que se obtuvieran múltiples *clusters* de TSS para una misma región promotora, solo se utilizó aquel que solapaba con la posición establecida por la EPD.

Coocurrencia de las características de los promotores basales

Las características de los promotores, tanto aquellas basadas en la secuencia como aquellas que surgen de estudiarlos en sus contextos endógenos, no son completamente independientes entre sí (Fig. M12). Esto implica que frecuentemente, sea complejo poder asignar a un efecto observado, una característica particular, con confianza de que no se trate de un efecto confusor de otra característica con alto grado de coocurrencia. Si bien esto se puede resolver en ciertos casos con una estratificación por a potencial característica confusora, en los casos más extremos de co-presencia, esta tarea resulta prácticamente imposible y es un limitante en este tipo de enfoques experimentales, basados en secuencias naturales.  

Fig. M12 co-ocurrencia

Algoritmos de predicción de la actividad promotora basados en la secuencia

Puffin es un modelo interpretable de aprendizaje automático que predice la señal de inicio transcripcional a resolución de base a partir de la secuencia del promotor, descomponiendo su predicción en la contribución aditiva de tres tipos de elementos de secuencia: motivos, iniciadores y trinucleótidos (85). A diferencia de otros modelos de *deep learning* de tipo "caja negra", Puffin permite atribuir la señal predicha a posiciones y motivos específicos de la secuencia, facilitando la interpretación mecanística de los determinantes de la actividad promotora. En este trabajo se utilizó un modelo pre-entrenado con datos de CAGE en humano (GRCh38).

El modelo predice actividad a partir de regiones de al menos 650 pb. Dado que los promotores de la *library* poseen únicamente 252 pb, resultó imprescindible incluir la región lindante en su contexto de inserción, incorporando los *primers* y parte del gen EGFP hacia uno de los extremos. En el extremo opuesto, el promotor se encuentra próximo al sitio LoxP y, dado que se desconoce la identidad del sitio de inserción en el genoma, no se dispuso de esas bases para completar el tamaño de secuencia requerido. Para compensar esta limitación, se generaron 20 secuencias aleatorias que se añadieron a cada promotor evaluado, y la predicción final para cada uno se promedió entre estas 20 variantes. Si bien se evaluaron varias métricas de salida del modelo, para estimar la actividad promotora se utilizó el *score* obtenido en la posición 0 de cada promotor, que en términos generales coincide con el TSS de mayor actividad en los contextos endógenos evaluados por CAGE, según EPD.

Adicionalmente, se incorporó la métrica de selectividad al contexto genómico ("*selectivity*"), reportada originalmente por Dudnyk et al. (85) a partir de predicciones con Puffin-D, una variante del modelo orientada a la predicción cuantitativa de expresión que, a diferencia de Puffin, admite hasta 100kb de secuencia como input, permitiendo así incorporar el contexto genómico circundante al promotor. Esta métrica cuantifica el grado de variación en el nivel de expresión predicho de un promotor al insertarlo *in silico* en miles de ubicaciones genómicas distintas, reflejando así su dependencia del contexto regulatorio circundante. En nuestro caso, no se recalculó esta métrica, sino que se utilizaron los valores ya reportados por los autores para el subconjunto de secuencias coincidentes con los promotores de nuestra *library*. 

Análisis de datos masivos de ChIP-seq

Se obtuvieron picos de ChIP-seq no redundantes de la base de datos ReMap y se intersectaron con la *library* de promotores para determinar el estado de unión de cada proteína. La intersección positiva de al menos un pico, fue evidencia suficiente para considerar la presencia de un CRM. El análisis se restringió a los factores de transcripción (TFs) que presentaron picos en al menos 100 promotores en ambas réplicas biológicas.  
De manera paralela, se llevó a cabo un análisis para modificaciones de histonas y marcas epigenéticas obtenidas del conjunto de datos de ChIP-Atlas (58).

**Evaluación del efecto de las características del promotor sobre los patrones transcripcionales medidos**

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

Análisis de enriquecimiento funcional de la unión de TFs

Las asociaciones entre la ocupación proteica de TFs y tanto el ruido como la actividad transcripcional se evaluaron utilizando el mismo marco de trabajo aplicado a todas las características binarias.   
Se realizó un análisis de enriquecimiento de conjuntos de genes (*GSEA*) (89) sobre el subconjunto de factores de transcripción identificados utilizando términos de *Gene Ontology* (*GO*), ordenados ya sea por el tamaño del efecto del test de suma de rangos de *Wilcoxon* (*Wilcoxon rank-sum test*) o por el *AUC-ROC*.Estos análisis se implementaron mediante *scripts* propios utilizando los paquetes de *R/Bioconductor* *AnnotationDbi*, *clusterProfiler* (90) y *enrichplot*.

**Promotores alternativos**

Clasificación de promotores alternativos

Para clasificar los promotores de acuerdo a su relación con otros regulando el mismo gen, unimos todos los promotores anotados en la base EPD junto con los datos de CAGE de FANTOM5 (tejidos y cultivos primarios).  Como ya se mencionó, las muestras de FANTOM5 fueron unificadas en base a su ontología para evitar redundancias y normalizadas por el tamaño de la *library*. En primer lugar, se clasificaron como “Sin actividad detectada” a aquellos promotores que no alcanzaran 1 TPM en ninguna muestra. Si bien podría llamar la atención que dichos promotores estuvieran incluidos en la base de datos de EPD, que también está basada en datos de FANTOM5, esto asumo que se debe a que son promotores que se registraron activos en líneas celulares, que estoy excluyendo del presente análisis. Se agruparon los promotores restantes por gen, y se separaron aquellos ”Promotores únicos”, para los genes sin promotores alternativos con actividad detectable.   
Para aquellos genes con múltiples promotores, se buscaron aquellos con mayor actividad para clasificar como “Promotores principales”. Para ello se consideraron aquellos promotores que:  
1- Concentran el máximo porcentaje de counts normalizadas sumando todas las muestras.   
2- Contienen el máximo número de muestras con la actividad más alta.   
Un pequeño grupo de promotores cumplió únicamente una de aquellas condiciones, que fueron considerados “No clasificables”. El resto se clasificaron como “Promotores secundarios”.  
En segunda instancia cada par Principal-Secundario fue clasificado en base a la correlación en sus actividades endógenas. Para ello, se inició por la identificación de casos de “Alternancia” (o *switch*), término con el que nos referimos a muestras donde se observa una clara alternancia en el rol Principal/Secundario de los promotores. Se definió como tal cuando el Secundario supera los 5 TPM y, o bien la relación de actividad en escala logarítmica fuera 50% mayor para el Secundario, o bien el Principal no tuviese actividad detectable en la muestra (\<1 TPM). Los pares de promotores sin un caso de alternancia, se clasificaron como “Correlacionados” o “Independientes” según si la correlación (o bien de Spearman o de Pearson) fuese mayor a 0,5. 

# Resultados 1

Estimación masiva de las propiedades transcripcionales de promotores basales humanos

Una vez obtenidos los datos de secuenciación, en sus dos réplicas, el primer paso fue el alineamiento de las lecturas obtenidas con las secuencias de la *library*. Como se ha mencionado previamente, una preocupación a lo largo del proceso experimental fue evitar un cuello de botella que implique que, al analizar los datos, las lecturas observadas correspondieran todas a un subgrupo pequeño de las  secuencias totales. Afortunadamente, este no fue el caso, y contamos con lecturas del 80,9% de las secuencias, en al menos una réplica, y del 67,3%, en ambas réplicas (Fig. R1.1 venn) . 

Fig. R1.1 venn  rep1-rep2-tot.library

Si bien los análisis a realizar posteriormente serán todos a nivel comparativo entre las secuencias de las que podemos extraer datos confiablemente, nos planteamos la posibilidad de que haya, más allá de un cuello de botella aleatorio, un sesgo en la representatividad de las secuencias en las células. Comprender este aspecto nos serviría tanto para poder dimensionar las limitaciones en la extrapolación de los datos a la totalidad de los promotores anotados en humanos, así como para identificar puntos clave a mejorar en la metodología experimental a futuro. Con este objetivo, tomamos datos de un experimento previo realizado en el laboratorio de la Dra. Fiszbein con el reportero fluorescente bicromático (Fig. M1A), en el cual utilizó la misma *library* con los adaptadores ya ligados, y con un proceso idéntico al ya expuesto, generó células HEK293T-A2 con una variante del reportero integrado. En dicho ensayo, a diferencia del nuestro, se realizó una secuenciación de los promotores en dichas células previo a cualquier *sorting* celular o filtro por expresión. Encontramos una fuerte asociación entre los promotores seleccionados en ambos ensayos (OR \= 13.7, p \< 2.2e-16, *Fisher’s exact test*) (Fig. R2A venn\_presort). Esto indicaría que la incorporación de una secuencia a las células en el ensayo previo aumenta drásticamente la probabilidad de ser detectado en el ensayo actual, sugiriendo un sesgo común. Ante esto, surgen dos hipótesis: o bien hay un diferencia de partida en la representatividad de las secuencias, o alguno/s de los pasos experimentales presentan un sesgo sistemático.

Una hipótesis en este último sentido es que el contenido de G y C en la secuencia podría tener asociado pequeños cambios de eficiencia en la amplificación de los promotores por PCR. Evidentemente, las secuencias no observadas en nuestros resultados, tienen valores de contenido G+C más extremos, en ambos sentidos (Fig. R1.2B gc\_bias\_violin). A su vez, este patrón se repite si se tiene en cuenta la distribución en el contenido de G+C en todas las lecturas: al comparar con la distribución hipotética e ideal en la que todos los promotores tuvieran igual representatividad en los resultados, se evidencia una depleción en promotores de alto y bajo contenido de G+C (Fig. R1.2C gc\_bias\_violin). Si bien en las amplificaciones por PCR se utilizó un *kit* que ha sido probado como eficiente para proporciones extremas de AT y GC en los amplicones, las eficiencias para dichos casos podrían diferir, evidenciando, luego de muchos procesos amplificadores, el patrón que observamos. Una posible solución para ensayos futuros podría ser la combinación de productos de PCR con protocolos optimizados para distintas proporciones de AT/GC. De cualquier manera, este sesgo, si bien limita el universo que podrá ser abordado en el presente trabajo, no afecta los resultados internos del mismo.

Fig. R2 Sesgo de representatividad

Dado que una parte esencial del análisis implica la comparación precisa entre las muestras correspondientes a los distintos *bins* de células por expresión, realizamos un control por *spike-in* celulares (ver Métodos). Evaluamos la abundancia de los *counts* (normalizando únicamente por *library size*) de las tres secuencias incorporadas post-*sorting* (Fig. R1.3). Afortunadamente, encontramos una relación relativamente constante para cada una de ellas, sin grandes desvíos en ninguna muestra que nos hubieran hecho plantearnos la posibilidad de reescalar los datos o tener que descartar alguna muestra.

Fig. R3. Spike-in

A partir de los datos de secuenciación obtenidos, en sus dos réplicas, y luego de los ya mencionados controles y procesamiento de los datos, se logró la reconstrucción de las distribuciones subyacentes de expresión asociadas a cada uno de las secuencias regulatorias. Esto se traduce visualmente en histogramas indicando la frecuencia de lecturas normalizadas con que se observó un cierto promotor en cada fracción de células (Fig R1.4). Fue posible con estos datos obtener la media y la dispersión de los conteos, indicadores de actividad promotora y ruido transcripcional respectivamente (Tabla 3, Anexo).

Fig. R1.4 ejemplo histograma

Para evaluar el grado de concordancia entre las distribuciones discretas reconstruidas y las distribuciones subyacentes, se seleccionaron una serie de secuencias de la *library* para replicar el procedimiento experimental pero con promotores individuales aislados. En estos casos, se generaron líneas celulares estables con el reportero integrado y regulado por cada una de dichas secuencias, a las que se le midió la fluorescencia de EGFP en cada célula por citometría de flujo (Fig. R1.5), con tres réplicas técnicas. 

Fig. R1.5 densities citometria (todos)

La cuantificación en este trabajo habría sido más precisa posiblemente si hubiéramos podido comparar la distribución de todos los promotores por citometría de flujo, sin la necesidad de fraccionar y secuenciar, que generan una métrica más indirecta y discretizada, tal como se puede observar en la Fig. R1.6A. Sin embargo, la masividad de los promotores a evaluar, entre otros factores, fueron un claro freno a esta metodología, por lo que las comparaciones entre las distribuciones medidas tienen un formato ejemplificado en la Fig. R1.6B. Pero contar con un cierto número de secuencias evaluadas por citometría, permite evaluar la consistencia de nuestra métrica de actividad: en la Fig. R1.6C se observa la asociación obtenida entre dichas métricas en el contexto de medición específica (por citometría) y masiva (por fraccionamiento y secuenciación). La correlación positiva, afortunadamente, es consistente con lo esperado. A nivel ruido, sin embargo, las métricas no son lo suficientemente robustas como para permitir una buena determinación a nivel individual de los promotores y replicar esta prueba de consistencia.

Fig. R1.6 \- comparacion histo-density \+ tendencia media vs media

La correlación entre las medias obtenidas para cada promotor en ambas réplicas tiene un coeficiente de Pearson de 0.72 (valor-p\<2e-16) (Fig. R1.7A mean\_replicates\_prefilter.jpg), contabilizando solamente aquellas secuencias con un umbral mínimo de *counts* y presentes en ambas réplicas, reflejando la robustez técnica del experimento realizado. Una vez incorporados los filtros de consistencia a nivel mediana y de exclusión de los patrones bimodales, que serán tenidos en cuenta para los análisis subsiguientes, se llega a un coeficiente de 0.85 (valor-p \< 2e-16)  (Fig. R1.7B mean\_replicates\_postfilter.jpg). En este último subconjunto de secuencias, de alta confianza, se encuentran los 12214 promotores que serán utilizados en los siguientes análisis, presentándose 10401 de ellos en ambas réplicas (Fig. R1.7C fig venn o bar post filter). 

En términos de ruido, tal como es esperable para una métrica de dispersión, la correlación entre la varianza medida entre ambas réplicas fue considerablemente menor, de 0.31 (coeficiente de Pearson, valor-*p*\<2e-16) luego de los filtrados (Fig. R1.7D). De las distribuciones observadas en las Fig. R1.7, se desprende que las poblaciones de ambas réplicas presentan distribuciones distintas a lo largo de las subpoblaciones predefinidas. La distribución de un promotor es en última instancia una relación respecto a la distribución de todos los promotores de esa réplica, y sus valores de media y ruido obtenidos no tienen mucho sentido sin dicho contexto. Es por ello que se consideró, para la evaluación de efectos estadísticos en los análisis subsiguientes, a cada réplica por separado y de forma independiente. 

Fig. R1.7 Filtros

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

[^1]:   Esta misma lógica fue llevada, mucho más recientemente, a una formulación cuantitativa por el físico Jeremy England, quien propuso que los sistemas de materia sometidos a un flujo externo de energía tienden, con el tiempo, a reorganizarse de manera que absorben y disipan esa energía cada vez más eficazmente —un fenómeno que llamó adaptación disipativa (England, 2015). Bajo esta lente, la aparición de estructuras capaces de sostener orden lejos del equilibrio, como las que caracterizan a la vida, dejaría de ser una rareza estadística para volverse, en cierto sentido, una tendencia esperable de la física de los sistemas dirigidos por energía. 

[^2]:  https://diagnostics.medgenome.com/research-service/

[^3]:  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[^4]:  https://epd.expasy.org/epd/

[^5]:  https://www.repeatmasker.org/

[^6]:  https://bioconductor.org/packages/release/bioc/html/Biostrings.html

[^7]:  Las bases subrayadas refieren a la posición del TSS anotado

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAgAAAAJBAMAAAD9fXAdAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVKvN74lEInYy3WYQmbv8EmWgAAAAEUlEQVR4XmP8z8DAwMRAHAEALZYBEce9Dw0AAAAASUVORK5CYII=>

[image2]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAALBAMAAABBvoqbAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAECJEVHaJq7vd75nNZjIrqulnAAAAEklEQVR4XmP8z8DwkYkBCEgiAGhhAgZfdY9rAAAAAElFTkSuQmCC>

[image3]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAANBAMAAACX52mGAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVLvvmRCrIt0yRGbNdol/YymBAAAAEklEQVR4XmP8z8DwkYkBCMgkAHy7AgpX6bqBAAAAAElFTkSuQmCC>

[image4]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAANBAMAAACX52mGAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVLvvmRCrIt0yRGbNdol/YymBAAAAEklEQVR4XmP8z8DwkYkBCMgkAHy7AgpX6bqBAAAAAElFTkSuQmCC>

[image5]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAKEAAAAwBAMAAACLRQdjAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMARO/dzburiWYQVHaZMiKgQ028AAAAOklEQVR4Xu3MoRHAIADAQGD/ZZkAfC1vepeXEZlnWHt9y7OORkejo9HR6Gh0NDoaHY2ORkejo/GH4wVwAAJQRw8HJwAAAABJRU5ErkJggg==>

[image6]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAcAAAANBAMAAACX52mGAAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMAVLvvmRCrIt0yRGbNdol/YymBAAAAEklEQVR4XmP8z8DwkYkBCMgkAHy7AgpX6bqBAAAAAElFTkSuQmCC>

[image7]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAALMAAAAwBAMAAACoHla2AAAAMFBMVEX///8AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAv3aB7AAAAD3RSTlMARO/dzburiWYQVHaZMiKgQ028AAAAOUlEQVR4Xu3MoQEAIACAMPX/Z71A+7I2FgnMMz7Zy/JOa7RGa7RGa7RGa7RGa7RGa7RGa7RGa3xcX/9kAlCB4WOCAAAAAElFTkSuQmCC>