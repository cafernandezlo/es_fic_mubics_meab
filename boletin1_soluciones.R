# Ejercicios resueltos Boletín 1

# carga de datos 
# disponibles en el siguiente enlace: https://udcgal-my.sharepoint.com/:u:/g/personal/carlos_fernandez_udc_es/IQDA-m4KvHdeRpYEcFmGcSJBAWjyzYRBmP05IaUZTMKO0s4?e=Tzs2fM

load("~/datos/GSE21779/gse21779.rda")
meta <- readRDS(file = '~/datos/TCGA-BRCA/meta.RDS')
data.conteos <- readRDS(file = '~/datos/TCGA-BRCA/data.conteos.RDS')

############################################################
#1. Con los datos gse21779 muestra las cinco primera filas y las cinco primeras columnas.
############################################################

# La función `exprs()` extrae del objeto `AffyBatch` la **matriz de intensidades crudas**: filas = sondas individuales del chip, columnas = muestras (archivos CEL). Indexamos `[1:5, 1:5]` para ver un subconjunto manejable. `intensity()` es sinónimo de `exprs()` para objetos Affy. `geneNames()` devuelve los identificadores de probe sets (los nombres de fila tras el resumen, 54675 probe sets).

exprs(gse21779)[1:5, 1:5]
intensity(gse21779)[1:5, 1:5]

############################################################
#2. ¿Qué son los valores que aparecen en las nombres de columna? Ejemplo GSM542488.CEL.gz (revisa lo que descargaste si tienes dudas)
############################################################

# Los nombres de columna son los **nombres de los archivos CEL** descargados de GEO. Cada archivo CEL corresponde a una muestra biológica
# escaneada con el chip Affymetrix. El nombre se descompone así:
# - **GSM542488**: identificador único de muestra en GEO (*GEO Sample*). Permite localizar la muestra en [https://www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/).
# - **.CEL**: formato propietario de Affymetrix que almacena las intensidades de fluorescencia píxel a píxel de cada celda del chip.

colnames(exprs(gse21779))[1:5]

############################################################
#3. ¿Qué es lo que hay en las filas? ¿Cómo puedes obtener el nombre? Saca los 5 primeros. Ejemplo `1007_s_at` ¿Qué tamaño tiene? 54675 
############################################################

# Las filas de `exprs(gse21779)` (en su estado crudo, antes de normalizar) contienen **sondas individuales**.Pero los **nombres de fila accesibles via `featureNames()`** son los **probe set IDs**: identificadores como `1007_s_at` que agrupan ~11-20 sondas del mismo gen. El chip hgu133plus2 tiene **54675 probe sets**.

#Los sufijos del identificador tienen significado:

#| Sufijo | Significado |
#|--------|-------------|
#| `_at` | Sonda estándar anti-sense |
#| `_s_at` | Reconoce múltiples transcritos del mismo gen |
#| `_x_at` | Posible hibridación cruzada con genes parálogos |
#| `_a_at` | Representación alternativa del mismo gen |

# Primeros 5 probe set IDs
featureNames(gse21779)[1:5]

# Número total de probe sets
length(featureNames(gse21779))

############################################################
# 4. Mapea el nombre de las sondas a genes, por ejemplo usa Gene symbol (nombre oficial del gen) y EntrezID (identificador numérico único asignado a cada gen dentro de la base de datos de Entrez Gene del NCBI). **Ayuda**, busca en Bioconductor el paquete que te coincida con la salida de `annotation(gse21779)`. Imprime por pantalla los 6 primeros elementos de un dataframe (llámalo gene_info) que tenga tres columnas PROBEID, SYMBOL y ENTREZID. Después mira a ver cuántos NAs hay en cada columna
############################################################

# PISTA: ¿te ha salido este mensaje? `'select()' returned 1:many mapping between keys and columns`. 
# Piensa un poco por qué, pero ahora tienes un gene_info con 57151 filas y se supone que estás mapeando el nombre de las sondas a genes. 
# un mismo `PROBEID` se asocia con múltiples `SYMBOL` o `ENTREZID`.
# Para resolverlo, tendrías que decidir cómo quieres manejar estas relaciones múltiples
# No hay una única opción, depende de las necesidades de cada momento.
# Además, tienes 10025 sondas que no tienen nombre del gen ni entrezID. Obten las primeras 15 sondas que no tienen valor en SYMBOL, verás que el segundo es la sonda 1552563_a_at que no se corresponde a ningún gen, sino a C8orf6 (chromosome 8 open reading frame 6).

gse21779  # hgu133plus2
# https://bioconductor.org/packages/release/data/annotation/html/hgu133plus2.db.html

pacman::p_load(hgu133plus2.db)

probe_names <- featureNames(gse21779)

# Mapeamos probe sets → SYMBOL + ENTREZID
gene_info <- select(hgu133plus2.db,
                    keys    = probe_names,
                    columns = c("SYMBOL", "ENTREZID"),
                    keytype = "PROBEID")
detach("package:hgu133plus2.db", unload = TRUE)

# Primeras 6 filas
print(head(gene_info))

#El paquete de anotación que coincide con `annotation(gse21779)` = `"hgu133plus2"` es # #**`hgu133plus2.db`** de Bioconductor. La función `select()` realiza la consulta. El mensaje #`'select()' returned 1:many mapping` indica que algunas sondas mapean a **más de un gen**: el #resultado tiene **57151 filas** aunque solo hay 54675 probe sets, porque algunos probe sets están #anotados a varias entidades. Las **10025 sondas sin SYMBOL ni ENTREZID** corresponden a: sondas de #control interno, regiones intergénicas, ORFs no caracterizados y lncRNAs descubiertos después del #diseño del chip.

# Dimensiones: observamos las 57151 filas por el mapeo 1:many
dim(gene_info)

# Conteo de NAs por columna
na_counts <- colSums(is.na(gene_info))
print(na_counts)

# Primeras 15 sondas sin anotación SYMBOL
library(dplyr)
na_symbol_info <- gene_info %>%
  filter(is.na(SYMBOL)) %>%
  head(15)
print(na_symbol_info)

# Sonda 
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

probe_1552 <- getBM(attributes = c('affy_hg_u133_plus_2', 'chromosome_name', 'start_position', 'end_position'),
      filters = 'affy_hg_u133_plus_2',
      values = "1552563_a_at",
      mart = ensembl)
print(probe_1552)
detach("package:dplyr", unload = TRUE)

############################################################
# 5. Busca aquellas filas que tengan duplicados. Utiliza dplyr
############################################################

# Un gen puede aparecer en múltiples filas de `gene_info` porque el chip incluye **varias sondas para el mismo gen** (diferentes exones, isoformas, o redundancia de control de calidad). Agrupamos por SYMBOL y filtramos los grupos con `n() > 1`.

pacman::p_load(dplyr)

duplicated_gene_info <- gene_info %>%
  group_by(SYMBOL) %>%
  filter(n() > 1) %>%
  ungroup() 

# Primeras filas para ilustrar
print(head(duplicated_gene_info, 20))

# Número de filas duplicadas
nrow(duplicated_gene_info)

detach("package:dplyr", unload = TRUE)

############################################################
# 6. En el ejercicio anterior verás que el gen `DDR1` es uno de los muchos, muestra todas las filas de `gene_info` que contengan en la columna `SYMBOL` el gen `DDR1`
# PISTA: hay cuatro sondas (1007_s_at, 207169_x_at, 208779_x_at, 210749_x_at) para el mismo gen (DDR1)
############################################################

# DDR1 (*Discoidin Domain Receptor Tyrosine Kinase 1*) tiene **4 sondas** en el chip:

# - `1007_s_at` (`_s_at`): reconoce múltiples transcritos de DDR1.
# - `207169_x_at`, `208779_x_at`, `210749_x_at` (`_x_at`): pueden tener hibridación cruzada con DDR2 u otros miembros de la familia. Las tres cubren regiones distintas del gen y aportan redundancia para confirmar la señal.

pacman::p_load(dplyr)

ddr1_rows <- gene_info %>%
  filter(SYMBOL == "DDR1")

print(ddr1_rows)

detach("package:dplyr", unload = TRUE)

############################################################
# 7. En gse21779 tenemos 1354896 filas. ¿Por qué?
############################################################

# 1. Redundancia y Validación: Múltiples sondas permiten confirmar la expresión de un gen, aumentando la confiabilidad de los datos.
# 2. Detección de Variantes: Diferentes sondas pueden detectar distintas variantes o isoformas de ARN mensajero del mismo gen.
# 3. Control de Calidad: Tener varias sondas actúa como un control de calidad interno para evaluar la consistencia y precisión de las mediciones.
# 4. Cobertura Completa: Sondas dirigidas a diferentes regiones del gen (como exones alternativos) aseguran una cobertura completa del gen.
# 5. Variabilidad Técnica: Ayuda a mitigar la variabilidad técnica inherente al proceso de hibridación en microarreglos.
# 6. Revisa el documento Arrancando con... affyBatch

# `exprs(gse21779)` sobre un objeto `AffyBatch` sin normalizar devuelve las intensidades a nivel de **sonda individual**, no de probe set. El número surge de:

# - El chip hgu133plus2 tiene **54675 probe sets**, cada uno con aproximadamente **22 sondas** (Perfect Match y MisMatch).
# - 54675 × 22 ≈ **1 202 850** sondas de hibridación.
# - A eso se suman sondas de **control interno** (BioB, BioC, BioD, CreX, AFFX-controls) que suman las restantes hasta alcanzar **1 354 896**.

# Este es el nivel crudo; tras `rma()` o `mas5()` se colapsa a 54675 valores (uno por probe set).

# Confirmamos las dimensiones: sondas individuales × muestras
dim(intensity(gse21779))

# Número de probe sets (tras el resumen)
length(geneNames(gse21779))

# Ratio aproximado: sondas / probe sets
nrow(intensity(gse21779)) / length(geneNames(gse21779))

############################################################
# 8. Como en gse21779 tenemos muchas sondas, añade una cuarta columna (GENETYPE) a gene_info en el que indiques el tipo de cada sonda. Imprime la salida de la función `head` para gene_info. ¿hay algo que te llame la atención?
############################################################

# Usamos `org.Hs.eg.db` para obtener el tipo biológico de cada gen a partir de su ENTREZID. Lo **llamativo** es que tras el `left_join` el dataframe crece de forma explosiva (de ~57 000 a más de 100 millones de filas). Esto ocurre porque `select()` sobre `org.Hs.eg.db` puede devolver **múltiples filas por ENTREZID** (distintas entradas de GENETYPE), y el `left_join` realiza el producto cartesiano de todas las coincidencias. Además, los ~10025 ENTREZID que eran NA siguen siendo NA en GENETYPE.

pacman::p_load(org.Hs.eg.db)

# Obtenemos el tipo de gen para cada ENTREZID presente en gene_info
gene_types <- select(org.Hs.eg.db,
                     keys    = gene_info$ENTREZID,
                     columns = "GENETYPE",
                     keytype = "ENTREZID")

detach("package:org.Hs.eg.db", unload = TRUE)

library(dplyr)

# El left_join puede causar expansión si hay duplicados en gene_types
gene_info <- gene_info %>%
  left_join(gene_types, by = "ENTREZID")

detach("package:dplyr", unload = TRUE)

# Lo que llama la atención: el número de filas explota
dim(gene_info)

# Y la mayoría son NAs en GENETYPE
print(head(gene_info))

############################################################
# 9. Obtén el número de elementos de cada tipo (y los propios tipos, incluyendo los NAs) que hay en la columna GENETYPE del dataframe gene_info. Fíjate que ahora gene_info tiene 100649403 filas, de los cuáles 100500625 son NAs :)
############################################################
table(gene_info$GENETYPE, useNA = "ifany")

#Si filtramos las filas repetidas
library(dplyr)
gene_info_u <- gene_info %>%
  distinct()

table(gene_info_u$GENETYPE, useNA = "ifany")
detach("package:dplyr", unload = TRUE)

############################################################
# 10. Usando dplyr busca las 5 primeras filas que contengan NAs en la columna GENETYPE. Se corresponde con la sonda 1552258_at que es de C2orf59 long intergenic non-protein coding RNA.
############################################################
library(dplyr)

na_gene_info <- gene_info %>%
  filter(is.na(GENETYPE)) %>%
  head(5)

print(na_gene_info)

detach("package:dplyr", unload = TRUE)

############################################################
# 11. Cuando os indiqué que ejecutáis este comando `head(assay(parathyroidGenesSE),n=2)`, en las filas tenemos elementos de nombre ENSG00000000003. ¿Qué es esto?
############################################################

# ENSG00000000003 es un **identificador estable de Ensembl para un gen humano**
# La estructura es:
# - **ENS**: prefijo Ensembl
# - **G**: tipo de feature (Gene; T = Transcript, P = Protein, E = Exon)
# - **00000000003**: número secuencial de 11 dígitos, único para cada gen
# - Prefijo de especie implícito: ENSG = Homo sapiens; ENSMUSG = Mus musculus
#
# En RNA-seq se usan IDs de Ensembl porque son:
# 1. **Unívocos**: cada gen tiene un ID único e irrepetible
# 2. **Estables**: no cambian entre releases (a diferencia del Gene Symbol)
# 3. **Versionados**: ENSG00000000003.15 indica la versión 15 del gen
#
# Ventaja sobre Gene Symbol: Los símbolos pueden tener sinónimos, cambiar
# de nombre, o ser ambiguos. El ID de Ensembl identifica inequívocamente
# la entidad genómica usada en el pipeline de mapeo (ejercicio 15).

pacman::p_load(SummarizedExperiment, parathyroidSE)
data(parathyroidGenesSE, package = "parathyroidSE")

# Mostramos 2 filas: los nombres de fila son IDs Ensembl
head(assay(parathyroidGenesSE), n = 2)

# Los primeros 5 IDs de Ensembl
rownames(assay(parathyroidGenesSE))[1:5]

############################################################
# 12. Es decir, se trabaja con Gene Symbol, EntrezID, sondas affymetrix e indicadores únicos de Ensembl (y más!). Crea un primer diccionario mediante el paquete biomaRt que tenga Gene Symbol, entrezid e identificador único de ensembl a partir de todos los identificadores Ensembl que tenga el paquete.
# PISTA: obtendrás 78610 Ensembl ID, haz de nuevo un table para ver el `gene_biotype`. Por ejemplo verás que salen 24057 protein_coding
############################################################

pacman::p_load(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_data <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "entrezgene_id",   "gene_biotype"),
  mart = ensembl
)

print(head(gene_data))
dim(gene_data)                          # ~91688 filas
table(gene_data$gene_biotype)           # 24080 protein_coding

detach("package:biomaRt", unload = TRUE)

############################################################
# 13. Ahora crea un nuevo objeto a partir de los nombres de gen único de ensembl de parathyroidGenesSE que puedas usar como diccionario y que tenga asociado su external_gene_name, entrezgene_id y gene_biotype para cada ensembl_gene_id. Obtén el tamaño y haz de nuevo el table de la columna gene_biotype. ¿Cuántos id únicos de ensembl son de protein_coding y cuántos de lncRNA?
############################################################

# A diferencia del ejercicio 12 (que descargaba **todos** los genes del 
# genoma humano, ~91718), ahora filtramos `getBM()` para que solo devuelva 
# información de los **63193 genes presentes** en `parathyroidGenesSE`. 
# Usamos `filters = "ensembl_gene_id"` y `values = ensembl_ids` para 
# limitar la consulta a Ensembl.

# **ADVERTENCIA:** Aquí también ocurre **1:many mapping**. Algunos IDs de 
# Ensembl devuelven múltiples filas en `getBM()` por tener múltiples 
# anotaciones (entrezgene_id alternativo, gene_biotype duplicado, etc.). 
# El `merge()` propagará estas duplicaciones, inflando el resultado de 
# 63193 a ~65000+ filas. Para el análisis posterior, deberías quedarte 
# con **una fila por gen** usando `distinct()`.

# Sin dups tenemos: 20333 protein_coding y 11724 lncRNA
# Con dups tenemos: 20869 protein_coding y 11883 lncRNA


pacman::p_load(biomaRt, SummarizedExperiment, dplyr)

ensembl     <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(assay(parathyroidGenesSE)) 

gene_data <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "entrezgene_id",   "gene_biotype"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = ensembl
)

# PROBLEMA: getBM puede devolver múltiples filas para el mismo ensembl_gene_id
# Verificamos si hay duplicados
duplicados <- gene_data %>%
  group_by(ensembl_gene_id) %>%
  filter(n() > 1) %>%
  arrange(ensembl_gene_id)

# Eliminar duplicados manteniendo la primera ocurrencia
gene_data_unique <- gene_data %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Diccionario con todos los IDs originales (all.x = TRUE para no perder ninguno)
result <- data.frame(ENSEMBL_ID = ensembl_ids, stringsAsFactors = FALSE)
result <- merge(result, gene_data_unique,
                by.x = "ENSEMBL_ID", by.y = "ensembl_gene_id",
                all.x = TRUE)
dim(result)
print(head(result))

# Conteo de tipos de genes
table(result$gene_biotype, useNA = "ifany")


# Sin quitar dulpicados
result_con_duplicados <- data.frame(ENSEMBL_ID = ensembl_ids, 
                                     stringsAsFactors = FALSE)
result_con_duplicados <- merge(result_con_duplicados, gene_data,
                                by.x = "ENSEMBL_ID", 
                                by.y = "ensembl_gene_id",
                                all.x = TRUE)

# Conteo de tipos de genes
table(result_con_duplicados$gene_biotype, useNA = "ifany")

detach("package:biomaRt", unload = TRUE)
detach("package:SummarizedExperiment", unload = TRUE)
detach("package:dplyr", unload = TRUE)


############################################################
# 14. ¿Cuántas pacientes de la cohorte TCGA-BRCA tienen subtipo Her2 de acuerdo a la clasificación PAM50? Haz un table con todos los subtipos para ver números globales
############################################################
# PAM50 (Prediction Analysis of Microarray 50) clasifica los tumores de 
# mama en 5 subtipos intrínsecos basándose en la expresión de 50 genes:
# 
# | Subtipo | Caracterización molecular | Pronóstico |
# |---------|--------------------------|------------|
# | LumA    | ER+, PR+, HER2−, Ki67 bajo | Mejor |
# | LumB    | ER+, HER2± Ki67 alto | Intermedio |
# | Her2    | HER2 amplificado | Mejorado con Trastuzumab |
# | Basal   | Triple negativo (ER−, PR−, HER2−) | Peor, sensible a quimio |
# | Normal  | Perfil similar a epitelio normal | Variable |
#
# Her2 es un subtipo caracterizado por amplificación del gen HER2, 
# lo que indica un peor pronóstico pero mejor respuesta a terapia 
# dirigida (Trastuzumab/Herceptin). Representa ~15-20% de canceres mama.

# Tabla de subtipos PAM50
table(meta$paper_BRCA_Subtype_PAM50)

# Con proporciones
round(prop.table(table(meta$paper_BRCA_Subtype_PAM50)) * 100, 1)


############################################################
# 15. Obten los datos de las cinco primeras pacientes para ver sus conteos. ¿Qué son los nombre de columna y por qué tienen la forma ENSG00000000003.15?
############################################################

############################################################
# 16. Usando `dplyr` obten todas las columnas que tengan como nombre un gen que comience por `ENSG00000185960`. ¿Por qué tiene ENSG00000185960.14 y ENSG00000185960.14_PAR_Y?
############################################################

############################################################
# 17. Investiga acerca del objeto `ExpressionSet` del paquete `BioBase`, usado para almacenar matrices de expresión junto con datos fenotípicos y de características. Usos y métodos habituales.
############################################################
# https://www.bioconductor.org/packages/devel/bioc/vignettes/Biobase/inst/doc/ExpressionSetIntroduction.pdf

############################################################
# 18. Investiga acerca del objeto `SummarizedExperiment` objeto contenedor típico de `Bioconductor`, usado para almacenar matrices de conteos, junto con metadatos de fila y columna. Usos y métodos habituales.
############################################################
# https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html

############################################################
# 19. Utilizando `ggplot2` genera un boxplot donde representes la expresión de los genes TP53, PIK3CA, GATA3 y KMT2C, para cada subtipo PAM50, solo para pacientes que tengan `Primary solid Tumor`
############################################################

############################################################
# 20. Utilizando `ggplot2` genera un boxplot donde representes la expresión de los genes TP53 tanto de la copia general del gen, como de la región PAR del mismo, para cada `ajcc_pathologic_t`
############################################################

############################################################
# 21. Utilizando biomaRt, descarga las secuencias de algún gen tanto en homo sapiens como en mus musculus y realiza el alineamiento local y global para ver el grado de conservación del gen entre especies. Investiga qué valores podías obtener con [pwalign](https://bioconductor.org/packages/release/bioc/html/pwalign.html). Puedes apoyarte en el [UCSC Genome Browser](https://genome.ucsc.edu/index.html) para buscar genes compartidos entre ambas.
############################################################
