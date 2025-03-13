# Ejercicios resueltos Boletín 1

# carga de datos 
# disponibles en el siguiente enlace: https://udcgal-my.sharepoint.com/:u:/g/personal/carlos_fernandez_udc_es/EV5v2btIVYdEu2iMzNE3liEBnaRI6jNsx9rgkilFd5Fonw?e=9aT0pB

load("~/datos/GSE21779/gse21779.rda")
meta <- readRDS(file = '~/datos/TCGA-BRCA/meta.RDS')
data.conteos <- readRDS(file = '~/datos/TCGA-BRCA/data.conteos.RDS')

############################################################
#1. Con los datos gse21779 muestra las cinco primera filas y las cinco primeras columnas.
############################################################

exprs(gse21779)[1:5,1:5]

############################################################
#2. ¿Qué son los valores que aparecen en las nombres de columna? Ejemplo GSM542488.CEL.gz (revisa lo que descargaste si tienes dudas)
############################################################

# las muestras

############################################################
#3. ¿Qué es lo que hay en las filas? ¿Cómo puedes obtener el nombre? Saca los 5 primeros. Ejemplo `1007_s_at` ¿Qué tamaño tiene? 54675 
############################################################

featureNames(gse21779)

# revisa el documento Arrancando con... affyBatch

############################################################
# 4. Mapea el nombre de las sondas a genes, por ejemplo usa Gene symbol (nombre oficial del gen) y EntrezID (identificador numérico único asignado a cada gen dentro de la base de datos de Entrez Gene del NCBI). **Ayuda**, busca en Bioconductor el paquete que te coincida con la salida de `annotation(gse21779)`. Imprime por pantalla los 6 primeros elementos de un dataframe (llámalo gene_info) que tenga tres columnas PROBEID, SYMBOL y ENTREZID. Después mira a ver cuántos NAs hay en cada columna
############################################################

annotation(gse21779)

BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)

claves <- keys(hgu133plus2.db,keytype = "PROBEID")

gene_info <- select(hgu133plus2.db,
                    keys = claves,
                    columns = c("SYMBOL","ENTREZID"),
                    keytype = "PROBEID")

print(head(gene_info))
detach("package:hgu133plus2.db", unload = TRUE)

# PISTA: ¿te ha salido este mensaje? `'select()' returned 1:many mapping between keys and columns`. 
# Piensa un poco por qué, pero ahora tienes un gene_info con 57151 filas y se supone que estás mapeando el nombre de las sondas a genes. 
# un mismo `PROBEID` se asocia con múltiples `SYMBOL` o `ENTREZID`.
# Para resolverlo, tendrías que decidir cómo quieres manejar estas relaciones múltiples
# No hay una única opción, depende de las necesidades de cada momento.

# Opciones:
# Filtrar por una sola asignación (la primera)
# Mantener la estructura y analizar manualmente
# Ver lo que sucede con las PROBEID y analizar pasos a realizar
# ...

# Además, tienes 10025 sondas que no tienen nombre del gen ni entrezID. Obten las primeras 15 sondas que no tienen valor en SYMBOL, verás que el segundo es la sonda 1552563_a_at que no se corresponde a ningún gen, sino a C8orf6 (chromosome 8 open reading frame 6).

na_counts <- sapply(gene_info, function(x) sum(is.na(x)))
print(na_counts)

library(dplyr)
na_symbol_info <- gene_info %>%
  filter(is.na(SYMBOL)) %>%
  head(15)

print(na_symbol_info)
detach("package:dplyr", unload = TRUE)

############################################################
# 5. Busca aquellas filas que tengan duplicados. Utiliza dplyr
############################################################

# ¿Qué es un duplicado para ti? En función de eso, usa la columna adecuada
# En clase vimos diferentes opciones

pacman::p_load(dplyr)

df_duplicated <- gene_info %>%
  group_by(SYMBOL) %>%
  add_count(SYMBOL)%>%
  filter(n()>1)

df_duplicated <- gene_info %>%
  group_by(SYMBOL) %>%
  filter(n()>1)%>%
  arrange(SYMBOL)

gene_info %>%
  group_by(PROBEID) %>%
  filter(n()>1)%>%
  arrange(PROBEID)

detach("package:dplyr", unload = TRUE)

############################################################
# 6. En el ejercicio anterior verás que el gen `DDR1` es uno de los muchos, muestra todas las filas de `gene_info` que contengan en la columna `SYMBOL` el gen `DDR1`
# PISTA: hay cuatro sondas (1007_s_at, 207169_x_at, 208779_x_at, 210749_x_at) para el mismo gen (DDR1)
############################################################

pacman::p_load(dplyr)

# Filtrar filas con el símbolo DDR1
ddr1_rows <- gene_info %>%
  filter(SYMBOL == "DDR1")

# Mostrar los resultados
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

############################################################
# 8. Como en gse21779 tenemos muchas sondas, añade una cuarta columna (GENETYPE) a gene_info en el que indiques el tipo de cada sonda. Imprime la salida de la función `head` para gene_info. ¿hay algo que te llame la atención?
############################################################
detach("package:dplyr", unload = TRUE)
gene_info <- AnnotationDbi::select(hgu133plus2.db,
                    keys = claves,
                    columns = c("SYMBOL","ENTREZID","GENETYPE"),
                    keytype = "PROBEID")

head(gene_info)
detach("package:dplyr", unload = TRUE)

############################################################
# 9. Obtén el número de elementos de cada tipo (y los propios tipos, incluyendo los NAs) que hay en la columna GENETYPE del dataframe gene_info. Fíjate que ahora gene_info tiene 100649403 filas, de los cuáles 100500625 son NAs :)
############################################################

table(gene_info$GENETYPE, useNA = "ifany")

############################################################
# 10. Usando dplyr busca las 5 primeras filas que contengan NAs en la columna GENETYPE. Se corresponde con la sonda 1552258_at que es de C2orf59 long intergenic non-protein coding RNA.
############################################################

library(dplyr)

# Filtrar las filas con NA en GENETYPE y mostrar las primeras 5
na_gene_info <- gene_info %>%
  filter(is.na(GENETYPE)) %>%
  head(15)

print(na_gene_info)
detach("package:dplyr", unload = TRUE)

############################################################
# 11. Cuando os indiqué que ejecutáis este comando `head(assay(parathyroidGenesSE),n=2)`, en las filas tenemos elementos de nombre ENSG00000000003. ¿Qué es esto?
############################################################

# ENSG00000000003 es un identificador único de Ensembl para un gen en el sistema de anotación genómica Ensembl. Estos identificadores son parte del estándar para referencias genómicas y se utilizan para seguir la información detallada de genes a través de diversas bases de datos y herramientas bioinformáticas.

############################################################
# 12. Es decir, se trabaja con Gene Symbol, EntrezID, sondas affymetrix e indicadores únicos de Ensembl (y más!). Crea un primer diccionario mediante el paquete biomaRt que tenga Gene Symbol, entrezid e identificador único de ensembl a partir de todos los identificadores Ensembl que tenga el paquete.
# PISTA: obtendrás 78610 Ensembl ID, haz de nuevo un table para ver el `gene_biotype`. Por ejemplo verás que salen 24057 protein_coding
############################################################

pacman::p_load(biomaRt)

# Conectar con Ensembl usando biomaRt
# en caso de querer descargar de humano:
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# también se podría usar usar useEnsembl()

# Obtener la información necesaria
gene_data <- getBM(attributes = c("ensembl_gene_id", 
                                  "external_gene_name", 
                                  "entrezgene_id", 
                                  "gene_biotype"),
                   mart = ensembl)
print(head(gene_data))
detach("package:biomaRt", unload = TRUE)

############################################################
# 13. Ahora crea un nuevo objeto a partir de los nombres de gen único de ensembl de parathyroidGenesSE que puedas usar como diccionario y que tenga asociado su external_gene_name, entrezgene_id y gene_biotype para cada ensembl_gene_id. Obtén el tamaño y haz de nuevo el table de la columna gene_biotype. ¿Cuántos id únicos de ensembl son de protein_coding y cuántos de lncRNA?
############################################################

#20856 protein_coding y 11884 lncRNA
pacman::p_load(biomaRt,SummarizedExperiment,parathyroidSE)
data(parathyroidGenesSE,package="parathyroidSE")

# Conecta con Ensembl mediante biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Supongamos que tienes el objeto parathyroidGenesSE
# Extraer los rownames
ensembl_ids <- rownames(assay(parathyroidGenesSE))

# ojo:
# > tail(ensembl_ids)
# [1] "LRG_93" "LRG_94" "LRG_96" "LRG_97" "LRG_98" "LRG_99"
# Opción:
# Filtrar elementos que comienzan con 'ENSG'
# filtered_ensembl_ids <- ensembl_ids[grep("^ENSG", ensembl_ids)]
# > length(ensembl_ids)
# [1] 63193
# > length(filtered_ensembl_ids)
# [1] 62893

# Obtener la información deseada de Ensembl
gene_data <- getBM(attributes = c("ensembl_gene_id", 
                                  "external_gene_name", 
                                  "entrezgene_id", 
                                  "gene_biotype"),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = ensembl)

# Convertir en data frame con ENSG como primera columna
result <- data.frame(ENSEMBL_ID = ensembl_ids, stringsAsFactors = FALSE)

# Opción con dplyr:
result <- gene_data %>%
  filter(ensembl_gene_id %in% ensembl_ids)

# Opción con R base:
# con R base podemos seleccionar igual que con dplyr filas y columnas
# para ello usamos el operador [,] sabiendo que el filtrado se hace
# [filas,columnas]
result <- gene_data[gene_data$ensembl_gene_id %in% ensembl_ids,]

# Mostrar las primeras filas del resultado
print(head(result))

detach("package:biomaRt", unload = TRUE)
detach("package:SummarizedExperiment", unload = TRUE)
detach("package:parathyroidSE", unload = TRUE)

############################################################
# 14. ¿Cuántas pacientes de la cohorte TCGA-BRCA tienen subtipo Her2 de acuerdo a la clasificación PAM50? Haz un table con todos los subtipos para ver números globales
############################################################

table(meta$paper_BRCA_Subtype_PAM50)
# Basal   Her2   LumA   LumB Normal 
# 197     82    571    209     40 

############################################################
# 15. Obten los datos de las cinco primeras pacientes para ver sus conteos. ¿Qué son los nombre de columna y por qué tienen la forma ENSG00000000003.15?
############################################################

data.conteos[1:5,1:5]

############################################################
# 16. Usando `dplyr` obten todas las columnas que tengan como nombre un gen que comience por `ENSG00000185960`. ¿Por qué tiene ENSG00000185960.14 y ENSG00000185960.14_PAR_Y?
############################################################

pacman::p_load(dplyr)

data.conteos %>%
  select(starts_with("ENSG00000185960"))

detach("package:dplyr", unload = TRUE)

# Esto sugiere que esas columnas representan diferentes aspectos o versiones del mismo gen:
# 
# - **ENSG00000185960.14**: Este es el identificador estándar del gen en Ensembl con una versión específica.
# 
# - **ENSG00000185960.14_PAR_Y**: Esto probablemente indica que se trata de la región pseudoautosómica en el cromosoma Y. Las regiones pseudoautosómicas (PAR) son regiones donde los cromosomas X e Y comparten homología.
# 
# Tener ambas columnas puede significar que estás analizando datos que consideran tanto la copia general del gen como su comportamiento en la región PAR en sistemas de determinación sexual.

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

pacman::p_load(ggplot2, dplyr, tidyr)

genes_symbols <- c("TP53", "PIK3CA", "GATA3", "KMT2C")

gene_symbols_map <- gene_data %>%
  filter(external_gene_name %in% genes_symbols)

colnames(data.conteos) <- sub("\\..*$", "", colnames(data.conteos))

#Esta información está en gene_symbols_map, se hace para mostrar el uso del operador %in%
selected_columns <- colnames(data.conteos)[colnames(data.conteos) %in% gene_symbols_map$ensembl_gene_id]

data_conteos_filtered <- data.conteos %>%
  select(all_of(selected_columns))

colnames(data_conteos_filtered) <- gene_symbols_map$external_gene_name[match(colnames(data_conteos_filtered), gene_symbols_map$ensembl_gene_id)]

meta_filtered <- as.data.frame(meta) %>%
  filter(definition == "Primary solid Tumor")

data_conteos_filtered <- data_conteos_filtered %>%
  filter(rownames(data_conteos_filtered) %in% rownames(meta_filtered))

data_long <- data_conteos_filtered %>%
  mutate(PAM50_Subtype = meta_filtered$paper_BRCA_Subtype_PAM50[match(rownames(.), rownames(meta_filtered))])%>%
  pivot_longer(cols = KMT2C:TP53,names_to = "Gene", values_to = "Expression")

ggplot(data_long, aes(x = PAM50_Subtype, y = Expression, fill = Gene)) +
  geom_boxplot() +
  facet_wrap(~ Gene, scales = "free") +
  theme_minimal() +
  labs(title = "Gene Expression por subtipo PAM50",
       x = "Subtipo PAM50",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

detach("package:ggplot2", unload = TRUE)
detach("package:dplyr", unload = TRUE)
detach("package:tidyr", unload = TRUE)

############################################################
# 20. Utilizando `ggplot2` genera un boxplot donde representes la expresión de los genes TP53 tanto de la copia general del gen, como de la región PAR del mismo, para cada `ajcc_pathologic_t`
############################################################

pacman::p_load(ggplot2, dplyr, tidyr)

tp53_ensembl_id <- gene_data %>%
  filter(external_gene_name == "TP53") %>%
  pull(ensembl_gene_id)

data_conteos_filtered <- data.conteos %>%
  select(starts_with(tp53_ensembl_id))%>%
  mutate(ajcc_pathologic_t = meta_filtered$ajcc_pathologic_t[match(rownames(.), rownames(meta_filtered))])

# solamente hay uno

data_long <- pivot_longer(data_conteos_filtered, 
                          cols = ENSG00000141510,
                          names_to = "transcript",
                          values_to = "expression")

ggplot(data_long, aes(x = ajcc_pathologic_t, y = expression, fill = transcript)) +
  geom_boxplot() +
  labs(title = "Comparación de la expresión de transcritos de TP53",
       x = "Transcrito",
       y = "Expresión") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

detach("package:ggplot2", unload = TRUE)
detach("package:dplyr", unload = TRUE)
detach("package:tidyr", unload = TRUE)

############################################################
# 21. Utilizando biomaRt, descarga las secuencias de algún gen tanto en homo sapiens como en mus musculus y realiza el alineamiento local y global para ver el grado de conservación del gen entre especies. Investiga qué valores podías obtener con [pwalign](https://bioconductor.org/packages/release/bioc/html/pwalign.html). Puedes apoyarte en el [UCSC Genome Browser](https://genome.ucsc.edu/index.html) para buscar genes compartidos entre ambas.
############################################################

pacman::p_load(pwalign, biomaRt)
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Obtén la secuencia del gen TP53 en humanos
human_seq <- getSequence(id = "ENSG00000141510", 
                         type = "ensembl_gene_id", 
                         seqType = "coding", 
                         mart = human_mart)[1, "coding"]

# Configura biomaRt para Mus musculus
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Obtén la secuencia del gen TP53 en ratones
mouse_seq <- getSequence(id = "ENSMUSG00000059552", 
                         type = "ensembl_gene_id", 
                         seqType = "coding", 
                         mart = mouse_mart)[1, "coding"]
# Alineamiento global
global_align <- pairwiseAlignment(human_seq, mouse_seq, type = "global")
print(global_align)

# Alineamiento local
local_align <- pairwiseAlignment(human_seq, mouse_seq, type = "local")
print(local_align)

detach("package:pwalign", unload = TRUE)
detach("package:biomaRt", unload = TRUE)
