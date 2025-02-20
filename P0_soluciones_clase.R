# Ejercicios resueltos Boletín 1

# carga de datos 
# disponibles en el siguiente enlace: https://udcgal-my.sharepoint.com/:u:/g/personal/carlos_fernandez_udc_es/EV5v2btIVYdEu2iMzNE3liEBnaRI6jNsx9rgkilFd5Fonw?e=9aT0pB

load("~/datos/gse21779.rda")
meta <- readRDS(file = '~/datos/meta.RDS')
data.conteos <- readRDS(file = '~/datos/data.conteos.RDS')

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

# las sondas del gen

############################################################
# 4. Mapea el nombre de las sondas a genes, por ejemplo usa Gene symbol (nombre oficial del gen) y EntrezID (identificador numérico único asignado a cada gen dentro de la base de datos de Entrez Gene del NCBI). **Ayuda**, busca en Bioconductor el paquete que te coincida con la salida de `annotation(gse21779)`. Imprime por pantalla los 6 primeros elementos de un dataframe (llámalo gene_info) que tenga tres columnas PROBEID, SYMBOL y ENTREZID. Después mira a ver cuántos NAs hay en cada columna
############################################################

# PISTA: ¿te ha salido este mensaje? `'select()' returned 1:many mapping between keys and columns`. 
# Piensa un poco por qué, pero ahora tienes un gene_info con 57151 filas y se supone que estás mapeando el nombre de las sondas a genes. 
# Además, tienes 10025 sondas que no tienen nombre del gen ni entrezID. Obten las primeras 15 sondas que no tienen valor en SYMBOL, verás que el segundo es la sonda 1552563_a_at que no se corresponde a ningún gen, sino a C8orf6 (chromosome 8 open reading frame 6).

annotation(gse21779)

############################################################
# 5. Busca aquellas filas que tengan duplicados. Utiliza dplyr
############################################################

############################################################
# 6. En el ejercicio anterior verás que el gen `DDR1` es uno de los muchos, muestra todas las filas de `gene_info` que contengan en la columna `SYMBOL` el gen `DDR1`
# PISTA: hay cuatro sondas (1007_s_at, 207169_x_at, 208779_x_at, 210749_x_at) para el mismo gen (DDR1)
############################################################

############################################################
# 7. En gse21779 tenemos 1354896 filas. ¿Por qué?
############################################################

############################################################
# 8. Como en gse21779 tenemos muchas sondas, añade una cuarta columna (GENETYPE) a gene_info en el que indiques el tipo de cada sonda. Imprime la salida de la función `head` para gene_info. ¿hay algo que te llame la atención?
############################################################

############################################################
# 9. Obtén el número de elementos de cada tipo (y los propios tipos, incluyendo los NAs) que hay en la columna GENETYPE del dataframe gene_info. Fíjate que ahora gene_info tiene 100649403 filas, de los cuáles 100500625 son NAs :)
############################################################

############################################################
# 10. Usando dplyr busca las 5 primeras filas que contengan NAs en la columna GENETYPE. Se corresponde con la sonda 1552258_at que es de C2orf59 long intergenic non-protein coding RNA.
############################################################

############################################################
# 11. Cuando os indiqué que ejecutáis este comando `head(assay(parathyroidGenesSE),n=2)`, en las filas tenemos elementos de nombre ENSG00000000003. ¿Qué es esto?
############################################################

############################################################
# 12. Es decir, se trabaja con Gene Symbol, EntrezID, sondas affymetrix e indicadores únicos de Ensembl (y más!). Crea un primer diccionario mediante el paquete biomaRt que tenga Gene Symbol, entrezid e identificador único de ensembl a partir de todos los identificadores Ensembl que tenga el paquete.
# PISTA: obtendrás 91718 Ensembl ID, haz de nuevo un table para ver el `gene_biotype`. Por ejemplo verás que salen 24057 protein_coding
############################################################

############################################################
# 13. Ahora crea un nuevo objeto a partir de los nombres de gen único de ensembl de parathyroidGenesSE que puedas usar como diccionario y que tenga asociado su external_gene_name, entrezgene_id y gene_biotype para cada ensembl_gene_id. Obtén el tamaño y haz de nuevo el table de la columna gene_biotype. ¿Cuántos id únicos de ensembl son de protein_coding y cuántos de lncRNA?
############################################################

############################################################
# 14. ¿Cuántas pacientes de la cohorte TCGA-BRCA tienen subtipo Her2 de acuerdo a la clasificación PAM50? Haz un table con todos los subtipos para ver números globales
############################################################

table(meta$paper_BRCA_Subtype_PAM50)

# 82 pacientes

############################################################
# 15. Obten los datos de las cinco primeras pacientes para ver sus conteos. ¿Qué son los nombre de columna y por qué tienen la forma ENSG00000000003.15?
############################################################

############################################################
# 16. Usando `dplyr` obten todas las columnas que tengan como nombre un gen que comience por `ENSG00000185960`. ¿Por qué tiene ENSG00000185960.14 y ENSG00000185960.14_PAR_Y?
############################################################

############################################################
# 17. Investiga acerca del objeto `ExpressionSet` del paquete `BioBase`, usado para almacenar matrices de expresión junto con datos fenotípicos y de características. Usos y métodos habituales.
############################################################

############################################################
# 18. Investiga acerca del objeto `SummarizedExperiment` objeto contenedor típico de `Bioconductor`, usado para almacenar matrices de conteos, junto con metadatos de fila y columna. Usos y métodos habituales.
############################################################

############################################################
# 19. Utilizando `ggplot2` genera un boxplot donde representes la expresión de los genes TP53, PIK3CA, GATA3 y KMT2C, para cada subtipo PAM50, solo para pacientes que tengan `Primary solid Tumor`
############################################################

############################################################
# 20. Utilizando `ggplot2` genera un boxplot donde representes la expresión de los genes TP53 tanto de la copia general del gen, como de la región PAR del mismo, para cada `ajcc_pathologic_t`
############################################################

############################################################
# 21. Utilizando biomaRt, descarga las secuencias de algún gen tanto en homo sapiens como en mus musculus y realiza el alineamiento local y global para ver el grado de conservación del gen entre especies. Investiga qué valores podías obtener con [pwalign](https://bioconductor.org/packages/release/bioc/html/pwalign.html). Puedes apoyarte en el [UCSC Genome Browser](https://genome.ucsc.edu/index.html) para buscar genes compartidos entre ambas.
############################################################
