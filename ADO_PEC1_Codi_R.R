############# PAC 1. ANÀLISI DE DADES ÒMIQUES. UOC. 6 DE NOVEMBRE 2024. #############
# ---------------------------------------------------------------------------- #

## PRIMER BLOC. SELECCIONAR, CLONAR, I CARREGAR UN DATASET DE METABOLÒMICA:

# Trio el dataset de cachexia (caquèxia en català) del repositori de Github 
# d'A. Sánchez. La caquèxia és una síndrome complexa que s'associa a malalties 
# com el càncer o l'infart de miocardi. Es caracteritza per una pèrdua 
# significativa de pes, de massa muscular, i fatiga. En la cachexia es dona 
# inflamació sistèmica que pot causar canvis perjudicials en el metabolisme i 
# la composició corporal.

# L'adreça del dataset és la següent: 
# https://github.com/nutrimetabolomics/metaboData/tree/main/Datasets/2024-Cachexia

# Adreça del repositori:
url_repositori <- "https://github.com/nutrimetabolomics/metaboData.git"
# Desar el repositori al meu escriptori:
local_path <- "~/Desktop/repository"

# Fent ús de git, fem un clon del repositori:
if (!dir.exists(local_path)) {
  system(paste("git clone", url_repositori, local_path))
} else {
  # Actualitzem el repositori si ja existeix:
  system(paste("cd", local_path, "&& git pull"))
}

# Desem la direcció del fitxer d'interès a la variable "file_path":
file_path <- file.path(local_path, "Datasets/2024-Cachexia/human_cachexia.csv")
# Assignem el dataset d'interès a la variable "metabolomics_data", indiquem 
# que la primera columna s'emprarà com a noms de les fileres:
metabolomics_data <- read.csv(file_path, row.names = 1)

# ------------------------------------------------------------------------------
## SEGON BLOC. CREAR UN CONTENIDOR DEL TIPUS SUMMARIZEDEXPERIMENT:

# Instal·lem els paquets necessaris si encara no els tenim:
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
#}
#BiocManager::install("SummarizedExperiment")

# Carreguem la llibreria "SummarizedExperiment":
library(SummarizedExperiment)
# Mostrem les primeres fileres del dataset:
head(metabolomics_data)
# Carreguem les dades com a matriu a la variable "dades_cachexia":
dades_cachexia <- as.matrix(metabolomics_data)
# Preparem les metadades de les fileres per als metabolits:
row_metadata <- DataFrame(column_id = rownames(dades_cachexia))
# Preparem les metadades de les columnes per a les mostres:
col_metadata <- DataFrame(sample_id = colnames(dades_cachexia))
# Creem el contenidor de SummarizedExperiment:
se <- SummarizedExperiment(assays = list(counts = dades_cachexia), rowData = row_metadata, colData = col_metadata)
# Mostrem el contenidor de SummarizedExperiment:
print(se)

# Ens indica que les dimensions de la matriu són de 77 fileres, i 64 columnes. 
# Ens indica també els noms de les fileres, 77 (78 si comptéssim el "header"), 
# que es corresponen als pacients, i el de les columnes, que es corresponen als 
# metabòlits.


# Obtenim dos arxius que se'ns demanen per la PEC:
# Guardem l'objecte "se" en un fitxer ".Rda":
save(se, file = "~/Desktop/objecte_dades.Rda")
# Extraiem la matriu de dades de "se":
data_matrix <- assay(se)  # Assumes 'counts' is the default assay
# Guardem la matriu de dades en un fitxer .txt
write.table(data_matrix, file = "~/Desktop/cachexia_data.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


# -----------------------------------------------------------------------------
## TERCER BLOC. EXPLORACIÓ DEL DATASET:


# 1. Fem una altra revisió de les metadades del nostre dataset:
# Tornem a explorar part dels noms de les fileres (Patient_ID)
print(rowData(se))
# Tornem a explorar part dels noms de les columnes (Metabòlits):
print(colData(se))
# Obtenim un resum bàsic del dataset:
summary(se)


# 2. Carreguem les llibreries necessàries per l'exploració:
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("corrplot")

# Llibreria per visualització de dades:
library(ggplot2)
# Llibreria que permet passar dades de formats amples a llargs:
library(reshape2)
# Llibreria per a visualitzar matrius de correlació:
library(corrplot)

# 3. Comprovem si hi ha valors no disponibles (missing data):
missing_data <- sum(is.na(assay(se)))
missing_data
# No n'hi ha, tal i com ja s'indica al document "description.md" del repositori.


# Extraiem la matriu de dades principal i l'assignem a assay_data:
assay_data <- assay(se)
# Convertim a numèric si és necessari:
assay_data <- as.data.frame(apply(assay_data, 2, function(x) as.numeric(as.character(x))))


# 4. Calculem la mitjana i la desviació estàndard dels valors dels metabòlits 
# del nostre dataset:

# Primer sense normalitzar les dades:
# Calculem la mitjana i la desviació estàndard dels valors dels metabolits:
mean_values <- rowMeans(assay_data, na.rm = TRUE)
sd_values <- apply(assay_data, 1, sd, na.rm = TRUE)

# Creem un dataframe per a visualitzar les dades:
metabolite_stats <- data.frame(Metabolite = rownames(se), Mean = mean_values, SD = sd_values)
# Visualitzem la mitjana i la desviació estàndard sense normalitzar:
NN <- ggplot(metabolite_stats, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_bar(stat = "identity", fill = "#936025") + 
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Mitjana i SD dels metabolits sense normalitzar", x = "Metabolit", y = "Mitjana")
NN

# A continuació amb les dades normalitzades, fent servir Z-score:
# Normalitzem les dades amb Z-score:
normalized_data <- scale(assay_data, center = TRUE, scale = TRUE)
# Calculem la mitjana i la desviació estàndard dels valors normalitzats
mean_values_normalized <- rowMeans(normalized_data, na.rm = TRUE)
sd_values_normalized <- apply(normalized_data, 1, sd, na.rm = TRUE)
# Crear un dataframe per a la visualització normalitzada
metabolite_stats_normalized <- data.frame(Metabolite = rownames(se), Mean = mean_values_normalized, SD = sd_values_normalized)
# Visualitzar la mitjana i la desviació estàndard amb ggplot2 (normalitzat)
N <- ggplot(metabolite_stats_normalized, aes(x = reorder(Metabolite, Mean), y = Mean)) +
  geom_bar(stat = "identity", fill = "#936025") + 
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Mitjana i SD dels metabolits normalitzats", x = "Metabolit", y = "Mitjana")
N


# 5. Realitzem una exploració de la distribució dels valors dels 
# metabòlits fent servir boxplots:
# Carregar les llibreries necessàries

# Convertim les dades a "long format":
long_data <- melt(as.data.frame(assay_data), variable.name = "Sample", value.name = "Raw_Value")

plot1 <- ggplot(long_data, aes(x = Sample, y = Raw_Value)) +
  geom_boxplot(fill = "#E67E22") +
  theme_minimal() +
  labs(title = "Distribució dels metabolits no normalitzats", x = "Metabolit", y = "Valor dels metabolits no normalitzats") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "lightgrey")
  )
plot1

# 5.2. Boxplot per a cada metabòlit per a veure distribució normalitzada.
# Passem dades normalitzades a "long format":
normalized_long_data <- melt(as.data.frame(normalized_data), variable.name = "Sample", value.name = "Normalized_Value")

plot2 <- ggplot(normalized_long_data, aes(x = Sample, y = Normalized_Value)) +
  geom_boxplot(fill = "#E67E22") +
  theme_minimal() +
  labs(title = "Distribució dels metabolits normalitzats", x = "Metabolit", y = "Valor dels metabolits normalitzats") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "lightgrey")
  )
plot2


# 6. Ara podem realitzar un histograma de distribució global de tots els valors 
# dels metabolits per a veure com es distribueixen.
# 6.1. Primer sense normalitzar:
plot3 <- ggplot(melt(assay_data), aes(x = value)) +
  geom_histogram(bins = 200, fill = "springgreen4", color = "black") +
  theme_minimal() +
  labs(title = "Histograma de la distribució dels valors dels metabolits", x = "Valor del metabolit", y = "Freqüència")
plot3

# 6.2. I a continuació normalitzats:
plot4 <- ggplot(melt(normalized_data), aes(x = value)) +
  geom_histogram(bins = 200, fill = "springgreen2", color = "black") +
  theme_minimal() +
  labs(title = "Histograma de la distribució dels valors dels metabolits normalitzats", x = "Valor del Metabolit", y = "Freqüència")
plot4


# 7. Detecció d'observacions atípiques (outliers)

# Establim un límit per defecte per a la detecció d'outliers (1.5):
outlier_threshold <- 1.5
# Funció per a detectar outliers en base a IQR:
detect_outliers <- function(data) {
  # Primer quartil
  Q1 <- apply(data, 1, quantile, 0.25, na.rm = TRUE)
  # Tercer quartil
  Q3 <- apply(data, 1, quantile, 0.75, na.rm = TRUE)
  # Rang interquartils
  IQR <- Q3 - Q1
  # Matriu lògica per indicar outliers:
  return(data < (Q1 - outlier_threshold * IQR) | data > (Q3 + outlier_threshold * IQR))
}
# Detecció d'outliers en les dades normalitzades:
outliers_normalized <- detect_outliers(normalized_data)
# Mostrarem els outliers fent servir boxplots.
# Convertim les dades normalitzades en format llarg per a poder-hi treballar:
normalized_long_data <- melt(as.data.frame(normalized_data))
# Afegim la informació sobre outliers:
normalized_long_data$outlier <- as.vector(outliers_normalized)
# Generem un boxplot que resalti els outliers amb colors personalitzats:
ggplot(normalized_long_data, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA) +  # Omitir outliers per claredat
  geom_jitter(aes(color = outlier), width = 0.2) +
  theme_minimal() +
  labs(title = "Boxplot de metabolits normalitzats amb outliers marcats",
       x = "Metabolite",
       y = "Valor normalitzat",
       color = "Outlier") +
  scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "grey")) + 
  theme(plot.title = element_text(hjust = 0.5))


# També podem detectar outliers fent servit un "heatmap":
# install.packages("heatmaply")
# Llibreria que ens permet generar heatmaps interactius.
library(heatmaply)
# Creem un heatmap en el que ens destaquen els outliers:
heatmaply(normalized_data, 
          Rowv = NA, Colv = NA,
          colors = ifelse(outliers_normalized, "red", "blue"),
          xlab = "Metabolits", ylab = "Pacients", 
          main = "Heatmap de valors de metabolits, amb els outliers destacats")


## ANÀLISIS ESTADÍSTIQUES

# ------------------------------------------------------------------------------

# 1. Fem una anàlisi de PCA sense els outliers:

# Llibreria per a visualitzar els resultats d'anàlisis multivariants de dades:
library(factoextra)
# Realitzem PCA sense els outliers:
pca_result <- prcomp(data_without_outliers, center = TRUE, scale. = TRUE)
# Obtenim els valors dels metabòlits per a cada component principal 
# (només PC1 i PC2):
PCA <- as.data.frame(pca_result$rotation[, 1:2])
# Afegim els noms dels metabòlits;
PCA$Metabolite <- rownames(PCA)
# Creem el gràfic de PCA on cada punt representa un metabòlit:
plot_pca_metabolits <- ggplot(PCA, aes(x = PC1, y = PC2, label = Metabolite)) +
  geom_point(color = "#0073C2", size = 3, alpha = 0.7) +
  # Etiquetes amb el nom dels metabòlits
  geom_text(vjust = -1, size = 3) +                       
  theme_minimal() +
  labs(title = "Anàlisi de PCA dels metabòlits",
       x = "PC1", y = "PC2") +
  theme(plot.title = element_text(hjust = 0.5))
# Mostrem el gràfic
print(plot_pca_metabolits)


# 2. Després, podem comprovar les correlacions entre metabòlits:
# Calculem una matriu de correlació:
correlation_matrix <- cor(data_without_outliers, use = "pairwise.complete.obs")
# Visualitzem la matriu de correlació amb corrplot:
corrplot(correlation_matrix, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         title = "Matriu de correlació dels metabolits")


# 3. I finalment fem un primer test estadístic (t-student) per fer-nos una idea
# de quins metabòlits poden ser significativament diferents entre persones controls
# i afectades per la síndrome.

# Filtrem les dades per eliminar outliers i preparar la matriu:
assay_data_no_outliers <- as.data.frame(data_without_outliers)
# Afegim la columna 'Muscle.loss' que havíem eliminat a les dades sense outliers:
assay_data_no_outliers$Muscle.loss <- metabolomics_data$Muscle.loss[!outlier_rows]
# Inicialitzem una llista per emmagatzemar els resultats:
results <- list()
# I iterem sobre cada metabòlit (excloent 'Muscle.loss'):
for (metabolite in colnames(assay_data_no_outliers)[-ncol(assay_data_no_outliers)]) {
  # Separem les dades en els dos grups (cachexic i control):
  group1 <- assay_data_no_outliers[assay_data_no_outliers$Muscle.loss == "cachexic", metabolite]
  group2 <- assay_data_no_outliers[assay_data_no_outliers$Muscle.loss == "control", metabolite]
  # Realitzem el test t-Student:
  t_test_result <- t.test(group1, group2, var.equal = TRUE)
  # Emmagatzemar resultats:
  results[[metabolite]] <- c(t_test_result$p.value, t_test_result$estimate, t_test_result$conf.int)
}
# Convertim els resultats a un dataframe:
results_df <- do.call(rbind, results)
colnames(results_df) <- c("p_valor", "mitjana_grup1", "mitjana_grup2", "conf_int_inf", "conf_int_sup")
results_df <- as.data.frame(results_df)

# Filtrem els resultats significatius (p_value < 0.05):
significant_results_df <- subset(results_df, p_valor < 0.05)

# Mostrem la taula completa i la dels resultats significatius:
print("Tots els resultats:")
print(results_df)
print("Resultats significatius (p-value < 0.05):")
print(significant_results_df)
