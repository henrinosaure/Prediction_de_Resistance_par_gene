
# ╔══════════════════════════════════════════════════════════════════════════════════════╗
# ║                                  Paquets & fichiers                                  ║
# ╚══════════════════════════════════════════════════════════════════════════════════════╝

source("R/filtre_donnees.R")
source("R/plier.R")
source("R/acm.R")
source("R/knn.R")
source("R/infos_knn.R")
source("R/train_test.R")

library(FactoMineR)
library(factoextra)
library(dplyr)
library(future)
library(future.apply)
library(caret)

# ╔══════════════════════════════════════════════════════════════════════════════════════╗
# ║                                         Main                                         ║
# ╚══════════════════════════════════════════════════════════════════════════════════════╝

# Définition de la graine pour la reproductibilité
set.seed(123)

# Lecture des fichiers
AST <- read_tsv("Data/NCBI_AST.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))
AMR <- read_tsv("Data/NCBI_AMR.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))

# Filtrage des données
donnees_filtrees <- filtre_donnes(AST, AMR)

# Entrainement du modèle
modele_final <- train_test(
    donnees_filtrees,
    5,
    c(1, 3, 5, 7, 9, 11, 13, 15, 17),
    c("manhattan", "euclidean", "jaccard_coefficient"),
    TRUE
)

# Envoi à la console des performances du modèle
cat("\nLe modèle final et performances sur le jeu de validation:\n")
print(modele_final$parametres_modele)
cat("\nperformances sur le jeu de test:\n")
print(modele_final$performance)
