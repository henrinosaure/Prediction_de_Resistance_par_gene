# Main

set.seed(123)

source("R/filtre_donnees.R")
source("R/plier.R")
source("R/acm.R")
# source("R/knn.R")

library(FactoMineR)
library(factoextra)
library(dplyr)

# Définition du nombre de plis pour l'entrainement
nplis <- 5

# Lecture des fichiers
AST <- read_tsv("Data/NCBI_AST.tsv") %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))

AMR <- read_tsv("Data/NCBI_AMR.tsv") %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))

# Filtrage des données
donnees <- filtre_donnes(AST, AMR)

# Plier les données
idx_plis <- plier(dim(donnees)[1], 5, 0.2)

for (i in 1:nplis)
    
    train <- donnees[idx_plis[[i]]$train, ]
    val <- donnees[idx_plis[[i]]$val, ]
    
    acm <- acm_custom(train, val)
    
    # faire un knn, mais tester chaque voisin possible (cehck si deja implementé)