# ╔═════════════════════════════════════════╗
# ║            Paquets & fichiers           ║
# ╚═════════════════════════════════════════╝

source("R/filtre_donnees.R")
source("R/plier.R")
source("R/acm.R")
source("R/knn.R")

library(FactoMineR)
library(factoextra)
library(dplyr)
library(future)
library(future.apply)

# ╔═════════════════════════════════════════╗
# ║                   Main                  ║
# ╚═════════════════════════════════════════╝

# Définition de la graine pour la reproductibilité
set.seed(123)

# Définition du nombre de plis pour l'entrainement
nplis <- 5

# Lecture des fichiers
AST <- read_tsv("Data/NCBI_AST.tsv") %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))
AMR <- read_tsv("Data/NCBI_AMR.tsv") %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))

# Filtrage des données
donnees <- filtre_donnes(AST, AMR)

# Définition des indices des plis des données
idx_plis <- plier(dim(donnees)[1], 5, 0.2)

# ╔═════════════════════════════════════════╗
# ║           Boucle sur les plis           ║
# ╚═════════════════════════════════════════╝

plan(multisession, workers = parallel::detectCores() - 1)

# Boucle pour tous les plis en utilisant sapply du paquet future
variable_a_nommer_comme_du_monde <- future_sapply(
    
    1:nplis, 
    function(i) {
        
        train <- donnees[idx_plis[[i]]$train, ]
        val <- donnees[idx_plis[[i]]$val, ]
        
        acm <- acm_custom(train, val)
        
        knn <- knn_custom(acm$entrainement_acm, acm$test_acm)
        
    }
)
    
   