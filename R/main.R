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

# Définition des paramètres pour l'entrainement
nplis <- 5
k_totest <- c(1, 3, 5, 7, 9, 13, 17, 21)

# Définition de la graine pour la reproductibilité
set.seed(123)

# Lecture des fichiers
AST <- read_tsv("Data/NCBI_AST.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))
AMR <- read_tsv("Data/NCBI_AMR.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))

# Filtrage des données
donnees <- filtre_donnes(AST, AMR)

# Définition des indices des plis des données
idx_plis <- plier(dim(donnees)[1], 5, 0.2)

# ╔═════════════════════════════════════════╗
# ║           Boucle sur les plis           ║
# ╚═════════════════════════════════════════╝

plan(multisession, workers = parallel::detectCores())

# Boucle pour tous les plis en utilisant sapply du paquet future
error_par_pli <- future_sapply(
    
    1:nplis, 
    function(i) {
        
        # Extraction des set de données train et val
        train <- donnees[idx_plis[[i]]$train, ]
        val <- donnees[idx_plis[[i]]$val, ]
        
        # ACM -> KNN
        acm <- acm_custom(train, val)
        
        error_par_k <- sapply(
            k_totest,
            function(k) {
                knn <- knn_custom(acm$entrainement_acm, acm$test_acm, nb_voisins = k)
            
                # Calcul de l'erreur par rapport à la validation
                pred <- as.numeric(knn[["predictions"]])
                error <- sum(pred!=val$Resistance.phenotype) / length(pred)
            }
        )
    },
    future.seed=TRUE
)

error_par_pli <- as.data.frame(error_par_pli)
colnames(error_par_pli) <- paste0("Pli ", c(1:nplis))
rownames(error_par_pli) <- paste0("k=", k_totest)

    
   