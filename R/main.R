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
poss_k <- c(1, 3, 5, 7, 9, 13, 17, 21)
poss_acm <- c(TRUE, FALSE)
poss_distance <- c("manhattan", "euclidean", "jaccard_coefficient")
nplis <- 5


# Définition de la graine pour la reproductibilité
set.seed(122)

# Lecture des fichiers
AST <- read_tsv("Data/NCBI_AST.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))
AMR <- read_tsv("Data/NCBI_AMR.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))

# Filtrage des données
donnees_filtrees <- filtre_donnes(AST, AMR)

# TODO garder un test

# Définition des indices des plis des données
idx_plis <- plier(dim(donnees_filtrees)[1], 5, 0.2)

# ╔═════════════════════════════════════════╗
# ║           Boucle sur les plis           ║
# ╚═════════════════════════════════════════╝

# Initiation de la session de multithreading
plan(multisession, workers = parallel::detectCores())

# Création d'un grand data frame qui contient toues les combinaisons de k, ACM ou non, de
# calcul de distance et de pli
combinaisons <- expand.grid(
    k = poss_k,
    do_acm = poss_acm,
    distance = poss_distance,
    pli = c(1:nplis)
)

# La distance calculée avec l'indice de jaccard utilise des données binaires. On ne peut
# donc pas utiliser l'AXM avec Jaccard. On doit donc enlever les combinaisons qui font une
# ACM et qui utilisent l'indice de Jaccard.
combinaisons <- combinaisons %>% 
    filter(!(do_acm == TRUE & distance == "jaccard_coefficient"))

variable_a_nommer_comme_du_monde <- array(
    NA,
    dim = c(length(poss_k), length(poss_acm), length(poss_distance), nplis),
    dimnames = list(
        k = poss_k,
        acm = poss_acm,
        distance = poss_distance,
        pli = paste0("pli", 1:nplis)
    )
)

# Début de la boucle sur toutes les combinaisons
error_list <- future_lapply(
    
    1:dim(combinaisons)[1],
    function(i) {
        
        # Extraction des paramètres pour l'itération courante
        k <- combinaisons$k[i]
        do_acm <- combinaisons$do_acm[i]
        distance <- as.character(combinaisons$distance[i])
        pli <- combinaisons$pli[i]
        
        # Extraction des set de données train et val
        train <- donnees_filtrees[idx_plis[[pli]]$train, ]
        val <- donnees_filtrees[idx_plis[[pli]]$val, ]
        
        # ACM si nécessaire
        if (do_acm) {
            donnees_projetees <- acm_custom(train, val)
        } else {
            donnees_projetees <- list(train = train, test = val)
        }
        
        # KNN
        knn <- knn_custom(donnees_projetees$train, donnees_projetees$test, distance, k)
        
        # Calcul de l'erreur par rapport à la validation
        pred <- as.numeric(knn[["predictions"]])
        error <- sum(pred!=val$Resistance.phenotype) / length(pred)
        
        list(
            error = error,
            k = as.character(k),
            do_acm = as.character(do_acm),
            distance = as.character(distance),
            pli = paste0("pli", pli)
        )
    },
    future.seed=TRUE
)

for (i in 1:length(error_list)) {
    variable_a_nommer_comme_du_monde[
        error_list[[i]]$k, 
        error_list[[i]]$do_acm, 
        error_list[[i]]$distance, 
        error_list[[i]]$pli
    ] <- error_list[[i]]$error
}

# DEBUG EF
view (variable_a_nommer_comme_du_monde)

print(apply(variable_a_nommer_comme_du_monde, c(1,2,3), mean))

    
   