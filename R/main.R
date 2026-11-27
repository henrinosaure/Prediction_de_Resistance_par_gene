# ╔═════════════════════════════════════════╗
# ║            Paquets & fichiers           ║
# ╚═════════════════════════════════════════╝

source("R/filtre_donnees.R")
source("R/plier.R")
source("R/acm.R")
source("R/knn.R")
source("R/infos_knn.R")

library(FactoMineR)
library(factoextra)
library(dplyr)
library(future)
library(future.apply)
library(caret)

# ╔═════════════════════════════════════════╗
# ║                   Main                  ║
# ╚═════════════════════════════════════════╝

# Définition des paramètres pour l'entrainement
poss_k <- c(1, 3, 5, 7, 9, 11, 13, 15, 17)
poss_acm <- c(TRUE, FALSE)
poss_distance <- c("manhattan", "euclidean", "jaccard_coefficient")
nplis <- 5


# Définition de la graine pour la reproductibilité
set.seed(123)

# Lecture des fichiers
AST <- read_tsv("Data/NCBI_AST.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))
AMR <- read_tsv("Data/NCBI_AMR.tsv", show_col_types = FALSE) %>%
    rename_with(~ sub("^#BioSample$", "BioSample", .x))

# Filtrage des données
donnees_filtrees <- filtre_donnes(AST, AMR)

# Séparation en test / développement (train et val)
idx_dev <- createDataPartition(
    y = donnees_filtrees$Resistance.phenotype,
    p = 0.8,
    list = FALSE
)

donnees_dev <- donnees_filtrees[idx_dev,]
donnees_test <- donnees_filtrees[-idx_dev,] 

# Définition des indices des plis des données
idx_plis <- plier(dim(donnees_dev)[1], 5, 0.2)

# ╔═════════════════════════════════════════╗
# ║                 Boucle                  ║
# ╚═════════════════════════════════════════╝

# Initiation de la session de multithreading
plan(multisession, workers = parallel::detectCores())

# Création d'un grand data frame qui contient toues les combinaisons de k, ACM ou non, de
# calcul de distance et de pli
combinaisons <- expand.grid(
    k = poss_k,
    do_acm = poss_acm,
    distance = poss_distance
)

# --- FILTRAGE DE 2 COMBINAISONS ---
# 1. La distance calculée avec l'indice de jaccard utilise des données binaires. On ne
# peut donc pas utiliser l'ACM avec Jaccard. On doit donc enlever les combinaisons qui
# font une ACM et qui utilisent l'indice de Jaccard.
#
# 2. Dans les combinaisons ou on ne fait pas d'ACM, on travaille avec des distances sur
# des variables binaires. Cela signifie que la seule différence entre les distances de
# Mannhattan et Euclidienne est que la distance Euclidienne contient une racine carrée
# (Man = Euc²). Par conséquent, l'ordre de distance des K voisins est la même pour ces 2
# distances. Il est donc redondant de claculer une KNN pour ces 2 distances. On ne
# calculera que Manhattan (plus simple). On retire donc les combinaisons qui croisent
# Euclidienne et pas d'ACM.

combinaisons <- combinaisons %>% 
    filter(!(do_acm == TRUE & distance == "jaccard_coefficient")) %>%
    filter(!(do_acm == FALSE & distance == "euclidean"))

all_err <- array(
    NA,
    dim = c(length(poss_k), length(poss_acm), length(poss_distance), nplis, 4),
    dimnames = list(
        k = poss_k,
        acm = poss_acm,
        distance = poss_distance#,
        # metrique = c("exactitude", "precision", "sensibilite", "specificite")
    )
)

# Début de la boucle sur toutes les combinaisons
# Utilisation de future.apply pour le multithreading
# La boucle retourne une liste contenant les erreurs de toutes les combinaisons
err_list <- future_lapply(1:dim(combinaisons)[1], function(i) {
        
    # Extraction des paramètres pour l'itération courante
    k <- combinaisons$k[i]
    do_acm <- combinaisons$do_acm[i]
    distance <- as.character(combinaisons$distance[i])
    
    # Manipulations sur chaque pli. La variable metriques contiendra les valeurs pour
    # tous les plis
    metriques <- sapply(1:5, function(pli) {
        
        # Extraction des set de données train et val
        train <- donnees_dev[idx_plis[[pli]]$train, ]
        val <- donnees_dev[idx_plis[[pli]]$val, ]
        
        # ACM si nécessaire
        if (do_acm) {
            donnees_projetees <- acm_custom(train, val)
        } else {
            donnees_projetees <- list(train = train, val = val)
        }
        
        # KNN
        knn <- knn_custom(donnees_projetees$train, donnees_projetees$val, distance, k)
        
        # Calcul de l'erreur par rapport à la validation
        gt <- val$Resistance.phenotype
        pr <- as.numeric(knn[["predictions"]])
        metriques_du_pli <- infos_knn(pr, gt)
        as.numeric(metriques_du_pli[-1]) # On retire la matrice de confusion
    })
    
    # On fait la moyenne pour chaque pli
    metriques <- apply(metriques, 1, mean)
    names(metriques) <- c("exactitude", "precision", "sensibilite", "specificite")
    
    # Valeur de retour de chaque itération de la boucle
    data.frame(
        exactitude = metriques["exactitude"],
        precision = metriques["precision"],
        sensibilite = metriques["sensibilite"],
        specificite = metriques["specificite"],
        k = as.character(k),
        do_acm = as.character(do_acm),
        distance = as.character(distance)
    )
}, future.seed=TRUE)

# On forme un gros data.frame qui contient toutes les métriques
err_list <- bind_rows(err_list)

# Ciblage de la meilleure combinaison
print(filter(err_list, exactitude == max(exactitude)))

write.csv(err_list, "Data/err_list.csv", row.names = TRUE)