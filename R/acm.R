#install.packages("FactoMineR")
#install.packages("factoextra")

library(FactoMineR)
library(factoextra)
library(dplyr)

# INPUT : -entrainement, une matrice contenant les colonnes Resistance.phenotype, 
#         Biosample et tous les genes en valeurs binaires
#         - test, meme format
# OUTPUT: une liste contenant entrainement_acm et test_acm. Ce sont les matrices
#         entrainement et test, projetés sur la matrice de l'ACM. Ils ont la colonne
#         Resistance.phenotype
acm_custom <- function (entrainement, test){
  
  entrainement_original <- entrainement
  test_original <- test
  
  entrainement <- entrainement_original %>%
    select(-any_of(c("Resistance.phenotype", "BioSample", "X"))) %>%
    # enlève les colonnes où toutes les valeurs sont la même (dans ce cas, sont tous égaux à 0)
    select(where(~ n_distinct(.) > 1)) %>%
    mutate(across(everything(), as.factor))
  
  # Fait l'ACM
  acm1 <- MCA(entrainement, ncp = 650, graph = FALSE)
  
  # Lire banque_test
  test <- test_original %>%
    # garder exactement les mêmes gènes que pour l'ACM
    select(all_of(names(entrainement)))
  test_corrige <- test
  
  # Codé par IA : 
  for (col in names(entrainement)) {
    
    # Modalités originales de l'entraînement
    original <- as.character(entrainement[[col]])
    
    # Modalités enregistrées dans l'ACM
    interne <- as.character(acm1$call$X[[col]])
    
    # Construire la correspondance
    indices <- !duplicated(original)
    correspondance <- setNames(interne[indices], original[indices])
    
    # Vérifier les modalités inconnues
    valeurs_test <- as.character(test[[col]])
    inconnues <- setdiff(unique(valeurs_test), names(correspondance))
    
    if (length(inconnues) > 0) {
      stop(paste("Modalités inconnues dans", col,
                 paste(inconnues, collapse = ", ")))
    }
    
    # Appliquer le codage de l'ACM
    test_corrige[[col]] <- factor(
      unname(correspondance[valeurs_test]),
      levels = levels(acm1$call$X[[col]])
    )
  }
  
  # Projeter toutes les observations de test
  projection_test <- predict(acm1, newdata = test_corrige)
  
  # Critère % variance cumulée >1/p
  n_axes <- which(acm1$eig[, 2] <= 100/650)[1]
  
  
  
  # Matrice des coordonnées des observations (gènes seulement), entrainement
  X_entrainement_reduit <- acm1$ind$coord[, 1:n_axes, drop = FALSE]
  
  # Matrice des coordonnées des observations (gènes seulement), test
  X_test_reduit <- projection_test$coord[, 1:n_axes, drop = FALSE]
  
  entrainement_acm <- cbind(entrainement_original[, "Resistance.phenotype", drop = FALSE], X_entrainement_reduit)
    
  test_acm <- cbind(test_original[, "Resistance.phenotype", drop = FALSE], X_test_reduit)
  
  return(list(
    train = entrainement_acm,
    test = test_acm
  ))
}
