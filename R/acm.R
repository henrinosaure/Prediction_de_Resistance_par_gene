#install.packages("FactoMineR")
#install.packages("factoextra")

library(FactoMineR)
library(factoextra)
library(dplyr)

# INPUT : -entrainement, une matrice contenant les colonnes Resistance.phenotype, 
#         Biosample et tous les genes en valeurs binaires
#         - val, meme format
# OUTPUT: une liste contenant entrainement_acm et val_acm. Ce sont les matrices
#         entrainement et val, projetés sur la matrice de l'ACM. Ils ont la colonne
#         Resistance.phenotype
acm_custom <- function (entrainement, val){
  
  entrainement_original <- entrainement
  val_original <- val
  
  entrainement <- entrainement_original %>%
    select(-any_of(c("Resistance.phenotype", "BioSample", "X"))) %>%
    # enlève les colonnes où toutes les valeurs sont la même (dans ce cas, sont tous égaux à 0)
    select(where(~ n_distinct(.) > 1)) %>%
    mutate(across(everything(), as.factor))
  
  # Fait l'ACM
  acm1 <- MCA(entrainement, ncp = 650, graph = FALSE)
  
  # Lire banque_val
  val <- val_original %>%
    # garder exactement les mêmes gènes que pour l'ACM
    select(all_of(names(entrainement)))
  val_corrige <- val
  
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
    valeurs_val <- as.character(val[[col]])
    inconnues <- setdiff(unique(valeurs_val), names(correspondance))
    
    if (length(inconnues) > 0) {
      stop(paste("Modalités inconnues dans", col,
                 paste(inconnues, collapse = ", ")))
    }
    
    # Appliquer le codage de l'ACM
    val_corrige[[col]] <- factor(
      unname(correspondance[valeurs_val]),
      levels = levels(acm1$call$X[[col]])
    )
  }
  
  # Projeter toutes les observations de val
  projection_val <- predict(acm1, newdata = val_corrige)
  
  # Critère % variance cumulée >90%
  n_axes <- which(acm1$eig[, 3] >= 90)[1]
  
  
  
  # Matrice des coordonnées des observations (gènes seulement), entrainement
  X_entrainement_reduit <- acm1$ind$coord[, 1:n_axes, drop = FALSE]
  
  # Matrice des coordonnées des observations (gènes seulement), val
  X_val_reduit <- projection_val$coord[, 1:n_axes, drop = FALSE]
  
  entrainement_acm <- cbind(entrainement_original[, "Resistance.phenotype", drop = FALSE], X_entrainement_reduit)
    
  val_acm <- cbind(val_original[, "Resistance.phenotype", drop = FALSE], X_val_reduit)
  
  return(list(
    train = entrainement_acm,
    val = val_acm
  ))
}
