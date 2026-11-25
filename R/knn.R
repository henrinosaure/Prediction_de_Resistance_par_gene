# Author : Henry Cloutier
# Date : nov 2026
# 
# 
#    Description
#  ═══════════════
#  fonction customisable pour faire la knn
# 
# 
#    Arguments
#  ═════════════
# • entrainement (matrice)
#     matrice d'entrainement
#     - colonne Resistance.phenotype
#     - colonnes de présence de gènes (ou dimensions dans le cas de l'ACM)
# 
# • validation (matrice)
#     matrice de validation
#     - colonne Resistance.phenotype
#     - colonnes de présence de gènes (ou dimensions dans le cas de l'ACM)
# 
# • distance (string)
#     distance choisie :
#     "manhattan" = Manhattan
#     "euclidean" = Euclidienne
#     "jaccard_coefficient" = Jaccard
#     Utiliser Jaccard uniquement sur les données binaires (sans ACM).
#     Si jaccard est choisi, mettez kernel à NULL
#
# • nb_voisins (int)
#   nombre de voisins
#     
# • kernel (string)
#     noyau de pondération des voisins choisi :
#     "uniform" = Uniforme (poids égaux)
#     "triangular" = Triangulaire
#     "epanechnikov" = Epanechnikov
#     "biweight" = Biquadratique
#     "gaussian" = Gaussien
#     "inverse" = Inverse de la distance
#     "exponential" = Exponentiel
#     NULL si distance est jaccard
# 
#    Returns
#  ═══════════
# une liste contenant 2 choses:
# 1. model : une df avec 2 colonnes : 
#               colonne 1 = %chance susceptible
#               colonne 2 = %chance resistant
# 2. prediction : un vecteur de la prédiction de chaque observation:
#                 0 si susceptible
#                 1 si résistant
#

#install.packages("KernelKnn")
library(KernelKnn)
library(dplyr)

knn_custom <- function(entrainement, validation, distance = "manhattan", nb_voisins = 5, kernel = NULL) {
  
  # =============================================================================
  # ------------------------- préparer pour la knn -----------------------------
  # =============================================================================
  # Make sure that both entrainement and validation dont have BioSample 
  # and also have Resistance.phenotype as a factor
  entrainement <- entrainement %>%
    select(-any_of("BioSample")) %>%
    select(-any_of("X")) %>%
    mutate(Resistance.phenotype = as.factor(Resistance.phenotype))
  validation <- validation %>%
    select(-any_of("BioSample")) %>%
    select(-any_of("X")) %>%
    mutate(Resistance.phenotype = as.factor(Resistance.phenotype))
  
  # =============================================================================
  # -------------------------- librairie KernelKnn ------------------------------
  # =============================================================================
  
  # Cette librairie nécessite que les banques de données soient sans leur classe
  # genes_entrainement s'agit de tous les gênes sans la résistance
  # même affaire pour genes test
  genes_entrainement <- entrainement %>%
    select(-Resistance.phenotype)
  
  genes_test <- validation %>%
    select(-Resistance.phenotype)
  
  # resistance_entrainement c basically la colonne de Resistance.phenotype
  resistance_entrainement <- as.factor(entrainement$Resistance.phenotype)
  
  # model est un
  model <- KernelKnn(
    data = genes_entrainement,
    TEST_data = genes_test,
    y = as.integer(resistance_entrainement),
    k = nb_voisins,
    method = distance,
    weights_function = kernel,
    regression = FALSE,
    Levels = seq_along(levels(resistance_entrainement)),
    threads = max(1, parallel::detectCores() - 2)
    ) %>%
    
    as.data.frame() %>%
    rename(susceptible = class_1, resistant = class_2)
  
  predictions <- levels(resistance_entrainement)[max.col(model)]
  
  return(list(model = model, 
            predictions = predictions))
}
