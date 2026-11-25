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
#     "minkowski" = Minkowski (ordre défini par p)
#     "jaccard_coefficient" = Jaccard
#     Utiliser Jaccard uniquement sur les données binaires (sans ACM).
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
#
# • ordre_minkowski(int)
#     pour customiser l'ordre de minkowski si on veut aller à p=3, 4, 5, etc.
#     seulement utilisé si la distance est minkowski, c'est donc pas nécessaire
#     de l'appeller si on utilise une autre distance
# 
#    Returns
#  ═══════════
# bla_bla
#

#install.packages("KernelKnn")
library(KernelKnn)


knn_custom <- function(entrainement, validation, distance = "manhattan", nb_voisins = 5, kernel = "rectangular", ordre_minkowski=3) {
  # =============================================================================
  # ------------------------- préparer pour la knn -----------------------------
  # =============================================================================
  # Make sure that both entrainement and validation dont have BioSample 
  # and also have Resistance.phenotype as a factor
  entrainement <- entrainement %>%
    select(-any_of("BioSample")) %>%
    mutate(Resistance.phenotype = as.factor(Resistance.phenotype))
  validation <- validation %>%
    select(-any_of("BioSample")) %>%
    mutate(Resistance.phenotype = as.factor(Resistance.phenotype))
  
  # =============================================================================
  # ---------------------- librairie KernelKnn (Jaccard) ------------------------
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
    
    model <- KernelKnn(
      data = genes_entrainement,
      TEST_data = genes_test,
      y = as.integer(resistance_entrainement),
      k = nb_voisins,
      method = distance,
      p = ordre_minkowski,
      weights_function = kernel,
      regression = FALSE,
      Levels = seq_along(levels(resistance_entrainement))
    )
    
    predictions <- levels(resistance_entrainement)[max.col(model)]
  return(predictions)
}

knn_custom(entrainementtt, testtt, distance = "manhattan", nb_voisins = 20, kernel = "rectangular", ordre_minkowski=3)