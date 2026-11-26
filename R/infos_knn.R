# Author : Henry Cloutier
# Date : nov 2026
# 
# 
#    Description
#  ═══════════════
#  fonction qui trouve la matrice de confusion, sensibilite, specificite, precision et exactitude des resultats
# 
# 
#    Arguments
#  ═════════════
# • predictions_test (vecteur)
#     vecteur, contenant seulement des 0 et des 1, des predictions de resistance/susceptibilite
# 
# • resistance_test (vecteur)
#     vecteur, contenant seulement des 0 et des 1, des vraies valeurs du test
# 
# 
#    Returns
#  ═══════════
# une liste contenant 5 choses:
# 1. matrice_confusion (df) : une df avec 2 colonnes, 2 rangees
#               colonne 1 = predit pos
#               colonne 2 = predit neg
#               rangee 1 = reel pos
#               rangee 2 = reel neg
# 2. sensibilite (int)
# 3. specificite (int)
# 4. precisiion (int)
# 5. precision (int)
#




infos_knn <- function(predictions_test, resistance_test) {
  # vrais positifs
  vp <- sum((predictions_test == 1) & (resistance_test == 1))
  # faux positifs
  fp <- sum((predictions_test == 1) & (resistance_test == 0))
  # vrais negatifs
  vn <- sum((predictions_test == 0) & (resistance_test == 0))
  # faux negatifs
  fn <- sum((predictions_test == 0) & (resistance_test == 1))
  
  matrice_confusion <- data.frame(
    Predit_positif = c(vp, fp),
    Predit_negatif = c(fn, vn),
    row.names = c("Reel positif", "Reel negatif")
  )
  
  sensibilite <- vp/(vp+fn)
  specificite <- vn/(vn + fp)
  precision <- vp/(vp + fp)
  exactitude <- (vp + vn)/(vp + vn + fp + fn)
  
  return(list(matrice_confusion = matrice_confusion, 
              sensibilite = sensibilite, 
              specificite = specificite, 
              precision = precision, 
              exactitude = exactitude))
}
