# Author : Henry Cloutier
# Date : nov 2026
# 
# 
#    Description
#  ═══════════════
#  fonction qui trouve la matrice de confusion, sensibilite, specificite, precision et
#  exactitude des resultats
# 
# 
#    Arguments
#  ═════════════
# • predictions_test (vecteur)
#     vecteur, contenant seulement des 0 et des 1, des predictions de
#     resistance/susceptibilite
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
# 2. sensibilite (float)
# 3. specificite (float)
# 4. precision (float)
# 5. exactitude (float)
#

infos_knn <- function(pr, gt) {
  
  # vrais positifs
  vp <- sum((pr == 1) & (gt == 1))
  # faux positifs
  fp <- sum((pr == 1) & (gt == 0))
  # vrais negatifs
  vn <- sum((pr == 0) & (gt == 0))
  # faux negatifs
  fn <- sum((pr == 0) & (gt == 1))
  
  matrice_confusion <- data.frame(
    Predit_positif = c(vp, fp),
    Predit_negatif = c(fn, vn),
    row.names = c("Reel positif", "Reel negatif")
  )
  
  exactitude <- (vp + vn) / (vp + vn + fp + fn) # Combien de prédictions sont bonnes?
  precision <- vp / (vp + fp) # Parmi les prédictions positives, combien sont vraies?
  sensibilite <- vp / (vp + fn) # Combien de positifs ont été détéctés?
  specificite <- vn / (vn + fp) # Combien de négatifs ont été détéctés
  
  return(list(
    matrice_confusion = matrice_confusion, 
    exactitude = exactitude,
    precision = precision, 
    sensibilite = sensibilite, 
    specificite = specificite 
    )
  )
}
