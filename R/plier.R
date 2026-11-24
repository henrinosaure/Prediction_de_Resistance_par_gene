#
# Author : E. Fortin
# Date : month 20--
# 
# 
#    Description
#  ═══════════════
#  - Cette fonction divise un jeu de données en plis -
# 
# La fonction retourne les indices des observations d'entrainement et validation
# déterminés aléatoirement pour chaque pli. Elle n'utilise pas le jeu de données originel,
# mais seulement sa longueur et répartit les indices en fonction de cette longueur.
# 
# 
# 
#    Arguments
#  ═════════════
# • p_len_data (integer)
#     Le nombre d'observations dans le jeu de données
# 
# • p_nplis (integer)
#     Le nombre de plis souhaités
#
# • p_val_ratio (numeric)
#     Le ratio des observations qui seront placés en validation. Par exemple, si cette
#     valeur est égale à 0.2, 20% des indices seront placés en validation.
# 
# 
#    Retourne
#  ════════════
# La fonction retourne une liste contenant chaque pli. Les plis sont stockés sous forme de
# liste contenant les éléments `val` et `train` contenant respectivement les indices de
# validation et d'entrainement.


# ╔═════════════════════════════════════════╗
# ║                  Main                   ║
# ╚═════════════════════════════════════════╝

plier <- function(p_len_data, p_nplis, p_val_ratio) {
    
    # Vecteur contenant tous les indices dans un ordre aléatoire
    idx = sample(1:p_len_data)
    
    # On coupe le vecteur aléatoire idx en plis (train et val ne sont pas encore séparés)
    groupes <- split(idx, cut(seq_along(idx), p_nplis, labels = FALSE))
    
    # Liste qui contiendra les plis finaux
    plis = vector("list", p_nplis)
    
    # Pour chaque pli, on prend les premiers indices (selon p_val_ratio) qui iront dans la
    # validation, et on met le reste dans train
    for(i in 1:p_nplis) {
        
        groupe <- groupes[[i]] # groupe courant
        nb_obs <- length(groupe) # nombre d'observations (d'indices) dans le pli
        nb_val <- round(nb_obs * p_val_ratio) # nombre d'observations qui iront en val
        
        plis[[i]]$val <- groupe[1 : nb_val] # On place les premiers nb_val indices en val
        plis[[i]]$train <- groupe[(nb_val+1) : nb_obs] # On met le reste en train
    }
    
    plis
}