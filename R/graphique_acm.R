# Author : Henry Cloutier
# Date : nov 2026
# 
# 
#    Description
#  ═══════════════
#  fonction qui génère un graphique de la variance cumulée pour l'ACM en fonction du nombre de variables
# 
# 
#    Arguments
#  ═════════════
# • entrainement (matrice)
#     matrice contenant les valeurs d'entrainement
# 
# 
# 
#    Returns
#  ═══════════
#

#install.packages("ggplot2") # Une seule fois

library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggplot2)

# Prépare la banque d'entrainement pour l'ACM
entrainement = read.csv("Data/banque_entrainement.csv") 

entrainement <- entrainement %>%
  select(-any_of(c("Resistance.phenotype", "BioSample", "X"))) %>%
  # enlève les colonnes où toutes les valeurs sont la même (dans ce cas, sont tous égaux à 0)
  select(where(~ n_distinct(.) > 1)) %>%
  mutate(across(everything(), as.factor))

# Fait l'ACM
acm1 <- MCA(entrainement, ncp = 650, graph = FALSE)

# Valeurs pour le graphique : variance cumulée et indexes

dernier_index_choisi = which(acm1$eig[, 2] <= 100/650)[1]
variance_choisie <- acm1$eig[dernier_index_choisi, 3]

df <- data.frame(
  variance_cumulee = variance_cumulee,
  axes = seq_along(variance_cumulee)
)

graphique <- ggplot(data = df, aes(x = axes, y = variance_cumulee)) +
  
  # Courbe principale
  geom_line(color = "#455f51", linewidth = 1.5) +
  
  # Ligne pointillée horizontale
  annotate(
    "segment",
    x = 0, y = variance_choisie,
    xend = dernier_index_choisi, yend = variance_choisie,
    linetype = "dashed",
    color = "#2f3e36",
    linewidth = 0.9
  ) +
  
  # Ligne pointillée verticale
  annotate(
    "segment",
    x = dernier_index_choisi, y = variance_choisie,
    xend = dernier_index_choisi, yend = 0,
    linetype = "dashed",
    color = "#2f3e36",
    linewidth = 0.9
  ) +
  
  # Mise en évidence du nombre d'axes retenus
  annotate(
    "label",
    x = dernier_index_choisi,
    y = 0,
    label = dernier_index_choisi,
    fill = "#455f51",
    color = "white",
    fontface = "bold",
    size = 4.5,
    linewidth = 0
  ) +
  
  # Mise en évidence de l'inertie cumulée
  annotate(
    "label",
    x = 0,
    y = variance_choisie,
    label = paste0(round(variance_choisie, 1), "%"),
    fill = "#455f51",
    color = "white",
    fontface = "bold",
    size = 4.5,
    linewidth = 0
  ) +
  
  # Titres
  labs(
    title = "Sélection du nombre d'axes de l'ACM",
    x = "Nombre d'axes retenus",
    y = "Inertie cumulée expliquée (%)"
  ) +
  
  # Graduations de l'axe X
  scale_x_continuous(
    breaks = seq(0, max(df$axes), by = 100),
    expand = expansion(mult = c(0.05, 0.02))
  ) +
  
  # Graduations de l'axe Y
  scale_y_continuous(
    breaks = seq(0, 100, by = 10),
    limits = c(0, 100),
    labels = function(x) paste0(x, "%")
  ) +
  
  # Apparence
  theme_minimal(base_size = 13) +
  theme(
    plot.background = element_rect(
      fill = "#d1dfce", color = NA
    ),
    panel.background = element_rect(
      fill = "#d1dfce", color = NA
    ),
    panel.grid.major = element_line(
      color = "white", linewidth = 0.7
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(
      face = "bold", color = "#2f3e36"
    ),
    axis.text = element_text(
      color = "#2f3e36"
    ),
    plot.title = element_text(
      face = "bold",
      size = 18,
      color = "#2f3e36",
      hjust = 0.5
    ),
    plot.margin = margin(15, 15, 15, 15)
  )

# Affichage
print(graphique)

# Exportation
ggsave(
  "Graphiques/graphique_ACM.svg",
  plot = graphique,
  width = 25,
  height = 15,
  units = "cm",
  bg = "#d1dfce"
)