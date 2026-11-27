# Graphiques de choix de k

library(ggplot2)
library(dplyr)

err_list <- read.csv("Data/err_list.csv") %>%
  filter(distance == "jaccard_coefficient")

resultats <- data.frame(
  k = err_list$k,
  erreur = (1 - err_list$exactitude),
  sensibilite = err_list$sensibilite
)

# ==============================================================
# ------------------------- Thème ------------------------------
# ==============================================================

theme_projet <- theme_minimal(base_size = 13) +
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
    axis.title.x = element_text(
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      margin = margin(r = 12)
    ),
    axis.text = element_text(
      color = "#2f3e36"
    ),
    plot.title = element_text(
      face = "bold",
      size = 18,
      color = "#2f3e36",
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    plot.margin = margin(20, 25, 20, 25)
  )

# ==============================================================
# --------------- Taux d'erreur en fonction de k ---------------
# ==============================================================

# Déterminer le k minimisant l'erreur
index_min <- which.min(resultats$erreur)

k_erreur <- resultats$k[index_min]
erreur_min <- resultats$erreur[index_min] * 100

graphique_erreur <- ggplot(
  data = resultats,
  aes(x = k, y = erreur * 100)
) +
  
  # Courbe principale
  geom_line(color = "#455f51", linewidth = 1.5) +
  
  # Ligne pointillée horizontale
  annotate(
    "segment",
    x = 0, y = erreur_min,
    xend = k_erreur, yend = erreur_min,
    linetype = "dashed",
    color = "#2f3e36",
    linewidth = 0.9
  ) +
  
  # Ligne pointillée verticale
  annotate(
    "segment",
    x = k_erreur, y = erreur_min,
    xend = k_erreur, yend = 9,
    linetype = "dashed",
    color = "#2f3e36",
    linewidth = 0.9
  ) +
  
  # Point optimal
  annotate(
    "point",
    x = k_erreur,
    y = erreur_min,
    color = "#455f51",
    size = 4
  ) +
  
  # Mise en évidence du k optimal
  annotate(
    "label",
    x = k_erreur,
    y = 9,
    label = k_erreur,
    fill = "#455f51",
    color = "white",
    fontface = "bold",
    size = 4.5,
    linewidth = 0
  ) +
  
  # Mise en évidence de l'erreur minimale
  annotate(
    "label",
    x = 0,
    y = erreur_min,
    label = paste0(round(erreur_min, 2), "%"),
    fill = "#455f51",
    color = "white",
    fontface = "bold",
    size = 4.5,
    linewidth = 0
  ) +
  
  # Titres
  labs(
    title = "Taux d'erreur en fonction de k",
    x = "Nombre de voisins (k)",
    y = "Taux d'erreur (%)"
  ) +
  
  # Graduations de l'axe X (nombres impairs)
  scale_x_continuous(
    breaks = seq(1, max(resultats$k), by = 2),
    expand = expansion(mult = c(0.08, 0.05))
  ) +
  
  # Graduations de l'axe Y
  scale_y_continuous(
    breaks = seq(9, 17, by = 1),
    labels = function(x) paste0(x, "%"),
    expand = expansion(add = 0.4)
  ) +
  
  # Limites de l'axe Y
  coord_cartesian(ylim = c(9, 17)) +
  
  theme_projet

# Affichage
print(graphique_erreur)


# ==============================================================
# ---------------- Sensibilité en fonction de k ----------------
# ==============================================================

# Déterminer le k maximisant la sensibilité
index_max <- which.max(resultats$sensibilite)

k_sensibilite <- resultats$k[index_max]
sensibilite_max <- resultats$sensibilite[index_max] * 100

graphique_sensibilite <- ggplot(
  data = resultats,
  aes(x = k, y = sensibilite * 100)
) +
  
  # Courbe principale
  geom_line(color = "#455f51", linewidth = 1.5) +
  
  # Ligne pointillée horizontale
  annotate(
    "segment",
    x = 0, y = sensibilite_max,
    xend = k_sensibilite, yend = sensibilite_max,
    linetype = "dashed",
    color = "#2f3e36",
    linewidth = 0.9
  ) +
  
  # Ligne pointillée verticale
  annotate(
    "segment",
    x = k_sensibilite, y = sensibilite_max,
    xend = k_sensibilite, yend = 78,
    linetype = "dashed",
    color = "#2f3e36",
    linewidth = 0.9
  ) +
  
  # Point optimal
  annotate(
    "point",
    x = k_sensibilite,
    y = sensibilite_max,
    color = "#455f51",
    size = 4
  ) +
  
  # Mise en évidence du k optimal
  annotate(
    "label",
    x = k_sensibilite,
    y = 78,
    label = k_sensibilite,
    fill = "#455f51",
    color = "white",
    fontface = "bold",
    size = 4.5,
    linewidth = 0
  ) +
  
  # Mise en évidence de la sensibilité maximale
  annotate(
    "label",
    x = 0,
    y = sensibilite_max,
    label = paste0(round(sensibilite_max, 2), "%"),
    fill = "#455f51",
    color = "white",
    fontface = "bold",
    size = 4.5,
    linewidth = 0
  ) +
  
  # Titres
  labs(
    title = "Sensibilité en fonction de k",
    x = "Nombre de voisins (k)",
    y = "Sensibilité (%)"
  ) +
  
  # Graduations de l'axe X (nombres impairs)
  scale_x_continuous(
    breaks = seq(1, max(resultats$k), by = 2),
    expand = expansion(mult = c(0.08, 0.05))
  ) +
  
  # Graduations de l'axe Y
  scale_y_continuous(
    breaks = seq(78, 86, by = 1),
    labels = function(x) paste0(x, "%"),
    expand = expansion(add = 0.4)
  ) +
  
  # Limites de l'axe Y
  coord_cartesian(ylim = c(78, 86)) +
  
  theme_projet

# Affichage
print(graphique_sensibilite)


# ==============================================================
# ---------------- Exportation des graphiques ------------------
# ==============================================================

dir.create("Graphiques", showWarnings = FALSE)

ggsave(
  "Graphiques/graphique_erreur.svg",
  plot = graphique_erreur,
  width = 25,
  height = 15,
  units = "cm",
  bg = "#d1dfce"
)

ggsave(
  "Graphiques/graphique_sensibilite.svg",
  plot = graphique_sensibilite,
  width = 25,
  height = 15,
  units = "cm",
  bg = "#d1dfce"
)