# Author : Henry Cloutier
# Date : nov 2026
#
#
#    Description
#  ═══════════════
#  Fonction permettant de nettoyer, filtrer et combiner les banques
#  de données AST et AMR du NCBI afin de préparer les données pour
#  l'analyse de la résistance aux antibiotiques.
#
#  La fonction effectue les opérations suivantes :
#  - Combine les données phénotypiques et génotypiques par BioSample.
#  - Garde uniquement les phénotypes susceptible et resistant.
#  - Enlève les observations sans génotype AMR renseigné.
#  - Sépare les génotypes en gènes individuels.
#  - Enlève les gènes présents chez moins de 10 isolats.
#  - Garde les 25 antibiotiques les plus fréquents.
#  - Enlève les doublons de combinaisons BioSample/Antibiotic/gene.
#  - Garde uniquement les observations pour la ciprofloxacine.
#  - Convertit les données en format large et binaire.
#  - Enlève les observations ayant des profils génétiques identiques.
#
#
#    Arguments
#  ═════════════
# • AST (data.frame)
#     Banque de données phénotypiques de résistance aux antibiotiques.
#     Colonnes nécessaires :
#     - BioSample
#     - Organism group
#     - Antibiotic
#     - Resistance phenotype
#
# • AMR (data.frame)
#     Banque de données génotypiques de résistance aux antibiotiques.
#     Colonnes nécessaires :
#     - BioSample
#     - AMR genotypes
#       Gènes séparés par des virgules, sous la forme gene=detection.
#
#
#    Returns
#  ═══════════
#  Un tibble en format large contenant :
#
#  - BioSample : identifiant unique de l'échantillon.
#
#  - Resistance.phenotype : phénotype de résistance à la ciprofloxacine.
#       0 = susceptible
#       1 = résistant
#
#  - Colonnes de gènes : une colonne par gène retenu.
#       0 = absence du gène dans les annotations
#       1 = présence du gène dans les annotations
#
#  Chaque ligne correspond à une observation ayant un profil
#  génétique unique.
#
#
#    Dependencies
#  ════════════════
#  tidyverse
#

#install.packages("tidyverse")
library(tidyverse)
filtre_donnes <- function(AST, AMR, antibiotique = "ciprofloxacin") {
  # =============================================================================
  # ---------- Lecture et clean/réorganisation des fichiers de données ----------
  # =============================================================================
  
  # Filtre et combine les 2 banques de donnees, gardant les colonnes :
  # BioSample, Organism Group, Antibiotic, Resistance.phenotype et AMR Genotypes
  NCBI_filtered <- AST %>%
    filter(`BioSample` %in% AMR$`BioSample`) %>%
    filter(`Resistance phenotype` %in% c("susceptible", "resistant")) %>%
    select(`BioSample`, `Organism group`, Antibiotic, 
           `Resistance.phenotype`=`Resistance phenotype`) %>%
    left_join(AMR %>% select(`BioSample`, `AMR genotypes`), by = "BioSample", 
              relationship = "many-to-many") %>%
    filter(
      !is.na(`AMR genotypes`),
      trimws(`AMR genotypes`) != ""
    )
  
  # Split la banque de donnees pour chaque gene
  NCBI_genes <- NCBI_filtered %>%
    separate_longer_delim(`AMR genotypes`, delim = ",") %>%
    separate_wider_delim(
      `AMR genotypes`,
      delim = "=",
      names = c("gene", "detection")
    )
  
  # =============================================================================
  # -------------------- Filtration pour données de qualité ---------------------
  # =============================================================================
  
  # Enlever les donnees avec des gènes trop niches (apparaissent moins de 10 fois)
  # Garder seulement les 25 antibiotiques les plus fréquents
  
  # Trouve les genes a enlever en fonction de ceux qui sont les plus frequents 
  # (on garde seulement ceux qui ont 10 ou plus d'apparitions)
  genes_to_remove <- NCBI_genes %>%
    distinct(`BioSample`, gene) %>%
    count(gene) %>%
    # Changer parametre à côté de n pour seuil de fois que le gêne doit être la pour être gardé
    filter(n < 10)
  
  # Garde seulement les 25 antibiotiques les plus fréquents
  antibio_to_keep <- NCBI_genes %>%
    distinct(BioSample, Antibiotic) %>%
    count(Antibiotic) %>%
    # Changer n pour nombre d'antibiotiques à garder
    slice_max(n, n = 25, with_ties = FALSE)
  
  # Applique les 2 banques de filtration precedentes et enleve les observations en double 
  #(Qui ont tous le meme biosample, gene et antibiotique testé)
  NCBI_genes_filtered <- NCBI_genes %>%
    filter(!(gene %in% genes_to_remove$gene)) %>%
    filter(Antibiotic %in% antibio_to_keep$Antibiotic) %>%
    distinct(BioSample, Antibiotic, gene, .keep_all = TRUE)
  
  # =============================================================================
  # -------------------- préparer les données pour l'ACM ------------------------
  # =============================================================================
  
  # Enlève tous les antibiotiques qui ne sont pas celui choisi
  # ciprofloxacin par défaut
  genes_ciprofloxacin <- NCBI_genes_filtered %>%
    filter(Antibiotic == antibiotique)
  
  # On veut mettre la banque en mode large pour l'ACM : 
  # La colonne Resistance.phenotype est en binaire : susceptible = 0 et resistant = 1
  # Une colonne désigne la présence ou non du gene
  
  genes_ACM_ciprofloxacin_wide <- genes_ciprofloxacin %>%
    mutate(`Resistance.phenotype` = if_else(`Resistance.phenotype` == "susceptible", 0, 1)) %>%
    mutate(presence = 1) %>%
    pivot_wider(
      id_cols = c(BioSample, `Resistance.phenotype`),
      names_from = gene,
      values_from = presence,
      values_fill = 0
    ) %>%
    # Enlève les doublons de gènes (observations qui ont les mêmes valeurs pour tous les gènes)
    distinct(across(-c(BioSample, Resistance.phenotype)), .keep_all = TRUE)
  
  # Valeur de retour
  genes_ACM_ciprofloxacin_wide
}