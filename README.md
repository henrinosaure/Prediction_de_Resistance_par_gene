# analysededonnees

Question: **Pouvons-nous déterminer la résistance d'une bactérie à un antibiotique particulier, sachant ses gènes de résistance?**

# Données

**à mettre dans un fichier Data/ :**

**Antimicrobial Resistance (le gêne)**
NCBI\_AMR : https://www.ncbi.nlm.nih.gov/pathogens/isolates/
garder seulement colonnes : "BioSample"           "AMR genotypes"

**Antibiotic sensibility test**
NCBI\_AST : https://www.ncbi.nlm.nih.gov/pathogens/ast/
garder seulement colonnes : "BioSample"            "Organism group"       "Scientific name"     "Antibiotic"     "Resistance phenotype"



# Paquets
Pour installer les paquets nécessaires au projet, lancer la commande suivante:

install.packages(c("dplyr", "future", "KernelKnn", "tidyverse", "FactoMineR", "factoextra", "future.apply"))

