# NF1-long-read
NF1 Single-cell RNA isoform analysis

Arborescence :

.
├── 1_metadata\
│   └── 1_make_metadata.Rmd\
├── 2_individual\
│   ├── 1_make_individual.Rmd\
│   ├── dataset\
│   └── input\
│       ├── TITI\
│       │   ├── Illumina\
│       │   └── Nanopore\
│       └── TOTO\
│           ├── Illumina\
│           └── Nanopore\
├── 3_combined\
│   └── 1_combined_all.Rmd\
├── 4_isoforms\
│   ├── 1_analysis.Rmd\
│   └── 2_visualisation.Rmd\
└── README.md\


**1_metadata** :*
- nécessaire pour l'annotation des types cellulaires
- 1_metadata.Rmd : notebbok pour définir des couleurs pour les échantillons

**2_individual** : 
- **input** : un dossier par échantillon (TOTO et TITI), avec les matrices de comptage (un dossier "genes" et un "isoforms" pour chaque échantillon), faire des liens symboliques si les données sont ailleurs, pour éviter les duplications
- **dataset** : dossier de sortie avec les objets Seurat au format RDS
- 1_make_individual.Rmd : notebook pour générer les objets Seurat individuels

**3_combined** : 
- 1_combined_all.Rmd : notebook pour intégrer les deux jeux de données

**4_isoforms** : 
- 1_analysis.Rmd : notebook pour faire l'analyse des isoformes différentiellement exprimés
- 2_visualisation.Rmd : notebook pour faire beaucoup de figures d'exploration et analyses des résultats

Commande pour kniter dans le terminal, avec un paramètre :

cd ./2_individual
sample="TOTO" ; file="1_make_individual" ; Rscript -e "rmarkdown::render(input = '${file}.Rmd', output_file = '${sample}.html', params = list(sample_name = '${sample}'))"
