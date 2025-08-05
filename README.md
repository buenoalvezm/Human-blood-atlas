
# Pan-disease Human Blood Atlas 

![GitHub top language](https://img.shields.io/github/languages/top/buenoalvezm/Human-blood-atlas)
![License](https://img.shields.io/badge/license-Apache2.0-yellow)
![GitHub last commit](https://img.shields.io/github/last-commit/buenoalvezm/Human-blood-atlas)
![GitHub issues](https://img.shields.io/github/issues/buenoalvezm/Human-blood-atlas)

## Table of contents

1. [Project overview](#project-overview)
2. [Data description](#data-description)
3. [Methodology](#methodology)
4. [Usage](#usage)
5. [Directory structure](#directory-structure)
6. [Results & insights](#results--insights)
7. [License](#license)
8. [Contact](#contact)

## Project overview

This repository contains the code and resources for the study *A human pan-disease blood atlas of the circulating proteome*, which presents a large-scale analysis of blood protein profiles across various diseases and longitudinal healthy cohorts aimed at understanding human health states. The study leverages Proximity Extension Assay (PEA) proteomics to identify proteins that vary across human development, associated to demographic variables susch as age, sex and BMI, and disease-specific and shared biomarkers. An open-access **Disease Blood Resource** has been created as part of the Human Protein Atlas (www.proteinatlas.org/humanproteome/blood), enabling researchers to explore protein profiles across diseases and healthy individuals.

## Data description

This study is based on plasma proteome profiles from over 8,000 individuals, encompassing both healthy individuals and patients across a wide spectrum of diseases. Key characteristics of the dataset include:

- **Cohort composition**: samples were collected from multiple biobank cohorts, including healthy longitudinal cohorts (BAMSE, Wellness), cross-sectional disease cohorts (Human Disease Blood Atlas, U-CAN), and a prospective disease cohort (UKB-PPP).
- **Proteomic platform**: protein profiling was performed using Proximity Extension Assay (PEA), a high-throughput, antibody-based method that enables simultaneous quantification of thousands of proteins from minimal plasma volumes.

## Methodology

The analysis employed a combination of univariate and multivariate techniques to uncover key health and disease-related protein patterns. Key components include:

- **Exploratory Data Analysis (EDA)**: Uniform Manifold Approximation and Projection (UMAP) and correlation analyses were performed to uncover major trends in the dataset.
- **Analyses of variance**: statistical analyses were applied to identify proteins correlated with disease and demographic factors (e.g., age, sex, BMI). 
- **Differential abundance analysis**: applied to all measured proteins to identify significant changes in abundance across diseases.
- **Machine learning**: lasso-based models were applied for disease classification and prediction of biological age, sex, and BMI.

## Results & insights

Several key insights were derived from the analysis of the dataset:

- **Stability of the plasma proteome across human lifespan**: the most pronounced changes in protein levels occurred during puberty and the transition to adulthood, while the adult plasma proteome remained remarkably stable over time.
- **Age and disease correlations**: strong correlations were found between protein levels and chronological age, but the presence of disease emerged as the main contributor to variation in protein profiles.
- **Biological age prediction**: machine learning models trained on plasma protein profiles accurately predicted biological age, highlighting proteins with potential relevance for aging and age-related diseases.
- **Insights from a pan-disease perspective**: By integrating differential abundance and machine learning analyses across 59 diseases, we identified proteins with distinct disease-specific profiles as well as those consistently elevated across multiple disease categories, such as cancer and autoimmune disorders. This pan-disease approach enabled us to distinguish proteins with potential diagnostic relevance from those involved in broader systemic responses, such as inflammation.

 
## Usage

To explore the code used for the analysis described in the manuscript:

1. Clone this repository:
   ```bash
   git clone https://github.com/buenoalvezm/Human-blood-atlas.git
   cd Human-blood-atlas
   ```
2. Install R dependencies using renv:
   ```R
   renv::restore()
   ```

## Directory structure

📁 Pan-disease-profiling/  
├── 📂 scripts/            # Analysis scripts (R)  
├── 📂 server-scripts/     # Analysis scripts for HPC (R/shell)  
├── 📜 README.md           # Overview of the project  
└── 📜 LICENSE             # License file

## License

This project is licensed under the Apache License 2.0 - see the LICENSE file for details.

## Contact

For questions or inquiries, please contact:

María Bueno Álvez  
GitHub: [buenoalvezm](https://github.com/buenoalvezm)
 
