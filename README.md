
# Pan-Disease Blood Proteome Profiling

![GitHub top language](https://img.shields.io/github/languages/top/buenoalvezm/Pan-disease-profiling)
![License](https://img.shields.io/badge/license-Apache2.0-yellow)
![GitHub last commit](https://img.shields.io/github/last-commit/buenoalvezm/Pan-disease-profiling)
![GitHub issues](https://img.shields.io/github/issues/buenoalvezm/Pan-disease-profiling)

## Table of Contents

1. [Project Overview](#project-overview)
2. [Data Description](#data-description)
3. [Methodology](#methodology)
4. [Usage](#usage)
5. [Directory Structure](#directory-structure)
6. [Results & Insights](#results--insights)
7. [License](#license)
8. [Contact Information](#contact-information)

## Project Overview

This project presents a comprehensive blood proteome profiling atlas, aimed at understanding human health states by analyzing proteomic signatures across various diseases. The core objective of the project is to identify disease-specific biomarkers through machine learning and integrate these findings into an open-access resource.

This research integrates a pan-disease strategy across major disease categories such as cancer, autoimmune, cardiovascular, metabolic, and infectious diseases, leveraging state-of-the-art targeted proteomics using the Proximity Extension Assay (PEA). This work also incorporates longitudinal analysis of healthy individuals from childhood through adulthood, providing insights into the stability and variability of protein profiles.

This repository contains the code and resources for the study "Pan-Disease Blood Proteome Profiling", which presents a large-scale analysis of blood protein profiles across multiple human diseases. The study leverages Proximity Extension Assay (PEA) proteomics to identify disease-specific and shared biomarkers, contributing to a deeper understanding of disease biology and potential translational applications.

An open-access Human Disease Blood Atlas has been created as part of the Human Protein Atlas (www.proteinatlas.org), enabling researchers to explore protein profiles across diseases, healthy individuals, and different demographic factors.

## Data Description

The dataset includes plasma proteomic profiles from over 8,000 individuals, representing diverse disease categories as well as healthy individuals. Key characteristics of the dataset are:

- **Source of Data**: Collected from multiple biobank cohorts, covering 59 diseases and healthy samples.
- **Proteomics Approach**: Proximity Extension Assay (PEA) was used for high-throughput analysis, enabling quantification of thousands of proteins from small blood sample volumes.
- **Proteomics Dataset**: the dataset encompasses up to 5,400 protein measurements across 59 diseases and healthy cohorts, with longitudinal samples from healthy individuals to monitor protein variability over time and disease-related changes.
- **Data Processing**: quality control (QC) was performed to ensure data reliability and consistency.

## Methodology

The analysis employed a combination of univariate and multivariate techniques to uncover key health and disease-related protein patterns. Key components include:

- **Exploratory Data Analysis (EDA)**: Principal component analysis (PCA) and correlation analyses were performed to uncover major trends in the dataset.
- **Analyses of variance**: Statistical analyses  were applied to identify proteins correlated with disease and demographic factors (e.g., age, sex, BMI). 
- **Differential Expression Analysis**: Identified proteins that show significant changes in abundance across disease categories.
- **Machine Learning Classification**: Lasso-based models were applied for disease classification and biological age prediction.

## Results & Insights

Several key insights were derived from the analysis of the dataset:

- **Disease-specific biomarkers**: distinct protein profiles were identified for specific diseases. Some proteins were elevated across multiple diseases, suggesting their role as markers of inflammation rather than specific diseases.
- **Age and disease correlations**: strong correlations were found between protein levels and both biological and chronological age. Disease presence was found to be the primary factor contributing to variation in protein profiles.
- **Cross-disease patterns**: certain proteins exhibited elevated levels across multiple disease categories, such as cancer and autoimmune diseases, indicating shared biological mechanisms.
- **Disease-specific signatures**: we identified proteins associated to each of the analyzed diseases.
  
## Usage

To run the analysis and reproduce the results:

1. Clone this repository:
   ```bash
   git clone https://github.com/buenoalvezm/Pan-disease-profiling.git
   cd Pan-disease-profiling
   ```
2. Install R dependencies using renv:
   ```R
   renv::restore()
   ```

## Directory Structure

📁 Pan-disease-profiling/  
├── 📂 data/               # Processed datasets (or links to access them)  
├── 📂 scripts/            # Analysis scripts (R)  
├── 📂 server-scripts/     # Analysis scripts for HPC (R)  
├── 📜 README.md           # Overview of the project  
└── 📜 LICENSE             # License file

## License

This project is licensed under the Apache License 2.0 - see the LICENSE file for details.

## Contact

For questions or inquiries, please contact:

María Bueno Álvez  
GitHub: [buenoalvezm](https://github.com/buenoalvezm)
 
