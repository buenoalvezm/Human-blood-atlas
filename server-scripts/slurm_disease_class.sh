#!/bin/bash

#SBATCH -A sens2018115          
#SBATCH -p core                  
#SBATCH --cpus-per-task=16       
#SBATCH -t 30:00:00              
#SBATCH -J disease_class   
#SBATCH --array=1-100            
#SBATCH --mail-type=ALL          
#SBATCH --mail-user=mariba@kth.se
#SBATCH --output=slurm_output/%A_disease_class_%a.out   

# Load necessary modules
module load R_packages/4.2.1    

# Define the classes and seed
classes=("Cardiovascular" "Metabolic" "Cancer" "Autoimmune" "Infection")
seed=${SLURM_ARRAY_TASK_ID}

# Iterate over each class and run the R script for binary classification and multiclass (class)
for class in "${classes[@]}"; do
  echo "Processing class: $class, seed: $seed"
  Rscript R/ml_disease_class.R $seed $class
done
