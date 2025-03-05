#!/bin/bash

#SBATCH -A sens2018115         
#SBATCH -p core                 
#SBATCH --cpus-per-task=16       
#SBATCH -t 30:00:00              
#SBATCH -J ml_bmi       
#SBATCH --array=1-100            
#SBATCH --mail-type=ALL          
#SBATCH --mail-user=mariba@kth.se
#SBATCH --output=slurm_output/%A_bmi_%a.out   


echo $SLURM_ARRAY_TASK_ID

module load R_packages/4.2.1    # Load necessary modules

seed=${SLURM_ARRAY_TASK_ID}

# Run the R script with the seed parameter
Rscript R/ml_bmi.R $seed
