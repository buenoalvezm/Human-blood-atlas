#### Title: Machine learning within class
#### Author: María Bueno Álvez
#### Description: script to perform ML analyses comparing each disease against healthy and within class classification
#### Last edited : 11/10/2024

# Set locale
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
class <- as.character(args[2])

# Load libraries and functions
library(tidyverse)
library(tidymodels)
library(themis)
library(multiROC)
library(furrr)
library(future.callr)

# Read functions and data
source("functions_ml.R", encoding="utf-8")
resource_meta <- read_tsv("meta_phase1_pandisease_20250226.tsv")

# Create directory to save results and set seed 
dir.create(paste0("results_disease_class"), showWarnings = FALSE)
dir.create(paste0("results_disease_class/model/"), showWarnings = FALSE)
dir.create(paste0("results_disease_class/", class), showWarnings = FALSE)
set.seed(seed)

# Configure future plan for parallelization
n_workers <- parallel::detectCores() - 1
plan(multicore, workers = n_workers) 

# Load splits
random_splits <- readRDS("random_splits.rds")
current_split <- random_splits[[seed]]
rm(random_splits)
gc()

# Define diseases to run 
include_diseases <- 
  current_split$data_train |>
  bind_rows(current_split$data_test) |> 
  distinct(Disease) |> 
  left_join(resource_meta |>
              distinct(Disease, Class), by = "Disease") |>
  filter(Class == class) |>
  pull(Disease)

# Run binary ML for each disease against the class
cat(paste0("Starting binary (class) runs for seed ",seed))
t0 <-  Sys.time()

if(class != "Psychiatric") {
  
  class_ml_results_splits <-
    disease_against_class(
      class = class,
      split_train = current_split$data_train,
      split_test = current_split$data_test,
      path_save =  paste0("results_disease_class/", class),
      seed = seed
    )
  
  
} else {
  class_ml_results_splits <-  NULL
}

t1 <-  Sys.time()
gc()

# Save output and record run time
saveRDS(class_ml_results_splits, paste("results_disease_class/model/", class, "_multiclass_ml_results_splits_", seed, ".rds", sep = ""))
rm(class_ml_results_splits)
gc()
t1-t0

