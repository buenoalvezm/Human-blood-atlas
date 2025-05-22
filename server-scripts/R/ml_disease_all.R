#### Title: Machine learning with all diseases
#### Author: María Bueno Álvez
#### Description: script to perform multiclass classification including all diseases
#### Last edited : 12/08/2024

# Set locale
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])

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
dir.create(paste0("results_disease_all/"), showWarnings = FALSE)
dir.create(paste0("results_disease_all/model/"), showWarnings = FALSE)
set.seed(seed)

# Configure future plan for parallelization
n_workers <- parallel::detectCores() - 1
plan(multicore, workers = n_workers) 

# Load splits 
random_splits <- readRDS("random_splits.rds")
current_split <- random_splits[[seed]]
rm(random_splits)
gc()

# Define diseases included
include_diseases <- 
  current_split$data_train |>
  bind_rows(current_split$data_test) |> 
  distinct(Disease) |>
filter(Disease != "Healthy") |>
  pull()

# Run multiclass ML for all diseases
cat(paste0("Starting multiclassification (all diseases) run for seed ",seed))
t0 <-  Sys.time()
multiclass_ml_results_splits <-
  
  disease_against_all(
    split_train = current_split$data_train,
    split_test = current_split$data_test,
    path_save =  paste("results_disease_all/"),
seed = seed
  )

t1 <-  Sys.time()
gc()

# Save output and record run time
saveRDS(multiclass_ml_results_splits, paste("results_disease_all/model/all_ml_results_splits_", seed, ".rds", sep = ""))
rm(multiclass_ml_results_splits)
gc()
t1-t0

