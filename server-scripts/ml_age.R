# Set locale
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# Read in arguments
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])

# Load libraries and functions
library(tidyverse)
library(tidymodels)
library(furrr)
library(themis)
library(future.callr)

# Read functions and data
source("R/functions_ml.R", encoding="utf-8")
resource_meta <- read_csv("data/meta_phase1_pandisease_20250226.tsv")


# Create directory to save results and set seed 
dir.create("results_age", showWarnings = FALSE)
dir.create("results_age/model/", showWarnings = FALSE)
set.seed(seed)

# Configure future plan for parallelization
n_workers <- parallel::detectCores() - 1
plan(multicore, workers = n_workers) 


cat("Preparing Age splits")
age_ml_split <- readRDS("data/age_ml_split.rds")


age_model_split <- 
  generate_split(data = age_ml_split$data_train, 
                 proportion = 0.8,
                 seed = seed,
                 variable_stratify = "Age")

age_ml_splits <- list("master_split" = age_ml_split,
                      "sub_split" = age_model_split)
rm(age_ml_split, age_model_split)

cat("Starting Age prediction")

age_ml_results_splits <-
  
  continuous_prediction(
    split_train = age_ml_splits$sub_split$data_train,
    split_test = age_ml_splits$sub_split$data_test,
    variable_predict = "Age",
    path_save =  "results_age/",
    seed = seed
  )

saveRDS(age_ml_results_splits, paste("results_age/age_res_", seed, ".rds", sep = ""))
saveRDS(age_ml_splits, paste("results_age/age_split_", seed, ".rds", sep = ""))

rm(age_ml_results_splits, age_ml_splits)
gc()


