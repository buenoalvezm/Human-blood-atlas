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
resource_meta <- read_csv("data/disease_meta_phase1.tsv")


# Create directory to save results and set seed 
dir.create("results_sex", showWarnings = FALSE)
dir.create("results_sex/model/", showWarnings = FALSE)
set.seed(seed)

# Configure future plan for parallelization
n_workers <- parallel::detectCores() - 1
plan(multicore, workers = n_workers) 

#Sex
cat("Preparing Sex splits")

sex_ml_split <- readRDS("data/sex_ml_split.rds")

sex_model_split <-
  generate_split(data = sex_ml_split$data_train,
                 proportion = 0.8,
                 seed = seed,
                 variable_stratify = "Sex")

sex_ml_splits <- list("master_split" = sex_ml_split,
                      "sub_split" = sex_model_split)

rm(sex_ml_split, sex_model_split)

cat("Starting Sex prediction")

sex_ml_results_splits <-

  discrete_prediction(
    variable_predict = "Sex",
    variable_case = "F",
    split_train = sex_ml_splits$sub_split$data_train,
    split_test = sex_ml_splits$sub_split$data_test,
    path_save =  "results_sex/",
    seed = seed
  )

saveRDS(sex_ml_results_splits, paste("results_sex/sex_res_", seed, ".rds", sep = ""))
saveRDS(sex_ml_splits, paste("results_sex/model/sex_split_", seed, ".rds", sep = ""))

rm(sex_ml_results_splits, sex_ml_splits)
gc()
