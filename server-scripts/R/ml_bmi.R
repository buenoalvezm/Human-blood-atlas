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
dir.create("results_bmi", showWarnings = FALSE)
dir.create("results_bmi/model/", showWarnings = FALSE)
set.seed(seed)

# Configure future plan for parallelization
n_workers <- parallel::detectCores() - 1
plan(multicore, workers = n_workers)

# BMI
cat("Preparing BMI splits")

bmi_ml_split <- readRDS("data/bmi_ml_split.rds")

bmi_model_split <-
  generate_split(data = bmi_ml_split$data_train,
                 proportion = 0.8,
                 seed = seed,
                 variable_stratify = "BMI")

bmi_ml_splits <- list("master_split" = bmi_ml_split,
                      "sub_split" = bmi_model_split)

rm(bmi_ml_split, bmi_model_split)

cat("Starting BMI prediction")

bmi_ml_results_splits <-
  
  continuous_prediction(
    split_train = bmi_ml_splits$sub_split$data_train,
    split_test = bmi_ml_splits$sub_split$data_test,
    variable_predict = "BMI",
    path_save =  "results_bmi/",
    seed = seed
  )

saveRDS(bmi_ml_results_splits, paste("results_bmi/bmi_res_", seed, ".rds", sep = ""))
saveRDS(bmi_ml_splits, paste("results_bmi/model/bmi_split_", seed, ".rds", sep = ""))

rm(bmi_ml_results_splits, bmi_ml_splits)
gc()
