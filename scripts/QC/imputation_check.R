#### Title: Missing value imputation exploration
#### Author: María Bueno Álvez
#### Description: script to investigate the extent of missing values in the dataset 
#### Last edited : 14/05/2025

source("scripts/functions/functions_utility.R")
source("scripts/functions/functions_visualization.R")
source("scripts/functions/themes_palettes.R")

# Look at control samples - to explore batch effects

data_wellness <- read_tsv("../Human-disease-blood-atlas/data/final_data/HPA/v24_2/wellness_data_ht_phase2.tsv")
data_bamse <- read_tsv("../Human-disease-blood-atlas/data/final_data/HPA/v24_2/bamse_data_ht_phase2.tsv")
data_disease <- read_tsv("../Human-disease-blood-atlas/data/final_data/HPA/v24_2/disease_data_phase1.tsv")


# Function to look at extent of imputation
do_imputation_check <- function(data,
                                dataset) {
  
  message("Processing dataset: ", dataset)
  
  # Prepare data in wide format
  data_w <- 
    data |> 
    select(DAid, Assay, NPX) |> 
    pivot_wider(names_from = Assay,
                values_from = NPX)
 
   # Imputation recipe
  impute_rec <-
    recipe( ~ ., data = data_w) %>%
    update_role(DAid, new_role = "id")  |>
    step_normalize(all_predictors()) |>
    step_impute_knn(all_predictors()) 
  
  impute_rec_prep <- prep(impute_rec)
  
  data_imputed <- bake(impute_rec_prep, new_data = NULL)
  
  # Compare number of imputed values
  missing_before <- sapply(data_w, function(x) sum(is.na(x)))
  missing_after <- sapply(data_imputed, function(x) sum(is.na(x)))
  
  # Difference gives number of values imputed
  imputed_counts <- missing_before - missing_after
  
  # Get total number & percentage of imputed values
  total_imputed <- sum(imputed_counts)
  total_values <- prod(dim(data_imputed))
  percent_imputed <- total_imputed / total_values * 100
  
  # Print results
  cat("Total imputed values:", total_imputed, "\n")
  cat("Percentage of imputed values:", round(percent_imputed, 2), "%\n")
}

# Check for all three datasets
do_imputation_check(data = data_wellness, 
                    dataset = "Wellness")

do_imputation_check(data = data_bamse, 
                    dataset = "BAMSE")


do_imputation_check(data = data_disease, 
                    dataset = "HDBA")

