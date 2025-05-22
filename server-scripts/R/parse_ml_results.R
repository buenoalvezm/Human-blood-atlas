
library(tidyverse)
dir.create("server-scripts/parsed-results", showWarnings = FALSE)

n_seeds <- 100

# Parse age
common_path <- "server-scripts/results_age"

importances_age <- 
  map_df(c(1:n_seeds), \(seed) {
    read_tsv(paste0(common_path, "/important_proteins_Age_seed_", seed, ".tsv")) |> 
      mutate(Type = "Age",
             Seed = seed)
  }) 
      
write_tsv(importances_age, "server-scripts/parsed-results/importances_age.tsv")

predictions_age <- 
  map_df(c(1:n_seeds), \(seed) {
    
   # age_split <- 
    #  readRDS(paste0(common_path, "/age_split_", seed, ".rds")) 
    
    predictions <- 
      read_tsv(paste0(common_path, "/predictions_Age_seed_", seed, ".tsv")) |> 
      mutate(Type = "Sex",
             Seed = seed)#,
            # DAid = age_split$sub_split$data_test$DAid) |> 
     # relocate(DAid)
    
  }) 

write_tsv(predictions_age, "server-scripts/parsed-results/predictions_age.tsv")


# Parse sex
common_path <- "server-scripts/results_sex"

importances_sex <- 
  map_df(c(1:n_seeds), \(seed) {
    read_tsv(paste0(common_path, "/important_proteins_Sex_seed_", seed, ".tsv")) |> 
      mutate(Type = "Age",
             Seed = seed)
  }) 

write_tsv(importances_sex, "server-scripts/parsed-results/importances_sex.tsv")

predictions_sex <- 
  map_df(c(1:n_seeds), \(seed) {
    
    # age_split <- 
    #  readRDS(paste0(common_path, "/age_split_", seed, ".rds")) 
    
    predictions <- 
      read_tsv(paste0(common_path, "/predictions_Sex_seed_", seed, ".tsv")) |> 
      mutate(Type = "Sex",
             Seed = seed)#,
    # DAid = age_split$sub_split$data_test$DAid) |> 
    # relocate(DAid)
    
  }) 

write_tsv(predictions_sex, "server-scripts/parsed-results/predictions_sex.tsv")


# Parse BMI
common_path <- "server-scripts/results_age"

importances_age <- 
  map_df(c(1:n_seeds), \(seed) {
    read_tsv(paste0(common_path, "/important_proteins_Age_seed_", seed, ".tsv")) |> 
      mutate(Type = "Age",
             Seed = seed)
  }) 

write_tsv(importances_age, "server-scripts/parsed-results/importances_age.tsv")

predictions_age <- 
  map_df(c(1:n_seeds), \(seed) {
    
    # age_split <- 
    #  readRDS(paste0(common_path, "/age_split_", seed, ".rds")) 
    
    predictions <- 
      read_tsv(paste0(common_path, "/predictions_Age_seed_", seed, ".tsv")) |> 
      mutate(Type = "Age",
             Seed = seed)#,
    # DAid = age_split$sub_split$data_test$DAid) |> 
    # relocate(DAid)
    
  }) 

write_tsv(predictions_age, "server-scripts/parsed-results/predictions_age.tsv")


# Parse disease
common_path <- "server-scripts/results_disease_all/"

n_seeds <- 100

importances_disease_all <- 
  
  map_df(c(1:n_seeds), \(seed) {
    read_tsv(paste0(common_path, "important_proteins_multiclass_seed_", seed, ".tsv")) |> 
      mutate(Class = "Multiclass", 
             Type = "All other diseases",
             Seed = seed)
  }) 

write_tsv(importances_disease_all, "parsed-results/importances_disease_all.tsv")

predictions_disease_all <- 
  map_df(c(1:n_seeds), \(seed) {
    read_tsv(paste0(common_path, "predictions_multiclass_seed_", seed, ".tsv")) |> 
      mutate(Class = "Multiclass", 
             Type = "All other diseases",
             Seed = seed)
  }) 


write_tsv(predictions_disease_all, "parsed-results/predictions_disease_all.tsv")

roc_disease_all <- 
  
  map_df(c(1:n_seeds), \(seed) {
    read_tsv(paste0(common_path, "roc_multiclass_seed_", seed, ".tsv")) |> 
      mutate(Class = "Multiclass", 
             Type = "All other diseases",
             Seed = seed)
  }) 


write_tsv(roc_disease_all, "parsed-results/roc_disease_all.tsv")

