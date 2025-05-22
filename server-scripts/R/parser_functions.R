#### Title: Parser functions
#### Author: María Bueno Álvez
#### Description: script collecting general functions to parse machine learning results generated in the server
#### Last edited : 03/03/2025


parse_predictions <- function(n_seeds, type) {
  
  map_df(c(1:n_seeds), \(seed) {
    
    # age_split <-
    #  readRDS(paste0(common_path, "/age_split_", seed, ".rds"))
    
    predictions <- 
      read_tsv(paste0(common_path, "/predictions_", toupper(type), "_seed_", seed, ".tsv")) |> 
      mutate(Type = type,
             Seed = seed)#,
    # DAid = age_split$sub_split$data_test$DAid) |> 
    # relocate(DAid)
    
  }) 
  
}
 


pares_importances <- function(n_seeds, type) {
  
  map_df(c(1:n_seeds), \(seed) {
    read_tsv(paste0(common_path, "/important_proteins_", toupper(type), "_seed_", seed, ".tsv")) |> 
      mutate(Type = type,
             Seed = seed)
  }) 
  
  
} 
