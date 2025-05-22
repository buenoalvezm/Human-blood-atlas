library(tidyverse)

import_df <- function(file_path) {
  
  # Determine file extension from file path
  file_extension <- tools::file_ext(file_path)
  
  df <- switch(tolower(file_extension),
               csv = utils::read.csv(file_path, stringsAsFactors = FALSE),
               tsv = utils::read.delim(file_path, stringsAsFactors = FALSE),
               txt = utils::read.table(file_path, header = TRUE, stringsAsFactors = FALSE),
               rda = { load(file_path); get(ls()[1]) },
               rds = readRDS(file_path),
               xlsx = readxl::read_excel(file_path, guess_max=10000000),
               parquet = arrow::read_parquet(file_path),
               stop("Unsupported file type: ", file_extension))
  
  df <- tibble::as_tibble(df)
  return(df)
}

savepath <- 
  function(savename) { 
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_folder <- 
  function(folder, savename) { 
    result_folder <- paste0("results/", Sys.Date(), "/",folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_data <- 
  function(folder, savename) { 
    result_folder <- paste0("data/processed/", folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

savepath_results <- 
  function(folder, savename) { 
    
    dir.create("results/", showWarnings = FALSE)
    result_folder <- paste0("results/", folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    return(savename)
    
  }
