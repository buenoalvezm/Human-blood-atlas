
#### Title: Mixed model exploration
#### Author: María Bueno Álvez
#### Description: script to investigate more complex mixed models
#### Last edited : 14/05/2025

source("scripts/functions/functions_analyses.R")
source("scripts/functions/functions_utility.R")
source("scripts/functions/functions_visualization.R")
source("scripts/functions/themes_palettes.R")

dat_wellness <- readRDS(savepath_data("MM_test", "data_mm_Wellness.rds"))
dat_bamse <- readRDS(savepath_data("MM_test", "data_mm_BAMSE.rds"))

# Functions
test_mm_interaction <- function(df, protein) {
  
  df <- 
    df |> 
    filter(Assay == protein) |> 
    select(DAid, Age, Sex, Subject, NPX) |> 
    mutate(Age = as.numeric(Age),
           Sex = as.factor(Sex),
           Subject = as.factor(Subject))
  
  
  # Fit the mixed-effects model (no random slopes, random intercepts)
  m1 <- lmer(NPX ~ Age + Sex + (1 | Subject), data = df) 
  m2 <- lmer(NPX ~ Age * Sex + (1 | Subject), data = df) 
  
  models <- list(m1, m2)
  
  singular_flags <- sapply(models, function(m) {
    tryCatch(isSingular(m, tol = 1e-4), error = function(e) NA)
  })
  
  tidy(anova(m1,m2,refit=FALSE)) |> 
    mutate(Protein = protein,
           Singular = singular_flags) |> 
    relocate(Protein)
  
}


test_mm_slopes <- function(df, protein) {
  
  df <- 
    df |> 
    filter(Assay == protein) |> 
    select(DAid, Age, Sex, Subject, NPX) |> 
    mutate(Age = as.numeric(Age),
           Sex = as.factor(Sex),
           Subject = as.factor(Subject))
  
  
  # Fit the mixed-effects model (no random slopes, random intercepts)
  m1 <- lmer(NPX ~ Age * Sex + (1 | Subject), data = df) 
  m2 <- lmer(NPX ~ Age * Sex + (1 + Age | Subject), data = df) 
  m3 <- lmer(NPX ~ Age * Sex + (1 + Age + Sex | Subject), data = df) 
  # m4 <- lmer(NPX ~ Age * Sex + (1 + Age * Sex | Subject), data = df) 
  
  models <- list(m1, m2, m3) #, m4
  
  # Check for singular fit
  singular_flags <- sapply(models, function(m) {
    tryCatch(isSingular(m, tol = 1e-4), error = function(e) NA)
  })
  
  
  tidy(anova(m1,m2,m3,refit=FALSE)) |> #m4
    mutate(Protein = protein,
           Singular = singular_flags) |> 
    relocate(Protein)
  
} 

# Run interaction tests
test_interaction_wellness <-
  map_df(unique(dat_wellness$Assay), ~ test_mm_interaction(dat_wellness, .x))
saveRDS(test_interaction_wellness, savepath_data("MM_test", "test_interaction_wellness.rds"))

test_interaction_bamse <-
  map_df(unique(dat_bamse$Assay), ~ test_mm_interaction(dat_bamse, .x))
saveRDS(test_interaction_bamse, savepath_data("MM_test", "test_interaction_bamse.rds"))


# Run random slope tests
test_slopes_wellness <-
  map_df(unique(dat_wellness$Assay), ~ test_mm_slopes(dat_wellness, .x))
saveRDS(test_slopes_wellness, savepath_data("MM_test", "test_slopes_wellness.rds"))


test_slopes_bamse <-
  map_df(unique(dat_bamse$Assay), ~ test_mm_slopes(dat_bamse, .x))
saveRDS(test_slopes_bamse, savepath_data("MM_test", "test_slopes_bamse.rds"))


# Functions to summarise results (significance and singularity)
summary_signficance <- function(test_df, 
                                test,
                                dataset_name) {
  test_df |>
    mutate(p.value.adj = p.adjust(p.value, method = "BH")) |>
    filter(p.value.adj < 0.05) |> 
    count(term) |> 
    mutate(Dataset = dataset_name,
           Test = test) |> 
    relocate(Test, Dataset)
}

summary_singular <- function(test_df, 
                             test,
                             dataset_name) {
  test_df |>
    count(term, Singular) |> 
    mutate(Dataset = dataset_name,
           Test = test) |> 
    relocate(Test, Dataset) |> 
    filter(Singular == TRUE)
}

# Summary of interaction results
summary_signficance_interaction <- 
  summary_signficance(test_df  = test_interaction_wellness, 
                    test = "Interaction",
                    dataset = "Wellness") |> 
  bind_rows(summary_signficance(test_df = test_interaction_bamse, 
                                test = "Interaction",
                                dataset = "BAMSE"))

write_tsv(summary_signficance_interaction, savepath_data("MM_test", "summary_signficance_interaction.tsv"))


summary_singular_interaction <- 
  summary_singular(test_df  = test_interaction_wellness, 
                 test = "Interaction",
                 dataset = "Wellness") |> 
  bind_rows(summary_singular(test_df = test_interaction_bamse, 
                             test = "Interaction",
                             dataset = "BAMSE"))

write_tsv(summary_singular_interaction, savepath_data("MM_test", "summary_singular_interaction.tsv"))


# Summary of slope results
summary_signficance_slopes <- 
  summary_signficance(test_df  = test_slopes_wellness, 
                    test = "Interaction",
                    dataset = "Wellness") |> 
  bind_rows(summary_signficance(test_df = test_slopes_bamse, 
                                test = "Interaction",
                                dataset = "BAMSE"))

write_tsv(summary_signficance_slopes, savepath_data("MM_test", "summary_signficance_slopes.tsv"))


summary_singular_slopes <- 
  summary_singular(test_df  = test_slopes_wellness, 
                    test = "Interaction",
                    dataset = "Wellness") |> 
  bind_rows(summary_singular(test_df = test_slopes_bamse, 
                             test = "Interaction",
                             dataset = "BAMSE"))

write_tsv(summary_singular_slopes, savepath_data("MM_test", "summary_singular_slopes.tsv"))



