
# Functions for visualization
#library(broom.mixed)
library(lmerTest)
library(performance)
library(effectsize)
library(MatchIt)
library(limma)
library(tidymodels)
library(themis)
library(clusterProfiler)
#library(car)

do_mixed_effect_model <- function(df, protein, type) {
  
  df <- 
    df |> 
    filter(Assay == protein) |> 
    select(DAid, Age, Sex, Subject, NPX) |> 
    mutate(Age = as.numeric(Age),
           Sex = as.factor(Sex),
           Subject = as.factor(Subject))
  
  
  # Fit the mixed-effects model (no random slopes, random intercepts)
  mixed_model <- lmer(NPX ~ Age + Sex + (1 | Subject), data = df) 
  
  # Compute R² for fixed and random effects
  r2_values <- r2_nakagawa(mixed_model)
  
  # Check for singular fit
  if (isSingular(mixed_model, tol = 1e-4)) {
    message(paste("Protein", protein, ": Singular fit detected. Setting random variance explained to 0."))
    marginal_r2 <- r2_values$R2_marginal
    conditional_r2 <- marginal_r2  # No variance from random effects
    random_r2 <- 0  # Explicitly set to 0
  } else {
    marginal_r2 <- r2_values$R2_marginal
    conditional_r2 <- r2_values$R2_conditional
    random_r2 <- conditional_r2 - marginal_r2
  }
  
  # Store variance explained
  variance_explained <- tibble(
    Component = c("Fixed effects (age & sex)", "Random effects (subject)", "Residual"),
    Variance = c(marginal_r2, random_r2, 1 - conditional_r2)
  ) |> mutate(Protein = protein) |> relocate(Protein)
  
  
  # Tidy fixed effects
  fixed_effects <- 
    broom.mixed::tidy(mixed_model, effects = "fixed") |>  
    mutate(Assay = protein)
  
  if(type == "fixed_effects") {
    return(fixed_effects)
  } else if(type == "variance_explained") {
    return(variance_explained)
  }
  
}

# Function to run ANOVA for a given protein
do_anova <- function(df, protein) {
  
  df_protein <- 
    df |> 
    filter(Assay == protein)
  
  # Linear model with all variables
  model <- lm(NPX ~ Age + Sex + BMI + Disease , data = df_protein) 
  
  # Conduct ANOVA
  anova_res <- car::Anova(model, type = 3)
  
  # Calculate Eta-squared using effectsize package
  eta_squared_res <- 
    eta_squared(anova_res, partial = TRUE) |> 
    as.data.frame()
  
  # Get tidy model results
  tidy_anova <- broom::tidy(anova_res)
  
  # Add Eta-squared to the result
  tidy_anova <- 
    tidy_anova |> 
    left_join(eta_squared_res, by = c("term" = "Parameter")) |> 
    mutate(Eta2_perc = Eta2_partial * 100) |> 
    mutate(Protein = protein) |> 
    relocate(Protein) 
  
  return(tidy_anova)
  
}

# do_mixed_effect_model <- function(df, protein) {
#   
#   df <- 
#     df %>%
#     filter(Assay == protein) |> 
#     select(DAid, age, sex, subject, NPX) |> 
#     mutate(age = as.numeric(age))
#  
#   mixed_model <- lmer(NPX ~ age + sex + (1|subject), data = df)
#   #summary(mixed_model)
#   tidy(mixed_model, effects = "fixed") |> 
#     mutate(Assay = protein)
#   
#   # Tidy fixed and random effects
#   tidy(mixed_model) |> 
#     mutate(Protein = protein) |> 
#     relocate(Protein)
#   
# }


# anova(M1,M2) -> compare models

library(performance)



# Function to match case samples to controls
match_samples <-  function(metadata, 
                           case, 
                           control) {
  
  n_males <- 
    metadata |> 
    filter(Disease == case,
           Sex == "M") |> 
    nrow()
  
  n_females <- 
    metadata |> 
    filter(Disease == case,
           Sex == "F") |> 
    nrow()
  
  if(n_males == 0 & n_females > 0) {
    
    metadata <- 
      metadata |> 
      filter(Sex == "F")
    
  } else if (n_males > 0 & n_females == 0) {
    
    metadata <- 
      metadata |> 
      filter(Sex == "M")
    
  } else {
    metadata <- metadata
  }
  
  dat <-
    metadata |>
    filter(Disease %in% c(case, control),
           !is.na(Sex),
           !is.na(Age)) |> 
    mutate(Disease = factor(Disease, levels = c(control, case)))
  
  if(n_males == 0 | n_females == 0) {
     
    set.seed(213)
    output <- matchit(
      Disease ~ Age,
      data = dat,
      method = "nearest",
      distance = "glm"
    )
    
  } else {
    set.seed(213)
    output <- matchit(
      Disease ~ Age + Sex,
      data = dat,
      method = "nearest",
      distance = "glm"
    )
  }
  
  
  match.data(output)
  
}

# Function to run differential expression using limma
do_limma_disease <-
  function(data_wide, 
           metadata,
           disease,
           controls,
           correct = T,
           cutoff = 0) {
    
    # Select current disease
    dat <-
      data_wide %>% 
      inner_join(metadata %>% 
                   select(DAid, Sex, Age, BMI, Disease), by = "DAid") %>% 
      rename(Group = Disease) %>% 
      mutate(Group = ifelse(Group == disease, "1_Case", "0_Control")) 
    
    n_males <- 
      metadata |> 
      filter(Disease == disease,
             Sex == "M") |> 
      nrow()
    
    n_females <- 
      metadata |> 
      filter(Disease == disease,
             Sex == "F") |> 
      nrow()
    
    
    if(correct == T) {
      dat <- 
        dat |> 
        filter(!is.na(Sex),
               !is.na(Age))
    } else {
      
      dat <- dat
    }
    
    # Design a model
    if(correct == T) {
      
      if(n_males == 0 | n_females == 0) {
        design <- model.matrix(~0 + as.factor(dat$Group) + dat$Age) 
        colnames(design) <- c("control", "case", "Age")
      } else if(disease %in% pediatric_diseases & controls == "Healthy") {
        design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex)) 
        colnames(design) <- c("control", "case",  "Sex") 
      } else {
        design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age) 
        colnames(design) <- c("control", "case",  "Sex", "Age") 
      }
      
    } else {
      design <- model.matrix(~0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }
    
    # Make contrast
    contrast <- makeContrasts(Diff = case - control, levels = design)
    
    # Fit linear model to each protein assay
    dat_fit <- 
      dat %>% 
      select(-Sex, -Age, -BMI, -Group)  %>% 
      column_to_rownames("DAid") %>% 
      t()
    
    fit <- lmFit(dat_fit, design = design, maxit = 100000) #method = "robust", 
    
    # Apply contrast
    contrast_fit <- contrasts.fit(fit, contrast)
    
    # Apply empirical Bayes smoothing to the SE
    ebays_fit <- eBayes(contrast_fit, robust = T)
    
    # Extract DE results
    DE_results <-
      topTable(ebays_fit,
               n = nrow(ebays_fit$p.value), 
               adjust.method = "fdr", 
               confint = TRUE)
    
    DE_res <- 
      DE_results %>% 
      as_tibble(rownames = "Assay") %>% 
      mutate(Disease = disease,
             sig = case_when(adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
                             adj.P.Val < 0.05 & logFC > cutoff ~ "significant up", 
                             T ~ "not significant"),
             Control = controls)
    
    return(DE_res)
  }


generate_split <- function(data, 
                           proportion = 0.7,
                           seed = 213,
                           variable_stratify) {
  
  set.seed(seed)
  
  data_split <-
    data |> 
    initial_split(prop = proportion, strata = !!sym(variable_stratify))
  
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  return(list("data_split" = data_split, 
              "data_train" = data_train,
              "data_test" = data_test))
}

# Function to run lasso pipeline for prediction of discrete variables
discrete_prediction <-  
  function(variable_predict,
           variable_case,
           split_train, 
           split_test,
           seed) {
    
    cat(paste0("\nPreparing training and testing data for ", variable_predict))
    
    split_train <- 
      split_train |> 
      rename(Variable = !!sym(variable_predict)) |> 
      mutate(Variable = case_when(Variable == variable_case ~ paste0("1_", Variable),
                                  Variable != variable_case ~ paste0("0_Control"))) |> 
      mutate(Variable = factor(Variable))
    
    split_test <- 
      split_test |> 
      rename(Variable = !!sym(variable_predict)) |> 
      mutate(Variable = case_when(Variable == variable_case ~ paste0("1_", Variable),
                                  Variable != variable_case ~ paste0("0_Control"))) |> 
      mutate(Variable = factor(Variable))
    
    variable_split <- make_splits(split_train, 
                                  split_test)
    
    
    cat(paste0("\nDefining ML specs for ", variable_predict))
    
    # Recipe with ML steps
    discrete_recipe <- 
      recipe(Variable ~ ., data = split_train) |> 
      step_relevel(Variable, ref_level = "0_Control") |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
      # step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric())
    
    # LASSO model specifications
    glmnet_specs <- 
      logistic_reg() |> 
      set_mode("classification") |> 
      set_engine("glmnet") |> 
      set_args(penalty = tune(), 
               mixture = 1) 
    
    # ML workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(discrete_recipe) |> 
      add_model(glmnet_specs) 
    
    # Define glmnet grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 20)
    
    # Define the resamples (CV)
    set.seed(213)
    ml_rs <- vfold_cv(split_train, v = 10, strata = Variable)
    
    # Define the evaluation metrics (add brier)
    eval_metrics <- metric_set(roc_auc)
    
    # Define control_grid
    set.seed(213)
    ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything", event_level = "second") 
    
    cat(paste0("\nFitting glmnet model for ", variable_predict))
    
    # Glmnet grid search
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = ml_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics
      )
    
    cat(paste0("\nSelecting best performing model for ", variable_predict))
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    # Select best hyperparameter
    best_glmnet <- 
      select_best(glmnet_res, metric = "roc_auc") |> 
      select(-.config)
    
    #Finalize the workflow and fit the final model
    glmnet_wflow <- 
      glmnet_wflow |>  
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- last_fit(glmnet_wflow, variable_split, metrics = eval_metrics) 
    
    # Extract model performance
    performance <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      select(-.config, -.estimator)
    
    
    glmnet_auc <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      filter(.metric == "roc_auc") |> 
      pull(.estimate) |> 
      round(2)
    
    # Extract protein importance
    important_proteins <- 
      final_glmnet_fit |> 
      extract_fit_parsnip()  |> 
      vip::vi(lambda = best_glmnet$penalty, event_level = "second")  |> 
      mutate(
        Importance = abs(Importance),
        Variable = fct_reorder(Variable, Importance)
      )
    
    # Extract model predictions
    predictions <- 
      final_glmnet_fit |> 
      collect_predictions(summarize = F) 
    
    # Add performance
    # ROC curve
    roc <- 
      predictions |>
      roc_curve(truth = Variable, paste0(".pred_1_", variable_case))
    
    # AUC
    res <- 
      predictions |> 
      rename(prediction = paste0(".pred_1_", variable_case)) 
    
    auc <- pROC::auc(res$Variable, res$prediction)
    ci <- pROC::ci.auc(res$Variable, res$prediction) 
    
    # Combine
    combined_roc_auc <- 
      roc |> 
      mutate(AUC = as.numeric(auc),
             CI_lower = as.numeric(ci)[[1]],
             CI_upper = as.numeric(ci)[[3]])
    

    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_workflow" = glmnet_wflow,
                "final_fit" = final_glmnet_fit,
                "predictions" = predictions,
                "performance" = performance,
                "auc" = auc,
                "auc_ci" = ci,
                "roc_curve" = roc, 
                "important_proteins" = important_proteins))
  }


#Function to run lasso pipeline for prediction of continuous variables
continuous_prediction <-  
  function(split_train,
           split_test,
           variable_predict,
           seed#, 
           #path_save
           ) {
    
    
    cat(paste0("\nPreparing training and testing data for ", variable_predict))
    
    split_train <- 
      split_train |> 
      rename(Variable = !!sym(variable_predict))
    
    split_test <- 
      split_test |> 
      rename(Variable = !!sym(variable_predict))
    
    variable_split <- make_splits(split_train, 
                                  split_test)
    
    cat(paste0("\nDefining ML specs for ", variable_predict))
    
    # Recipe with ML steps
    ml_recipe <- 
      recipe(Variable ~ ., data = split_train) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric_predictors()) |> 
      step_nzv(all_numeric_predictors()) |> 
      #step_corr(all_numeric_predictors()) |> 
      step_impute_knn(all_numeric_predictors())
    
    # LASSO model specifications
    glmnet_specs <- 
      linear_reg() |> 
      set_mode("regression") |> 
      set_engine("glmnet") |> 
      set_args(penalty = tune(), 
               mixture = 1) 
    
    # ML workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_specs) 
    
    # Define glmnet grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 20)
    
    # Define the resamples (CV)
    set.seed(213)
    ml_rs <- vfold_cv(split_train, v = 10, strata = Variable)
    
    
    # Define control_grid
    set.seed(213)
    ctrl <- control_grid(save_pred = TRUE, parallel_over = "everything") 
    
    cat(paste0("\nFitting glmnet model for ", variable_predict))
    
    # Glmnet grid search
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = ml_rs,
        grid = glmnet_grid,
        control = ctrl
      )
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    cat(paste0("\nSelecting best performing model for ", variable_predict))
    
    # Select best hyperparameter
    best_glmnet <- 
      select_best(glmnet_res, metric = "rmse") |> 
      select(-.config)
    
    #Finalize the workflow and fit the final model
    glmnet_wflow <- 
      glmnet_wflow |>  
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- last_fit(glmnet_wflow, variable_split) 
    
    # Extract model performance
    performance <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      select(-.config, -.estimator)
    
    glmnet_metrics <- 
      final_glmnet_fit |> 
      collect_metrics() |> 
      filter(.metric == "rmse") 
    
    # Extract protein importance
    important_proteins <- 
      final_glmnet_fit |> 
      extract_fit_parsnip()  |> 
      vip::vi(lambda = best_glmnet$penalty, event_level = "second")  |> 
      mutate(
        Importance = abs(Importance),
        Variable = fct_reorder(Variable, Importance)
      )
    
    #write_tsv(important_proteins, paste(path_save, "/important_proteins_", variable_predict, "_seed_", seed, ".tsv", sep = ""))
    
    # Extract model predictions
    predictions <- 
      final_glmnet_fit |> 
      collect_predictions(summarize = F) 
    
    # Add performance
    #write_tsv(predictions, paste(path_save, "/predictions_", variable_predict, "_seed_", seed, ".tsv", sep = ""))
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_workflow" = glmnet_wflow,
                "final_fit" = final_glmnet_fit,
                "predictions" = predictions,
                "performance" = performance,
                "important_proteins" = important_proteins))
  }

biological_themes <- list(
  Adaptive_Immunity = c("adaptive", "T cell", "B cell", "antigen", "MHC", "lymphocyte", "immunoglobulin"),
  Innate_Immunity = c("innate", "macrophage", "neutrophil", "complement", "inflammatory", "pattern recognition", "TLR"),
  Cytokine_Signaling = c("cytokine", "interleukin", "interferon", "chemokine", "TNF"),
  Cell_Cycle = c("cell cycle", "mitosis", "proliferation", "DNA replication"),
  Apoptosis_Stress = c("apoptosis", "oxidative", "hypoxia", "stress", "cell death"),
  Lipid_Metabolism = c("lipid", "fatty acid", "cholesterol", "triglyceride"),
  Glucose_Insulin = c("glucose", "insulin", "glycolysis", "gluconeogenesis"),
  ECM_Adhesion = c("extracellular matrix", "collagen", "adhesion", "integrin", "fibrosis"),
  Angiogenesis = c("angiogenesis", "vascular", "blood vessel", "endothelial", "VEGF"),
  Hemostasis = c("coagulation", "platelet", "fibrin", "thrombosis"),
  Neuro = c("neuron", "neurotransmitter", "synapse", "axon"),
  Wound_Healing = c("wound", "repair", "regeneration"),
  Development = c("development", "morphogenesis", "embryo", "differentiation"),
  Metabolism_General = c("metabolic", "metabolism"),
  
  # Additional themes
  DNA_Repair = c("DNA repair", "genome integrity", "mismatch repair", "nucleotide excision", "double strand break"),
  Signal_Transduction = c("kinase", "phosphorylation", "MAPK", "PI3K", "signaling pathway"),
  Immune_Checkpoints = c("PD-1", "CTLA-4", "immune checkpoint", "regulatory T cell"),
  Autophagy = c("autophagy", "lysosome", "degradation"),
  Mitochondrial_Function = c("mitochondria", "oxidative phosphorylation", "respiratory chain"),
  Protein_Folding = c("chaperone", "folding", "heat shock protein", "proteostasis"),
  Epigenetic_Regulation = c("methylation", "histone", "chromatin", "epigenetic"),
  Extracellular_Vesicles = c("exosome", "extracellular vesicle", "microvesicle"),
  Fibrosis_Remodeling = c("fibrosis", "scar", "remodeling", "TGF-beta"),
  Hormone_Signaling = c("hormone", "estrogen", "androgen", "glucocorticoid")
)

df_shared_themes <- data.frame(
  Theme = c(
    "Adaptive_Immunity", "Innate_Immunity", "Cytokine_Signaling", "Immune_Checkpoints",
    "Cell_Cycle", "Apoptosis_Stress", "DNA_Repair", "Signal_Transduction",
    "Autophagy", "Protein_Folding", "Epigenetic_Regulation",
    "Lipid_Metabolism", "Glucose_Insulin", "Metabolism_General", "Mitochondrial_Function",
    "ECM_Adhesion", "Fibrosis_Remodeling", "Extracellular_Vesicles",
    "Angiogenesis", "Hemostasis", "Neuro", "Wound_Healing", "Development", "Hormone_Signaling"
  ),
  Shared_Broad_Category = c(
    "Immunity", "Immunity", "Immunity", "Immunity",
    "Cell Regulation", "Cell Regulation", "Cell Regulation", "Cell Regulation",
    "Cell Regulation", "Cell Regulation", "Cell Regulation",
    "Metabolism", "Metabolism", "Metabolism", "Metabolism",
    "Structural / ECM", "Structural / ECM", "Structural / ECM",
    "Vascular", "Vascular", "Nervous System", "Repair / Regeneration", "Developmental Biology", "Endocrine System"
  ),
  stringsAsFactors = FALSE
)

assign_theme <- function(description) {
  for (theme in names(biological_themes)) {
    keywords <- biological_themes[[theme]]
    if (any(str_detect(tolower(description), keywords))) {
      return(theme)
    }
  }
  return("Other")
}