
# Functions for visualization
library(limma)
library(tidymodels)
library(themis)
library(effectsize)
library(lme4)
library(multilevelmod)
library(broom.mixed)
library(lmerTest)
library(igraph)
library(ggraph)
library(FNN)
library(MatchIt)


# Function to run ANOVA for a given protein
do_anova <- function(df, protein) {
  
  df <- 
    df |> 
    filter(Assay == protein)
  
  # Linear model with all variables
  model <- lm(NPX ~ subject + age + sex + visit , data = df) #+ month_of_sampling
  
  # Conduct ANOVA
  anova_res <- anova(model)
  
  # Calculate Eta-squared using effectsize package
  eta_squared_res <- 
    eta_squared(model, partial = TRUE) |> 
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


library(performance)

do_mixed_effect_model <- function(df, protein, type) {
  
  df <- 
    df %>%
    filter(Assay == protein) %>% 
    select(DAid, age, sex, subject, NPX) %>% 
    mutate(age = as.numeric(age))
  
  # Fit the mixed-effects model
  mixed_model <- lmer(NPX ~ age + sex + (1 | subject), data = df)
  
  # Extract variance components
  var_comp <- as.data.frame(VarCorr(mixed_model))
  random_effect_variance <- var_comp$vcov[var_comp$grp == "subject"]
  residual_variance <- sigma(mixed_model)^2
  total_variance <- sum(var_comp$vcov) + residual_variance
  
  # Marginal and conditional R²
  r2_values <- r2(mixed_model)
  marginal_r2 <- r2_values$R2_marginal  # Variance explained by fixed effects
  conditional_r2 <- r2_values$R2_conditional  # Variance explained by fixed + random effects
  
  # Variance explained by each component
  fixed_effect_variance <- marginal_r2 * total_variance
  random_effect_explained <- (conditional_r2 - marginal_r2) * total_variance
  residual_unexplained <- residual_variance / total_variance
  
  variance_explained <- data.frame(
    Component = c("Fixed Effects", "Random Effects (Subject)", "Residual"),
    Variance = c(fixed_effect_variance, random_effect_explained, residual_variance),
    Proportion = c(marginal_r2, conditional_r2 - marginal_r2, residual_unexplained)
  ) %>% 
    as_tibble() |> 
    mutate(Assay = protein)
  
  # Tidy fixed effects
  fixed_effects <- tidy(mixed_model, effects = "fixed") %>% 
    mutate(Assay = protein)
  
  if(type == "fixed_effects") {
    return(fixed_effects)
  } else if(type == "variance_explained") {
    return(variance_explained)
  }
  
}


# Function to match case samples to controls
match_samples <-  function(metadata, 
                           case, 
                           control) {
  
  if(case %in% female_diseases) {
    
    metadata <- 
      metadata |> 
      filter(Sex == "F")
    
  } else if (case %in% male_diseases) {
    
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
  
  if(case %in% c(male_diseases, female_diseases)) {
    
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
      
      if(disease %in% c(male_diseases, female_diseases)) {
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


# 
# do_limma <-
#   function(data, 
#            metadata,
#            case,
#            control,
#            correct = T,
#            cutoff = 0) {
#     
#     # Widen data
#     data_wide <- 
#       data |> 
#       select(DAid, Assay, NPX) |> 
#       pivot_wider(values_from = NPX, names_from = Assay)
# 
#     # Select current disease
#     dat <-
#       data_wide %>% 
#       inner_join(metadata %>% 
#                    select(DAid, Sex, Diagnose), by = "DAid") %>% 
#       rename(Group = Diagnose) %>% 
#       filter(Group %in% c(case, control)) |> 
#       mutate(Group = ifelse(Group == case, "1_Case", "0_Control")) 
#     
#     
#     if(correct == T) {
#       dat <- 
#         dat |> 
#         filter(!is.na(Sex))
#     } else {
#       
#       dat <- dat
#     }
#     
#     
#     # Design a model - add Group, and Sex, Age, BMI
#     if(correct == T ) {
#       design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex)) 
#       colnames(design) <- c("control", "case",  "Sex") 
#     } else {
#       design <- model.matrix(~0 + as.factor(dat$Group))
#       colnames(design) <- c("control", "case")
#     }
#     
#     # Make contrast
#     contrast <- makeContrasts(Diff = case - control, levels = design)
#     
#     # Fit linear model to each protein assay
#     dat_fit <- 
#       dat %>% 
#       select(-Sex, -Group)  %>% 
#       column_to_rownames("DAid") %>% 
#       t()
#     
#     fit <- lmFit(dat_fit, design = design,  method = "robust", maxit = 10000)
#     
#     # Apply contrast
#     contrast_fit <- contrasts.fit(fit, contrast)
#     
#     # Apply empirical Bayes smoothing to the SE
#     ebays_fit <- eBayes(contrast_fit)
#     
#     # Extract DE results
#     DE_results <-
#       topTable(ebays_fit,
#                n = nrow(ebays_fit$p.value), 
#                adjust.method = "fdr", 
#                confint = TRUE)
#     
#     DE_res <- 
#       DE_results %>% 
#       as_tibble(rownames = "Assay") %>% 
#       mutate(Comparison = paste0(case, " vs ", control),
#              sig = case_when(adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
#                              adj.P.Val < 0.05 & logFC > cutoff ~ "significant up", 
#                              T ~ "not significant")) 
#     
#     return(DE_res)
#   }
# 
# 
# 
# 
# 
# # Function to run differential expression using limma
# do_limma_disease <-
#   function(data_wide, 
#            metadata,
#            disease,
#            controls,
#            correct = T,
#            cutoff = 0) {
#     
#     # Select current disease
#     dat <-
#       data_wide %>% 
#       inner_join(metadata %>% 
#                    select(DAid, Sex, Age, BMI, Disease), by = "DAid") %>% 
#       rename(Group = Disease) %>% 
#       mutate(Group = ifelse(Group == disease, "1_Case", "0_Control")) 
#     
#     
#     if(correct == T) {
#       dat <- 
#         dat |> 
#         filter(!is.na(Sex),
#                !is.na(Age))
#     } else {
#       
#       dat <- dat
#     }
#     
#     # Design a model
#     if(correct == T) {
#       
#       if(disease %in% c(male_diseases, female_diseases)) {
#         design <- model.matrix(~0 + as.factor(dat$Group) + dat$Age) 
#         colnames(design) <- c("control", "case", "Age")
#       } else if(disease %in% pediatric_diseases & controls == "Healthy") {
#         design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex)) 
#         colnames(design) <- c("control", "case",  "Sex") 
#       } else {
#         design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age) 
#         colnames(design) <- c("control", "case",  "Sex", "Age") 
#       }
#       
#     } else {
#       design <- model.matrix(~0 + as.factor(dat$Group))
#       colnames(design) <- c("control", "case")
#     }
#     
#     # Make contrast
#     contrast <- makeContrasts(Diff = case - control, levels = design)
#     
#     # Fit linear model to each protein assay
#     dat_fit <- 
#       dat %>% 
#       select(-Sex, -Age, -BMI, -Group)  %>% 
#       column_to_rownames("DAid") %>% 
#       t()
#     
#     fit <- lmFit(dat_fit, design = design, maxit = 100000) #method = "robust", 
#     
#     # Apply contrast
#     contrast_fit <- contrasts.fit(fit, contrast)
#     
#     # Apply empirical Bayes smoothing to the SE
#     ebays_fit <- eBayes(contrast_fit, robust = T)
#     
#     # Extract DE results
#     DE_results <-
#       topTable(ebays_fit,
#                n = nrow(ebays_fit$p.value), 
#                adjust.method = "fdr", 
#                confint = TRUE)
#     
#     DE_res <- 
#       DE_results %>% 
#       as_tibble(rownames = "Assay") %>% 
#       mutate(Disease = disease,
#              sig = case_when(adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
#                              adj.P.Val < 0.05 & logFC > cutoff ~ "significant up", 
#                              T ~ "not significant"),
#              Control = controls)
#     
#     return(DE_res)
#   }
# 
# do_limma_continuous <-
#   function(data, 
#            metadata,
#            case,
#            control,
#            correct = T,
#            cutoff = 0) {
#     
#     # Widen data
#     data_wide <- 
#       data |> 
#       select(DAid, Assay, NPX) |> 
#       pivot_wider(values_from = NPX, names_from = Assay)
#     
#     # Select current disease
#     dat <-
#       data_wide %>% 
#       inner_join(metadata %>% 
#                    select(DAid, Sex, Diagnose), by = "DAid") %>% 
#       rename(Group = Diagnose) %>% 
#       filter(Group %in% c(case, control)) |> 
#       mutate(Group = ifelse(Group == case, "1_Case", "0_Control")) 
#     
#     
#     if(correct == T) {
#       dat <- 
#         dat |> 
#         filter(!is.na(Sex))
#     } else {
#       
#       dat <- dat
#     }
#     
#     
#     # Design a model - add Group, and Sex, Age, BMI
#     if(correct == T ) {
#       design <- model.matrix(~Difference + as.factor(dat$Sex)) 
#       colnames(design) <- c("control", "case",  "Sex") 
#     } else {
#       design <- model.matrix(~Difference)
#       colnames(design) <- c("control", "case")
#     }
#     
# 
#     # Fit linear model to each protein assay
#     dat_fit <- 
#       dat %>% 
#       select(-Sex, -Group)  %>% 
#       column_to_rownames("DAid") %>% 
#       t()
#     
#     fit <- lmFit(dat_fit, design = design,  method = "robust", maxit = 10000)
#     
#     
#     
#     # Apply empirical Bayes smoothing to the SE
#     ebays_fit <- eBayes(fit)
#     
#     # Extract DE results
#     de_results <- topTable(ebays_fit,
#                                   coef = variable,  # Extract results for Age specifically
#                                   n = nrow(ebays_fit$p.value),
#                                   adjust.method = "fdr",
#                                   confint = TRUE)
#     
#     DE_res <- 
#       DE_results %>% 
#       as_tibble(rownames = "Assay") %>% 
#       mutate(Comparison = paste0(case, " vs ", control),
#              sig = case_when(adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
#                              adj.P.Val < 0.05 & logFC > cutoff ~ "significant up", 
#                              T ~ "not significant")) 
#     
#     return(DE_res)
#   }


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
# 
# do_lasso_binary <-  
#   function(disease, 
#            metadata,
#            split_train, 
#            split_test,
#            seed = 213) {
#     
#     # Prepare data - make custom split for current disease 
#     cat(paste0("\nPreparing training and testing data for ", disease))
#     
#     training_dat <- 
#       split_train |> 
#       mutate(Disease = case_when(Disease == disease ~ paste0("1_", disease),
#                                  T ~ "0_Other")) |> 
#       mutate(Disease = factor(Disease))
#     
#     testing_dat <- 
#       split_test |> 
#       mutate(Disease = case_when(Disease == disease ~ paste0("1_", disease),
#                                  T ~ "0_Other"))|> 
#       mutate(Disease = factor(Disease))
#     
#     
#     # For male and female diseases, filter only male and female samples respectively 
#     if(disease %in% female_diseases) {
#       
#       training_dat <- 
#         training_dat |> 
#         left_join(metadata |> 
#                     select(DAid, Sex), by = "DAid") |> 
#         filter(Sex == "F") |> 
#         select(-Sex)
#       
#       testing_dat <- 
#         testing_dat |> 
#         left_join(metadata |> 
#                     select(DAid, Sex), by = "DAid") |> 
#         filter(Sex == "F") |> 
#         select(-Sex)
#       
#     } else if (disease %in% male_diseases) {
#       
#       training_dat <- 
#         training_dat |> 
#         left_join(metadata |> 
#                     select(DAid, Sex), by = "DAid") |> 
#         filter(Sex == "M") |> 
#         select(-Sex)
#       
#       testing_dat <- 
#         testing_dat |> 
#         left_join(metadata |> 
#                     select(DAid, Sex), by = "DAid") |> 
#         filter(Sex == "M") |> 
#         select(-Sex)
#       
#     } else {
#       
#       training_dat <- training_dat
#       testing_dat <- testing_dat
#       
#     }
#     
#     disease_split <- make_splits(training_dat, testing_dat)
#     
#     cat(paste0("\nDefining ML specs for ", disease))
#     
#     # Recipe with ML steps
#     disease_recipe <- 
#       recipe(Disease ~ ., data = training_dat) |> 
#       step_relevel(Disease, ref_level = "0_Other") |> 
#       update_role(DAid, new_role = "id") |> 
#       step_normalize(all_numeric()) |> 
#       step_nzv(all_numeric()) |> 
#       step_impute_knn(all_numeric()) |>   
#       step_downsample(Disease)
#     
#     # LASSO model specifications
#     glmnet_specs <- 
#       logistic_reg() |> 
#       set_mode("classification") |> 
#       set_engine("glmnet") |> 
#       set_args(penalty = tune(), 
#                mixture = 1) 
#     
#     # ML workflow
#     glmnet_wflow <-
#       workflow() |> 
#       add_recipe(disease_recipe) |> 
#       add_model(glmnet_specs) 
#     
#     # Define glmnet grid
#     set.seed(213)
#     glmnet_grid <-
#       glmnet_wflow |>
#       extract_parameter_set_dials() |>
#       grid_latin_hypercube(size = 20)
#     
#     # Define the resamples (CV)
#     set.seed(213)
#     cancer_rs <- vfold_cv(training_dat, v = 10, strata = Disease)
#     
#     # Define the evaluation metrics (add brier)
#     eval_metrics <- metric_set(roc_auc)
#     
#     # Define control_grid
#     set.seed(213)
#     ctrl <- control_grid(save_pred = TRUE, 
#                          event_level = "second",
#                          parallel_over = "everything")
#     
#     
#     cat(paste0("\nFitting glmnet model for ", disease))
#     
#     # Glmnet grid search
#     set.seed(213)
#     glmnet_res <-
#       glmnet_wflow |>
#       tune_grid(
#         resamples = cancer_rs,
#         grid = glmnet_grid,
#         control = ctrl,
#         metrics = eval_metrics
#       )
#     
#     cat(paste0("\nSelecting best performing model for ", disease))
#     
#     predictions_train <- 
#       glmnet_res |> 
#       collect_predictions()
#     
#     metrics_train <- 
#       glmnet_res |> 
#       collect_metrics()
#     
#     # Select best hyperparameter
#     best_glmnet <- 
#       select_best(glmnet_res, metric = "roc_auc") |> 
#       select(-.config)
#     
#     #Finalize the workflow and fit the final model
#     glmnet_wflow <- 
#       glmnet_wflow |>  
#       finalize_workflow(best_glmnet)
#     
#     final_glmnet_fit <- last_fit(glmnet_wflow, disease_split, metrics = eval_metrics) 
#     
#     # Extract model performance
#     performance <- 
#       final_glmnet_fit |> 
#       collect_metrics() |> 
#       select(-.config, -.estimator)
#     
#     glmnet_auc <- 
#       final_glmnet_fit |> 
#       collect_metrics() |> 
#       filter(.metric == "roc_auc") |> 
#       pull(.estimate) |> 
#       round(2)
#     
#     # Extract protein importance
#     important_proteins <- 
#       final_glmnet_fit |> 
#       extract_fit_parsnip()  |> 
#       vip::vi(lambda = best_glmnet$penalty, event_level = "second")  |> 
#       mutate(
#         Importance = abs(Importance),
#         Variable = fct_reorder(Variable, Importance)
#       )
#     
#     # Extract model predictions
#     predictions <- 
#       final_glmnet_fit |> 
#       collect_predictions(summarize = F)  |> 
#       mutate(DAid = testing_dat$DAid) |> 
#       relocate(DAid)
#     
#     # ROC curve
#     roc <- 
#       predictions |>
#       roc_curve(truth = Disease, paste0(".pred_1_", disease), event_level = "second")
#     
#     # AUC
#     res <- 
#       predictions |> 
#       rename(prediction = paste0(".pred_1_", disease)) |> 
#       mutate(Disease = factor(Disease, levels = c("0_Other", paste0("1_", disease))))
#     
#     auc <- pROC::auc(res$Disease, res$prediction)
#     ci <- pROC::ci.auc(res$Disease, res$prediction) 
#     
#     
#     # Combine
#     combined_roc_auc <- 
#       roc |> 
#       mutate(AUC = as.numeric(auc),
#              CI_lower = as.numeric(ci)[[1]],
#              CI_upper = as.numeric(ci)[[3]])
#     
#     # Collect all data and save as rds
#     complete_model <- 
#       list("penalty" = best_glmnet,
#            "glmnet_model" = glmnet_res,
#            "predictions_train" = predictions_train, 
#            "performance_train" = metrics_train,
#            "final_workflow" = glmnet_wflow,
#            "final_fit" = final_glmnet_fit,
#            "predictions" = predictions,
#            "performance" = performance,
#            "auc" = auc,
#            "auc_ci" = ci,
#            "roc_curve" = roc, 
#            "important_proteins" = important_proteins)
#     
#     return(complete_model)
#   }
# 
# do_limma_continuous <- function(join_data,
#                                 variable,
#                                 correct = c("Sex"),
#                                 correct_type = c("factor"),
#                                 pval_lim = 0.05,
#                                 logfc_lim = 0) {
#   nrows_before <- nrow(join_data)
#   
#   join_data <- join_data |>
#     dplyr::filter(!dplyr::if_any(dplyr::all_of(c(variable, correct)), is.na))  # Remove NAs from columns in formula
#   
#   nrows_after <- nrow(join_data)
#   if (nrows_before != nrows_after){
#     warning(paste0(nrows_before - nrows_after,
#                    " rows were removed because they contain NAs in ",
#                    variable,
#                    " or ",
#                    paste(correct, collapse = ", "),
#                    "!"))
#   }
#   
#   # Design a model
#   # formula <- paste("~0 +" , variable)
#   
#   # if (!is.null(correct)) {
#   #   for (i in 1:length(correct)) {
#   #     if (correct_type[i] == "factor") {
#   #       cofactor = paste("as.factor(", correct[i], ")")
#   #     } else {
#   #       cofactor = correct[i]
#   #     }
#   #     formula <- paste(formula, "+", cofactor)
#   #   }
#   # }
#   #design <- stats::model.matrix(stats::as.formula(formula), data = join_data)
#   
#   # Design a model
#   formula <- paste("~" , variable)  # Include intercept
#   
#   if (!is.null(correct)) {
#     for (i in 1:length(correct)) {
#       if (correct_type[i] == "factor") {
#         cofactor = paste("as.factor(", correct[i], ")")
#       } else {
#         cofactor = correct[i]
#       }
#       formula <- paste(formula, "+", cofactor)
#     }
#   }
#   design <- stats::model.matrix(stats::as.formula(formula), data = join_data)
#   
#   
#   # Fit linear model to each protein assay
#   data_fit <- 
#     join_data |>
#     dplyr::select(-dplyr::any_of(c(variable, correct))) |>
#     tibble::column_to_rownames("DAid") |>
#     t()
#   
#   fit <- limma::lmFit(data_fit, design = design, method = "robust", maxit = 10000)
#   
#   
#   # Apply empirical Bayes smoothing to the SE
#   ebays_fit <- limma::eBayes(fit)
#   
#   # Extract DE results
#   de_results <- limma::topTable(ebays_fit,
#                                 coef = variable,  # Extract results for Age specifically
#                                 n = nrow(ebays_fit$p.value),
#                                 adjust.method = "fdr",
#                                 confint = TRUE)
#   
#   de_res <- de_results |>
#     tibble::as_tibble(rownames = "Assay") |>
#     dplyr::rename(logFC = colnames(de_results)[1]) |>
#     dplyr::mutate(sig = dplyr::case_when(
#       adj.P.Val < pval_lim & logFC < -logfc_lim ~ "significant down",
#       adj.P.Val < pval_lim & logFC > logfc_lim ~ "significant up",
#       T ~ "not significant")
#     ) |>
#     dplyr::arrange(adj.P.Val)
#   
#   return(de_res)
# }
# 
# 
