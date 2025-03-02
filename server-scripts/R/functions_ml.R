#### Title: Sever functions
#### Author: María Bueno Álvez
#### Description: script collecting general functions to be used when running pipelines on the server
#### Last edited : 12/08/2024

# Levels
female_diseases <- c("Breast cancer", "Breast ductal carcinoma in situ", "Cervical cancer", "Endometrial cancer", "Ovarian cancer")
male_diseases <-  c("Prostate cancer", "Abdominal aortic aneurysm")

#Palette
library(RColorBrewer)
getPalette3 = colorRampPalette(brewer.pal(8, "Set2"))
pal_class<-getPalette3(8)
names(pal_class)<-c("Psychiatric","Cardiovascular","Cancer","Autoimmune","Pediatric","Infection","Metabolic","Healthy") 
class_order <- c("Healthy", "Cardiovascular","Metabolic","Cancer","Psychiatric","Autoimmune","Infection","Pediatric")

# HPA theme
theme_hpa <- 
  function(angled = F, axis_x = T, axis_y = T, facet_title = T) {
    t <- 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = rel(1), hjust = 0.5),
        plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.8)),
        legend.key.size= unit(0.7, "cm"),
        legend.title = element_text(size=rel(1)),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="grey90",fill="grey90"),
        strip.text = element_text(face="bold")
      )
    
    if(angled) {
      t <- t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
    }
    
    if(axis_x == F) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())
    } 
    
    if(axis_y == F) {
      t <- t +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank())
    }
    if(facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }
# Function to generate splits
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

# Function to run binary classification against healthy
disease_against_healthy <-  
  function(disease, 
           split_train, 
           split_test,
           seed, 
           path_save) {
    
    # Prepare data - make custom split for current disease (+ Healthy samples) based on master split (including all diseases)

cat(paste0("\nPreparing training and testing data for ", disease))

    training_dat <- 
      split_train |> 
      filter(Disease %in% c(disease, "Healthy")) |> 
      mutate(Disease = case_when(Disease == disease ~ paste0("1_", disease),
                                 Disease == "Healthy" ~ "0_Healthy")) |> 
      mutate(Disease = factor(Disease))
    
    testing_dat <- 
      split_test |> 
      filter(Disease %in% c(disease, "Healthy")) |> 
      mutate(Disease = case_when(Disease == disease ~ paste0("1_", disease),
                                 Disease == "Healthy" ~ "0_Healthy"))|> 
      mutate(Disease = factor(Disease))
    
    
    # For male and female diseases, filter only male and female samples respectively 
    if(disease %in% female_diseases) {
      
      training_dat <- 
        training_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "F") |> 
        select(-Sex)
      
      testing_dat <- 
        testing_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "F") |> 
        select(-Sex)
    } else if (disease %in% male_diseases) {
      
      training_dat <- 
        training_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "M") |> 
        select(-Sex)
      
      testing_dat <- 
        testing_dat |> 
        left_join(resource_meta |> 
                    select(DAid, Sex), by = "DAid") |> 
        filter(Sex == "M") |> 
        select(-Sex)
    } else {
      training_dat <- training_dat
      testing_dat <- testing_dat
      
    }
    
    disease_split <- make_splits(training_dat, testing_dat)
    
cat(paste0("\nDefining ML specs for ", disease))

    # Recipe with ML steps
    disease_recipe <- 
      recipe(Disease ~ ., data = training_dat) |> 
      step_relevel(Disease, ref_level = "0_Healthy") |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
      #step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric()) |>   
      step_downsample(Disease)
    
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
      add_recipe(disease_recipe) |> 
      add_model(glmnet_specs) 
    
    # Define glmnet grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 20)
    
    # Define the resamples (CV)
    set.seed(213)
    cancer_rs <- vfold_cv(training_dat, v = 10, strata = Disease)
    
    # Define the evaluation metrics (add brier)
    eval_metrics <- metric_set(roc_auc)
    
    # Define control_grid
    set.seed(213)
    ctrl <- control_grid(save_pred = TRUE, 
                         event_level = "second",
                         parallel_over = "everything") # look at extract = identity
    

cat(paste0("\nFitting glmnet model for ", disease))

    # Glmnet grid search
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = cancer_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics
      )
    #autoplot(glmnet_res)
    
cat(paste0("\nSelecting best performing model for ", disease))

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
    
    final_glmnet_fit <- last_fit(glmnet_wflow, disease_split, metrics = eval_metrics) 
    
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
    
    write_tsv(important_proteins, paste(path_save, "/important_proteins_", disease, "_healthy_seed_", seed, ".tsv", sep = ""))
   
    # Extract model predictions
    predictions <- 
      final_glmnet_fit |> 
      collect_predictions(summarize = F)  |> 
      mutate(DAid = testing_dat$DAid) |> 
      relocate(DAid)
    
    write_tsv(predictions, paste(path_save, "/predictions_", disease, "_healthy_seed_", seed, ".tsv", sep = ""))
    
    # Confusion matrix
    cm <-
      predictions |>
      mutate(pred = ifelse(.pred_0_Healthy > 0.5, "Healthy", disease),
             Disease = ifelse(Disease == "0_Healthy", "Healthy", disease)) |>
      mutate(Disease = factor(Disease, levels = c(disease, "Healthy")),
             pred = factor(pred, levels = c(disease, "Healthy"))) |> 
      conf_mat(Disease, pred)
    
    # ROC curve
    roc <- 
      predictions |>
      roc_curve(truth = Disease, paste0(".pred_1_", disease), event_level = "second")

# AUC
    res <- 
      predictions |> 
      rename(prediction = paste0(".pred_1_", disease)) |> 
      mutate(Disease = factor(Disease, levels = c("0_Healthy", paste0("1_", disease))))
    
    auc <- pROC::auc(res$Disease, res$prediction)
    ci <- pROC::ci.auc(res$Disease, res$prediction) 
    
   
    # Combine
    combined_roc_auc <- 
      roc |> 
      mutate(AUC = as.numeric(auc),
             CI_lower = as.numeric(ci)[[1]],
             CI_upper = as.numeric(ci)[[3]])
     
    write_tsv(combined_roc_auc, paste(path_save, "/roc_", disease, "_healthy_seed_", seed, ".tsv", sep = ""))
    
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_workflow" = glmnet_wflow,
                "final_fit" = final_glmnet_fit,
                "predictions" = predictions,
                "performance" = performance,
                "confusion_matrix" = cm,
                "auc" = auc,
                "auc_ci" = ci,
                "roc_curve" = roc, 
                "important_proteins" = important_proteins))
  }


# Functinon to run multiclassification against class
disease_against_class <- 
  function(class, 
           split_train, 
           split_test,
           seed, 
           path_save) {
    
    
    cat(paste0("\nPreparing training and testing data for ", class))
    
    disease_train <- 
      split_train |> 
      left_join(resource_meta |> 
                  select(DAid, Class), by = "DAid") |>
      filter(Class == class,
             Disease %in% include_diseases) |> 
      select(-Class)
    
    disease_test <-   
      split_test |> 
      left_join(resource_meta |> 
                  select(DAid, Class), by = "DAid") |>
      filter(Class == class,
             Disease %in% include_diseases) |> 
      select(-Class)
    
    
    disease_split <- make_splits(disease_train, disease_test)
    
    
    cat(paste0("\nDefining ML specs for ", class))
    
    ### Define general recipe
    ml_recipe <- 
      recipe(Disease ~ ., data = disease_train) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
     # step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric()) 
    
    # Generate resamples
    set.seed(213)
    disease_rs <- vfold_cv(disease_train, v = 10, strata = Disease)
    
    # Define evaluation metrics for all workflows
    eval_metrics <- metric_set(roc_auc)
    
    # Define control grid
    set.seed(213)
    ctrl <- control_grid(verbose = TRUE, 
                         allow_par = TRUE,
                         save_pred = TRUE, 
                         parallel_over = "everything") 
    
    # Tidymodels lasso multiclassification recipe
    glmnet_lasso_specs <-
      multinom_reg() |>
      set_mode("classification") |>
      set_engine("glmnet") |>
      set_args(penalty = tune(),
               mixture = 1)
    
    # Set up lasso workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_lasso_specs) 
    
    # Define hyperparameter tuning grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 20)
    
    cat(paste0("\nFitting glmnet model for ", class))
    
    # Hyperparameter tuning
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = disease_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics)
    
    # Explore results and select the best performing hyperparameter combination
    autoplot(glmnet_res)
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    cat(paste0("\nSelecting best performing model for ", class))
    best_glmnet <- 
      glmnet_res |> 
      select_best(metric = "roc_auc")
    
    # Final fit
    final_glmnet <- 
      glmnet_wflow |> 
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- 
      last_fit(final_glmnet, disease_split)
    
    # Important proteins
    important_proteins <- 
      final_glmnet_fit  |> 
      extract_fit_parsnip() %>%
      tidy() |> 
      filter(term != "(Intercept)") |> 
      arrange(-abs(estimate)) |> 
      filter(abs(estimate) > 0)
    
    write_tsv(important_proteins, paste(path_save, "/important_proteins_", class, "_seed_", seed, ".tsv", sep = ""))
    
      # Explore performance
    final_predictions <- 
      final_glmnet_fit |> 
      collect_predictions() |> 
      mutate(DAid = disease_test$DAid) |> 
      relocate(DAid) |> 
      select(-id)
    
    write_tsv(final_predictions, paste(path_save, "/predictions_", class, "_seed_", seed, ".tsv", sep = "")) 
    
    final_metrics <- 
      final_glmnet_fit |> 
      collect_metrics()
    
    # ROC
    dat <-
      disease_test %>% 
      select(DAid, Disease) %>% 
      mutate(value = 1) %>% 
      spread(Disease,value, fill= 0) 
    
    true_dat <- 
      dat %>% 
      set_names(paste(names(dat), "_true", sep = "")) %>%
      rename(DAid = `DAid_true`)
    
    dat_prob <- 
      final_predictions %>% 
      rename_all(~stringr::str_replace_all(.,".pred_","")) |> 
      select(-".row", -"class", -"Disease", -".config") 
    
    prob_data <- 
      dat_prob %>% 
      set_names(paste(names(dat_prob), "_pred_glmnet", sep = "")) %>% 
      rename(DAid = DAid_pred_glmnet)
    
    final_df <- 
      true_dat %>% 
      left_join(prob_data, by = "DAid") %>% 
      select(-DAid) |> 
      as.data.frame()
    
    roc_res <- multi_roc(final_df, force_diag=T)
    plot_roc_df <- plot_roc_data(roc_res)
    
#    auc_ci <- roc_auc_with_ci(final_df)
    
    roc_dat <- 
      plot_roc_df %>%
      filter(!Group %in% c("Macro","Micro")) %>% 
      mutate(Performance = paste(Group, ": ", round(AUC, 4), sep = "")) %>% 
      arrange(-AUC)
    
 #   roc_dat_ext <- 
  #    roc_dat |> 
   #   left_join(auc_ci |> 
    #              mutate(Group = gsub("glmnet.", "", Var)) |> 
     #             rename(AUC_CI = AUC), by = "Group")
    
    write_tsv(roc_dat, paste(path_save, "/roc_", class, "_seed_", seed, ".tsv", sep = "")) 
    
    roc_plot <- 
      roc_dat %>% 
      arrange(Performance) %>% 
      ggplot(aes(x = 1-Specificity, y=Sensitivity)) +
      geom_path(size=1, show.legend = F, color = pal_class[class]) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   color = "grey",
                   linetype = 'dotdash',
                   show.legend = F) +
      geom_text(aes(label = round(AUC, 2), x = 0.75, y = 0.25), show.legend = F) +
      theme_hpa() +
      scale_y_continuous(breaks = c(0, 1)) +
      scale_x_continuous(breaks = c(0, 1)) +
      facet_wrap(~Group, nrow = 5) 
    
    # Confusion matrix
    cm <- 
      final_glmnet_fit %>% 
      collect_predictions() %>%
      conf_mat(truth = Disease, estimate = .pred_class) 
    
    cm_plot <- 
      cm |>
      autoplot(type = "heatmap") +
      theme_hpa(angled = T)
  
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_fit" = final_glmnet_fit,
                "predictions" = final_predictions,
                "final_metrics" = final_metrics,
                "confusion_matrix" = cm,
                "confusion_matrix_plot" = cm_plot,
                "roc_curve" = roc_dat,
                "roc_plot" = roc_plot,
                "important_proteins" = important_proteins))
    
  }


# Function to run multiclassification against all other diseases
disease_against_all<-  
  function(split_train, 
           split_test,
           seed, 
           path_save) {                    
    
    disease_train <- 
      split_train |> 
      filter(Disease %in% include_diseases)
    
    disease_test <-   
      split_test |>
      filter(Disease %in% include_diseases) 
    
    disease_split <- make_splits(disease_train, disease_test)
    
    
    ### Define general recipe
    ml_recipe <- 
      recipe(Disease ~ ., data = disease_train) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
      #step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric()) 
    
    # Generate resamples
    set.seed(213)
    disease_rs <- vfold_cv(disease_train, v = 10, strata = Disease)
    
    # Define evaluation metrics for all workflows
    eval_metrics <- metric_set(roc_auc)
    
    # Define control grid
    set.seed(213)
    ctrl <- control_grid(verbose = TRUE, 
                         allow_par = TRUE,
                         save_pred = TRUE, 
                         parallel_over = "everything") 
    
    # Tidymodels lasso multiclassification recipe
    glmnet_lasso_specs <-
      multinom_reg() |>
      set_mode("classification") |>
      set_engine("glmnet") |>
      set_args(penalty = tune(),
               mixture = 1)
    
    # Set up lasso workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_lasso_specs) 
    
    # Define hyperparameter tuning grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 20)
    
    # Hyperparameter tuning
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = disease_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics)
    
    # Explore results and select the best performing hyperparameter combination
    autoplot(glmnet_res)
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    best_glmnet <- 
      glmnet_res |> 
      select_best(metric = "roc_auc")
    
    # Final fit
    final_glmnet <- 
      glmnet_wflow |> 
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- 
      last_fit(final_glmnet, disease_split)
    
    # Explore performance
    final_predictions <- 
      final_glmnet_fit |> 
      collect_predictions() |> 
      mutate(DAid = disease_test$DAid) |> 
      relocate(DAid) |> 
      select(-id)
    
    write_tsv(final_predictions, paste(path_save, "predictions_multiclass_seed_", seed, ".tsv", sep = ""))
    
    final_metrics <- 
      final_glmnet_fit |> 
      collect_metrics()
    
    # ROC
    dat <-
      disease_test %>% 
      select(DAid, Disease) %>% 
      mutate(value = 1) %>% 
      spread(Disease,value, fill= 0) 
    
    true_dat <- 
      dat %>% 
      set_names(paste(names(dat), "_true", sep = "")) %>%
      rename(DAid = `DAid_true`)
    
    dat_prob <- 
      final_predictions %>% 
      rename_all(~stringr::str_replace_all(.,".pred_","")) |> 
      select(-".row", -"class", -"Disease", -".config") 
    
    prob_data <- 
      dat_prob %>% 
      set_names(paste(names(dat_prob), "_pred_glmnet", sep = "")) %>% 
      rename(DAid = DAid_pred_glmnet)
    
    final_df <- 
      true_dat %>% 
      left_join(prob_data, by = "DAid") %>% 
      select(-DAid) |> 
      as.data.frame()
    
    roc_res <- multi_roc(final_df, force_diag=T)
    plot_roc_df <- plot_roc_data(roc_res)
    
#    auc_ci <- roc_auc_with_ci(final_df)
    
    roc_dat <- 
      plot_roc_df %>%
      filter(!Group %in% c("Macro","Micro")) %>% 
      mutate(Performance = paste(Group, ": ", round(AUC, 4), sep = "")) %>% 
      arrange(-AUC)
    
 #   roc_dat_ext <- 
  #    roc_dat |> 
   #   left_join(auc_ci |> 
    #              mutate(Group = gsub("glmnet.", "", Var)) |> 
     #             rename(AUC_CI = AUC), by = "Group")
    
    write_tsv(roc_dat, paste(path_save, "roc_multiclass_seed_", seed, ".tsv", sep = "")) 
    #write_tsv(final_predictions, paste(path_save, "/predictions_multiclass_seed_", seed, ".tsv", sep = ""))
    
    roc_plot <- 
      roc_dat %>% 
      left_join(resource_meta |> 
                  distinct(Disease, Class), by = c("Group" = "Disease")) |> 
      # mutate(Group = factor(Group, levels = disease_class_order_roc)) |> 
      ggplot(aes(x = 1-Specificity, y=Sensitivity, color= Class)) +
      geom_path(size=1, show.legend = F) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   color = "grey",
                   linetype = 'dotdash',
                   show.legend = F) +
      geom_text(aes(label = round(AUC, 2), x = 0.75, y = 0.25), show.legend = F) +
      theme_hpa() +
      scale_y_continuous(breaks = c(0, 1)) +
      scale_x_continuous(breaks = c(0, 1)) +
      scale_color_manual(values = pal_class) +
      facet_wrap(~Group, nrow = 5) +
      coord_fixed()
    
    # Confusion matrix
    cm <- 
      final_glmnet_fit %>% 
      collect_predictions() %>%
      conf_mat(truth = Disease, estimate = .pred_class) 
    
    cm_plot <- 
      cm |>
      autoplot(type = "heatmap") +
      theme_hpa(angled = T)
 
    # Important proteins
    important_proteins <- 
      final_glmnet_fit  |> 
      extract_fit_parsnip() %>%
      tidy() |> 
      filter(term != "(Intercept)") |> 
      arrange(-abs(estimate)) |> 
      filter(abs(estimate) > 0)

      write_tsv(important_proteins, paste(path_save, "important_proteins_multiclass_seed_", seed, ".tsv", sep = ""))
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_workflow" = final_glmnet,
                "final_fit" = final_glmnet_fit,
                "predictions" = final_predictions,
                "final_metrics" = final_metrics,
                #"performance" = performance,
                "confusion_matrix" = cm,
                "confusion_matrix_plot" = cm_plot,
                "roc_curve" = roc_dat,
                "roc_plot" = roc_plot,
                "important_proteins" = important_proteins))
    
  }

# Function to run lasso pipeline for prediction of discrete variables
discrete_prediction <-  
  function(variable_predict,
           variable_case,
           split_train, 
           split_test,
           seed, 
           path_save) {
    
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
    
    write_tsv(important_proteins, paste(path_save, "/important_proteins_", variable_predict, "_seed_", seed, ".tsv", sep = ""))
    
    # Extract model predictions
    predictions <- 
      final_glmnet_fit |> 
      collect_predictions(summarize = F) 
    
    # Add performance
    write_tsv(predictions, paste(path_save, "/predictions_", variable_predict, "_seed_", seed, ".tsv", sep = ""))
    
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
    
    write_tsv(combined_roc_auc, paste(path_save, "/roc_", variable_predict, "_seed_", seed, ".tsv", sep = ""))

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


# Function to run lasso pipeline for prediction of continuous variables
continuous_prediction <-  
  function(split_train,
           split_test,
           variable_predict,
           seed, 
           path_save) {
    

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
    
    write_tsv(important_proteins, paste(path_save, "/important_proteins_", variable_predict, "_seed_", seed, ".tsv", sep = ""))
    
    # Extract model predictions
    predictions <- 
      final_glmnet_fit |> 
      collect_predictions(summarize = F) 
    
    # Add performance
    write_tsv(predictions, paste(path_save, "/predictions_", variable_predict, "_seed_", seed, ".tsv", sep = ""))
    
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



