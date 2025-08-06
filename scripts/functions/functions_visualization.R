
# Functions for visualization
library(ggrepel)
library(tidytext)
library(embed)
library(ggbeeswarm)
library(patchwork)
library(ggsci)
library(eulerr)
library(ggplotify)
library(pheatmap)
library(ggridges)
library(ggraph)
library(tidygraph)
library(ggupset)

# Generate PCA
do_pca <- function(data,
                   meta = NULL,
                   variable = NULL,
                   wide = T,
                   impute = T,
                   plots = F) {
  if (wide) {
    data_w <-
      data |>
      rename(Sample = DAid)
  } else {
    data_w <-
      data |>
      select(Sample = DAid, Assay, NPX) |>
      pivot_wider(values_from = NPX,
                  names_from = Assay)
  }
  
  if (impute) {
    pca_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_impute_knn(all_predictors()) |>
      step_pca(all_predictors())
    
    pca_prep <- prep(pca_rec)
    
    tidied_pca <- tidy(pca_prep, 3)
    
  } else {
    pca_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_pca(all_predictors())
    
    pca_prep <- prep(pca_rec)
    
    tidied_pca <- tidy(pca_prep, 2)
  }
  loadings_data <-
    tidied_pca |>
    rename(Assay = terms,
           Value = value,
           PC = component)
  
  pca_res <-  juice(pca_prep)
  
  if (plots) {
    # Loadings plot
    loadings_plot <-
      tidied_pca %>%
      filter(component %in% paste0("PC", 1:4)) %>%
      group_by(component) %>%
      top_n(8, abs(value)) %>%
      ungroup() %>%
      mutate(terms = reorder_within(terms, abs(value), component)) %>%
      ggplot(aes(abs(value), terms, fill = value > 0)) +
      geom_col() +
      facet_wrap( ~ component, scales = "free_y") +
      scale_y_reordered() +
      labs(x = "Absolute value of contribution",
           y = NULL, fill = "Positive?") +
      theme_hpa()
    
    # PCA plot
    pca_plot <-
      pca_res %>%
      left_join(meta |>
                  rename(Sample = DAid), by = "Sample") %>%
      ggplot(aes(PC1, PC2)) +
      geom_point(aes(color = !!sym(variable)), alpha = 0.7, size = 2) +
      labs(color = NULL) +
      theme_hpa() +
      labs(color = variable)
    
    return(
      list(
        "pca_res" = pca_res,
        "loadings" = loadings_data,
        "pca_plot" = pca_plot,
        "loadings_plot" = loadings_plot
      )
    )
  } else {
    return(list("pca_res" = pca_res,
                "loadings" = loadings_data))
  }
  
}

# Generate UMAP
do_umap <- function(data,
                    meta = NULL,
                    variable = NULL,
                    wide = T,
                    impute = T,
                    plots = F,
                    seed = 213,
                    n_neighbors = 15) {
  if (wide) {
    data_w <-
      data |>
      rename(Sample = DAid)
  } else {
    data_w <-
      data |>
      select(Sample = DAid, Assay, NPX) |>
      pivot_wider(values_from = NPX,
                  names_from = Assay)
  }
  
  if (impute) {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_impute_knn(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors)
    
    set.seed(seed)
    umap_prep <- prep(umap_rec)
    
  } else {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors)
    
    set.seed(seed)
    umap_prep <- prep(umap_rec)
    
  }
  
  umap_res <-  juice(umap_prep)
  
  if (plots) {
    # Loadings plot
    umap_plot <-
      umap_res |>
      left_join(meta |>
                  rename(Sample = DAid), by = "Sample") |>
      ggplot(aes(UMAP1, UMAP2, color = !!sym(variable))) +
      geom_point(alpha = 0.7, size = 2) +
      theme_hpa()
    
    return(list("umap_res" = umap_res,
                "umap_plot" = umap_plot))
  } else {
    return(umap_res)
  }
  
}

# Generate volcano plot from differential expression results
plot_volcano <-
  function(de_results,
           cutoff = 0,
           labels_balanced = F) {
    if (labels_balanced) {
      labels <-
        de_results |>
        filter(sig == "significant up") |>
        top_n(n = 10, wt = -log10(adj.P.Val))  |>
        bind_rows(de_results |>
                    filter(sig == "significant down") |>
                    top_n(n = 10, wt = -log10(adj.P.Val)))
      
    } else {
      labels <-
        de_results |>
        top_n(n = 10, wt = -log10(adj.P.Val))
    }
    
    
    volcano_plot <-
      de_results |>
      ggplot(aes(
        x = logFC,
        y = -log10(adj.P.Val),
        color = sig,
        label = Assay
      )) +
      geom_point(size = 1,
                 alpha = 0.4,
                 show.legend = F) +
      geom_text_repel(data = labels,
                      size = 2,
                      show.legend = F) +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = 'dashed',
        color = "darkgrey"
      ) +
      geom_vline(
        xintercept = -cutoff,
        linetype = 'dashed',
        color = "darkgrey"
      ) +
      geom_vline(
        xintercept = cutoff,
        linetype = 'dashed',
        color = "darkgrey"
      ) +
      scale_color_manual(values = pal_de) +
      theme_hpa() +
      theme(axis.text = element_text(size = 8),
            legend.position = "top")
    
    return(volcano_plot)
  }

# Plot confusion matrix
plot_cm <- function(confusion_matrix,
                    percentage = F) {
  ann_row <-
    meta_disease |>
    distinct(Class, Disease) |>
    filter(Disease %in% include_diseases) |>
    rename(`True class` = Class) |>
    column_to_rownames("Disease")
  
  ann_col <-
    meta_disease |>
    distinct(Class, Disease) |>
    filter(Disease %in% include_diseases) |>
    rename(`Predicted class` = Class) |>
    column_to_rownames("Disease")
  
  
  cm_dat <-
    confusion_matrix$table |>
    as.data.frame() |>
    group_by(Truth) |>
    mutate(
      Truth = factor(Truth, levels = disease_class_order_ml),
      Prediction = factor(Prediction, levels = disease_class_order_ml)
    )
  
  if (percentage == F) {
    dat <-
      cm_dat |>
      arrange(Truth, Prediction) |>
      pivot_wider(names_from = Truth,
                  values_from = Freq) |>
      column_to_rownames("Prediction") |>
      t() |>
      as.data.frame()
    
    labels <-
      cm_dat |>
      mutate(Freq = ifelse(Freq > 0, Freq, ""),
             Freq = as.character(Freq)) |>
      arrange(Truth, Prediction) |>
      pivot_wider(names_from = Truth,
                  values_from = Freq) |>
      column_to_rownames("Prediction") |>
      t() |>
      as.data.frame()
    
    title <-  "Confusion matrix (n)"
    
    
  } else {
    dat <-
      cm_dat |>
      mutate(Freq = round((Freq / sum(Freq)) * 100, 0)) |>
      arrange(Truth, Prediction) |>
      pivot_wider(names_from = Truth,
                  values_from = Freq) |>
      column_to_rownames("Prediction") |>
      t() |>
      as.data.frame()
    
    labels <-
      cm_dat |>
      group_by(Truth) |>
      mutate(
        Freq = round((Freq / sum(Freq)) * 100, 0),
        Freq = ifelse(Freq > 10, Freq, ""),
        Freq = as.character(Freq)
      ) |>
      arrange(Truth, Prediction) |>
      pivot_wider(names_from = Truth,
                  values_from = Freq) |>
      column_to_rownames("Prediction") |>
      t() |>
      as.data.frame()
    
    title <- "Confusion matrix (%)"
  }
  
  dat |>
    pheatmap(
      cluster_rows = F,
      annotation_row = ann_row,
      annotation_col = ann_col,
      cellwidth = 9,
      cellheight = 9,
      annotation_colors = list("True class" = pal_class, "Predicted class" = pal_class),
      color = c("white", pal_heat),
      display_numbers = labels,
      cluster_cols = F
    ) |>
    as.ggplot() +
    coord_fixed() +
    theme(plot.title = element_text(face = "bold", size = rel(1))) +
    ggtitle(title)
}

# Generate boxplot from Olink data
plot_boxplot <- function(proteins,
                         data,
                         metadata,
                         platform = "HT",
                         title = "") {
  if (platform == "HT") {
    data_filtered <-
      data |>
      filter(Assay %in% proteins) |>
      select(DAid, Assay, NPX) |>
      left_join(metadata |>
                  select(DAid, Cohort, `PI - sample owner`, Disease, Diagnose),
                by = "DAid") |>
      mutate(
        Cohort = paste(Cohort, `PI - sample owner`, sep = "_"),
        Diagnose = ifelse(
          Diagnose == "healthy",
          paste(Diagnose, Cohort, sep = "_"),
          Diagnose
        )
      ) |>
      filter(!is.na(Diagnose))
    
    order <-
      data_filtered |>
      distinct(Cohort, Diagnose) |>
      mutate(Cohort = factor(Cohort, levels = names(pal_phase2))) |>
      arrange(Cohort) |>
      pull(Diagnose)
    
    if (length(proteins) > 1) {
      boxplot <-
        data_filtered |>
        mutate(Diagnose = factor(Diagnose, levels = order)) |>
        filter(!grepl("back-up", Diagnose)) |>
        ggplot(aes(Diagnose,
                   NPX,
                   fill = Cohort,
                   color = Cohort)) +
        geom_quasirandom(alpha = 0.7) +
        geom_boxplot(
          alpha = 0.3,
          outlier.color = NA,
          color = "grey20"
        ) +
        scale_color_manual(values = pal_phase2) +
        scale_fill_manual(values = pal_phase2) +
        facet_wrap( ~ Assay, scales = "free_y") +
        theme_hpa(angled = T) +
        xlab("") +
        ggtitle(title)
      
    } else {
      boxplot <-
        data_filtered |>
        mutate(Diagnose = factor(Diagnose, levels = order)) |>
        filter(!grepl("back-up", Diagnose)) |>
        ggplot(aes(Diagnose,
                   NPX,
                   fill = Cohort,
                   color = Cohort)) +
        geom_quasirandom(alpha = 0.7) +
        geom_boxplot(
          alpha = 0.3,
          outlier.color = NA,
          color = "grey20"
        ) +
        scale_color_manual(values = pal_phase2) +
        scale_fill_manual(values = pal_phase2) +
        theme_hpa(angled = T) +
        ggtitle(title) +
        xlab("")
    }
    
    
  } else if (platform == "1.5K") {
    boxplot <- ggplot()
  } else {
    stop("Platform not recognized")
  }
  
  return(boxplot)
}

# Generate pie overview
plot_donut <-
  function(proteins, type, legend = T) {
    if (type == "secretome") {
      donut_data <-
        tibble(Assay = proteins) |>
        left_join(secretome_hpa |>
                    select(Gene, `Secretome location`),
                  by = c("Assay" = "Gene")) |>
        mutate(
          `Secretome location` = ifelse(
            is.na(`Secretome location`),
            "Not secreted",
            `Secretome location`
          ),
          `Secretome location` = ifelse(
            `Secretome location` == "Intracellular and membrane",
            "Not secreted",
            `Secretome location`
          )
        ) |>
        count(`Secretome location`) |>
        mutate(
          percentage = n / sum(n) * 100,
          label = paste0(`Secretome location`, " (", n, ")")
        )
      
      # Donut plot
      donut_data |>
        ggplot(aes(x = 2, y = n, fill = `Secretome location`)) +
        geom_bar(
          stat = "identity",
          width = 1,
          color = "white",
          show.legend = legend
        ) +
        coord_polar(theta = "y") +
        xlim(1, 2.5) +
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_manual(values = c(pal_secreted, "Not secreted" = "grey80")) +
        geom_text(
          aes(label = ifelse(n > 25, n, "")),
          position = position_stack(vjust = 0.5),
          color = "white",
          size = 4
        ) +
        theme(legend.position = "right")
      
    } else if (type == "secretome-condensed") {
      donut_data <-
        tibble(Assay = proteins) |>
        left_join(secretome_hpa |>
                    select(Gene, `Secretome location`),
                  by = c("Assay" = "Gene")) |>
        mutate(
          `Secretome location` = case_when(
            is.na(`Secretome location`) ~ "Not secreted",
            `Secretome location` == "Intracellular and membrane" ~ "Not secreted",
            `Secretome location` == "Secreted to blood" ~ "Actively secreted to blood",
            T ~ "Secreted to other locations"
          )
        ) |>
        mutate(`Secretome location` = factor(
          `Secretome location`,
          levels = c(
            "Actively secreted to blood",
            "Secreted to other locations",
            "Not secreted"
          )
        )) |>
        count(`Secretome location`) |>
        mutate(
          percentage = n / sum(n) * 100,
          label = paste0(`Secretome location`, " (", n, ")")
        )
      
      # Donut plot
      donut_data |>
        ggplot(aes(x = 2, y = n, fill = `Secretome location`)) +
        geom_bar(
          stat = "identity",
          width = 1,
          color = "white",
          show.legend = legend
        ) +
        coord_polar(theta = "y") +
        xlim(1, 2.5) +
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_manual(values = pal_secretome_condensed) +
        geom_text(
          aes(label = ifelse(n > 25, n, "")),
          position = position_stack(vjust = 0.5),
          color = "white",
          size = 4
        ) +
        theme(legend.position = "right")
    }
    else if (type == "platform") {
      donut_data <-
        tibble(Assay = proteins) |>
        left_join(olink_targets, by = "Assay") |>
        mutate(Platform = ifelse(
          Platform == "Olink Explore 1463",
          "Olink Explore 1536",
          Platform
        )) |>
        count(Platform) |>
        mutate(percentage = n / sum(n) * 100,
               label = paste0(Platform, " (", n, ")"))
      
      # Donut plot
      donut_data |>
        ggplot(aes(x = 2, y = n, fill = Platform)) +
        geom_bar(
          stat = "identity",
          width = 1,
          color = "white",
          show.legend = legend
        ) +
        coord_polar(theta = "y") +
        xlim(1, 2.5) +
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_manual(values = pal_platforms) +
        geom_text(
          aes(label = ifelse(n > 25, n, "")),
          position = position_stack(vjust = 0.5),
          color = "white",
          size = 4
        ) +
        theme(legend.position = "right")
      
    } else if (type == "concentration") {
      donut_data <-
        tibble(Assay = proteins) |>
        left_join(concentrations, by = c("Assay" = "Gene")) |>
        count(unit) |>
        mutate(percentage = n / sum(n) * 100,
               label = paste0(unit, " (", n, ")"))
      
      # Donut plot
      donut_data |>
        ggplot(aes(x = 2, y = n, fill = unit)) +
        geom_bar(stat = "identity",
                 width = 1,
                 color = "white") +
        coord_polar(theta = "y") +
        xlim(1, 2.5) +
        theme_void() +
        theme(legend.position = "none") +
        geom_text(
          aes(label = ifelse(n > 25, n, "")),
          position = position_stack(vjust = 0.5),
          color = "white",
          size = 4
        ) +
        theme(legend.position = "right")
      
      
    } else if (type == "protein class") {
      donut_data <-
        tibble(Assay = proteins) |>
        left_join(secretome_hpa |>
                    select(Gene, `Protein class`),
                  by = c("Assay" = "Gene")) |>
        separate_longer_delim(`Protein class`, delim = ", ") |>
        count(`Protein class`) |>
        mutate(
          percentage = n / sum(n) * 100,
          label = paste0(`Protein class`, " (", n, ")")
        )
      
      # Donut plot
      donut_data |>
        ggplot(aes(x = 2, y = n, fill = `Protein class`)) +
        geom_bar(
          stat = "identity",
          width = 1,
          color = "white",
          show.legend = legend
        ) +
        coord_polar(theta = "y") +
        xlim(1, 2.5) +
        theme_void() +
        theme(legend.position = "none") +
        #  scale_fill_manual(values = c(pal_secreted, "Not secreted" = "grey80")) +
        geom_text(
          aes(label = ifelse(n > 25, n, "")),
          position = position_stack(vjust = 0.5),
          color = "white",
          size = 4
        ) +
        theme(legend.position = "right")
      
      
    } else if (type == "block") {
      donut_data <-
        tibble(Assay = proteins) |>
        left_join(ht_blocks, by = "Assay") |>
        count(Block) |>
        mutate(percentage = n / sum(n) * 100,
               label = paste0(Block, " (", n, ")"))
      
      # Donut plot
      donut_data |>
        mutate(Block = as.factor(Block)) |>
        ggplot(aes(x = 2, y = n, fill = Block)) +
        geom_bar(
          stat = "identity",
          width = 1,
          color = "white",
          show.legend = legend
        ) +
        coord_polar(theta = "y") +
        xlim(1, 2.5) +
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_brewer() +
        geom_text(
          aes(label = ifelse(n > 25, n, "")),
          position = position_stack(vjust = 0.5),
          color = "black",
          size = 4
        ) +
        theme(legend.position = "right")
      
      
    } else {
      warning("Type not recognized")
    }
  }

# Generate boxplot for UKB data
plot_boxplot_ukb <- function(cancer, protein) {
  meta <-
    ukb_dat |>
    left_join(ukb_meta)  |>
    mutate(
      time_to_diagnosis = Age_sample_collect - Cancer_age_diagnosis,
      Cancer_name = ifelse(
        Cancer_name %in% c("Colon_cancer", "Rectum_cancer"),
        "Colorectal_cancer",
        Cancer_name
      ),
      Group = case_when(
        time_to_diagnosis <= -7 ~ "> 7 years before",
        time_to_diagnosis > -7 &
          time_to_diagnosis <= -5 ~ "5-7 years before",
        time_to_diagnosis > -5 &
          time_to_diagnosis <= -3 ~ "3-5 years before",
        time_to_diagnosis > -3 &
          time_to_diagnosis <= -1 ~ "1-3 years before",
        time_to_diagnosis > -1 &
          time_to_diagnosis <= 1 ~ "1 year before/after",
        time_to_diagnosis > 1 &
          time_to_diagnosis <= 3 ~ "1-3 years after",
        time_to_diagnosis > 3 ~ "> 3 years after"
      ),
      Group_2 = factor(Group_2, levels = names(pal_ukb_2))
    )
  
  
  plot_dat <-
    ukb_data |>
    left_join(meta) |>
    filter(Assay == protein,
           Cancer_name == cancer)
  
  plot_dat |>
    ggplot(aes(Group_2, NPX, color = Group_2, fill = Group_2)) +
    geom_boxplot(show.legend = F) +
    #geom_violin(alpha = 0.8) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 95,
      size = 8,
      color = "white",
      show.legend = F
    ) +
    #facet_wrap(Cancer_name~Assay, scales = "free_y", nrow = 1) +
    scale_color_manual(values = pal_ukb_2) +
    scale_fill_manual(values = pal_ukb_2) +
    theme_hpa(axis_x_title = T) +
    ggtitle(paste0(protein, " - ", cancer))
}

# Plot relationship between importance and frequency of selection (ML models)
plot_importance_frequency <-
  function(importances, type, color = "NULL") {
    dat_imp_freq <-
      importances |>
      filter(Importance != 0) |>
      group_by(Variable) |>
      summarise(avg_importance = mean(abs(Importance))) |>
      arrange(-avg_importance)  |>
      left_join(importances |>
                  filter(Importance > 0) |>
                  count(Variable), by = "Variable") |>
      mutate(Type = type)
    
    
    if (color == "coefficient_proportion") {
      dat_imp_freq <-
        dat_imp_freq |>
        left_join(
          importances |>
            mutate(Sign = ifelse(Importance == 0, "NOT IMPORTANT", Sign)) |>
            filter(Sign %in% c("POS", "NEG")) |>
            distinct(Variable, Sign, Seed) |>
            group_by(Variable, Sign) |>
            count(Sign) |>
            ungroup() |>
            group_by(Variable) |>
            mutate(
              total_n = sum(n),
              prop = (n / total_n) * 100,
              max_prop = max(prop)
            ) |>
            distinct(Variable, max_prop),
          by = "Variable"
        )
      
      
      plot <-
        dat_imp_freq |>
        mutate(Proportion_diagreement = 100 - max_prop) |>
        ggplot(aes(avg_importance, n)) +
        geom_point(aes(color = Proportion_diagreement), alpha = 0.5) +
        geom_text_repel(
          data = filter(dat_imp_freq, avg_importance > 0.5),
          aes(label = Variable),
          size = 3
        ) +
        scale_color_viridis_c() +
        ggtitle(type) +
        theme_hpa() +
        xlab("Average importance") +
        ylab("Number of seeds")
      
      
    } else {
      plot <-
        dat_imp_freq |>
        ggplot(aes(avg_importance, n)) +
        geom_point(color = pal_anova[type]) +
        geom_text_repel(
          data = filter(dat_imp_freq, avg_importance > 0.5),
          aes(label = Variable),
          size = 3
        ) +
        theme_hpa() +
        ggtitle(type)
    }
    
    return(plot)
    
  }

# Plot importances for top proteins across seeds (ML models)
plot_top_proteins_seed <- function(protein_importance,
                                   n = 25) {
  top_proteins <-
    protein_importance |>
    filter(Importance > 0) |>
    group_by(Variable) |>
    summarise(avg_importance = mean(abs(Importance)),
              n = n_distinct(Seed)) |>
    arrange(-avg_importance) |>
    head(n)
  
  protein_importance |>
    filter(Variable %in% top_proteins$Variable) |>
    mutate(
      Variable = factor(Variable, levels = rev(top_proteins$Variable)),
      Sign = case_when(Sign == "POS" ~ "Positive",
                       Sign == "NEG" ~ "Negative")
    ) |>
    ggplot(aes(fct_reorder(Variable, abs(Importance)), abs(Importance))) +
    geom_quasirandom(size = 0.5, aes(color = Sign)) +
    geom_boxplot(fill = NA, outlier.color = NA) +
    scale_color_manual(values = c(
      "Positive" = "#FF7176",
      "Negative" = "#92C9DA"
    )) +
    coord_flip() +
    xlab("")  +
    ylab("Protein importance") +
    theme_hpa()
  
}
