
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
library(plotly)
library(ggridges)
library(ggraph)
library(tidygraph)

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

# Themes
theme_hpa <-
  function(angled = F,
           axis_x = T,
           axis_y = T,
           facet_title = T) {
    t <-
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        plot.title = element_text(
          face = "bold",
          size = rel(1),
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          face = "bold",
          hjust = 0.5,
          size = rel(1),
          vjust = 1
        ),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "cm"),
        legend.title = element_text(size = rel(1)),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = "grey90", fill = "grey90"),
        strip.text = element_text(face = "bold")
      )
    
    if (angled) {
      t <-
        t + theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ))
    }
    
    if (axis_x == F) {
      t <- t +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title.x = element_blank()
        )
    }
    
    if (axis_y == F) {
      t <- t +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank()
        )
    }
    if (facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }
## Generate volcano plot from differential expression results                                                                                                                                  
plot_volcano <- function(de_results, cutoff = 0, labels_balanced = F) {
  
  if (labels_balanced) {
    labels <- 
      de_results |> 
      filter(sig == "significant up") |> 
      top_n(n = 10, wt = -log10(adj.P.Val))  |> 
      bind_rows(de_results |> 
                  filter(sig == "significant down") |> 
                  top_n(n = 10, wt = -log10(adj.P.Val)) 
      )
    
  } else {
    labels <- 
      de_results |> 
      top_n(n = 10, wt = -log10(adj.P.Val)) 
  }
  
  
  volcano_plot <- 
    de_results |> 
    ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = sig, label = Assay)) +
    geom_point(size = 1, alpha = 0.4, show.legend = F) + 
    geom_text_repel(data = labels, size = 2, show.legend = F) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = -cutoff, linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = cutoff, linetype = 'dashed', color = "darkgrey") +
    scale_color_manual(values = pal_de) +
    theme_hpa() +   
    theme(axis.text = element_text(size = 8),
          legend.position = "top") 
  
  return(volcano_plot)
}

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
    mutate(Truth = factor(Truth, levels = disease_class_order_ml),
           Prediction = factor(Prediction, levels = disease_class_order_ml)) 
  
  if(percentage == F) {
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
      mutate(Freq = round((Freq/sum(Freq)) * 100, 0)) |> 
      arrange(Truth, Prediction) |> 
      pivot_wider(names_from = Truth, 
                  values_from = Freq) |>
      column_to_rownames("Prediction") |> 
      t() |> 
      as.data.frame() 
    
    labels <- 
      cm_dat |>
      group_by(Truth) |> 
      mutate(Freq = round((Freq/sum(Freq)) * 100, 0),
             Freq = ifelse(Freq > 10, Freq, ""),
             Freq = as.character(Freq)) |> 
      arrange(Truth, Prediction) |> 
      pivot_wider(names_from = Truth, 
                  values_from = Freq) |>
      column_to_rownames("Prediction") |> 
      t() |> 
      as.data.frame() 
    
    title <- "Confusion matrix (%)"
  }

  dat |> 
    pheatmap(cluster_rows = F, 
             annotation_row = ann_row,
             annotation_col = ann_col,
             cellwidth = 9,
             cellheight = 9, 
             annotation_colors = list("True class" = pal_class, "Predicted class" = pal_class),
             color = c("white", pal_heat),
             display_numbers = labels,
             cluster_cols = F) |> 
    as.ggplot() +
    coord_fixed() +
    theme(plot.title = element_text(face = "bold", size = rel(1))) +
    ggtitle(title) 
}

protein_summary <- function(protein) {
  
  # Overview HPA
  general_overview_table <- 
    hpa_info |> 
    filter(Gene == protein) |> 
    select(Gene, Ensembl, `Gene description`, `Protein class`, `Biological process`, `Molecular function`, `Disease involvement`,
           `Subcellular location`) |> 
    t() |> 
    as.data.frame() |> 
    rownames_to_column("A") |> 
    bind_rows(tibble(A = "Pancreatic cancer database", 
                     V1= ifelse(protein %in% pcd_all$Protein, "Yes", "No"))) |> 
    gt() |>
    tab_header(
      title = md("Protein summary based on HPA data")
    ) |> 
    tab_options(column_labels.hidden = TRUE) |> 
    tab_style(
      style = list(
        cell_fill(color = "grey90"),
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        columns = A)
    )
  
  gtsave(general_overview_table, paste0(protein, "_overview.png"), path = "results/gt_temp/")
  general_overview <-  ggdraw() + draw_image(paste0("results/gt_temp/", protein, "_overview.png"), scale = 0.8)
  
  # Literature
  ids <- get_pubmed_ids(paste0(protein, '[All Fields] AND "pancreatic cancer"[All Fields]'))
  abstracts_xml <- fetch_pubmed_data(pubmed_id_list = ids)
 
  
  if(!is.null(abstracts_xml)) {
    abstracts_list <- articles_to_list(abstracts_xml)
    abstract_table <- 
      abstracts_list |> 
      map_df(function(abstract) {
        article_to_df(pubmedArticle = abstract, autofill = FALSE) |> 
          head(1) |> 
          select(pmid, year, journal, firstname, lastname, title) |> 
          mutate(First_author = paste0(lastname, ", ", firstname)) |> 
          relocate(First_author, .before = year) |> 
          select(-firstname, -lastname)
      })
    
    table <-
      abstract_table |>
      head(10) |> 
      gt() |>
      tab_header(
        title = md("**PubMed search**")) |> 
      cols_label(
        pmid = md("**PubMed**"),
        First_author = md("**Author**"),
        year = md("**Year**"),
        journal = md("**Journal**"),
        title = md("**Title**")
      )
    
    gtsave(table, paste0(protein, "_table.png"), path = "results/")
    literature <-  ggdraw() + draw_image(paste0("results/", protein, "_table.png"), scale = 0.8)
    
  } else {
    
    literature <- ggplot() + theme_void()
  }
  
  # Tissue RNA
  healthy_rna <- 
    hpa_rna |> 
    filter(`Gene name` == protein) |> 
    mutate(Tissue = factor(Tissue, levels = names(pal_tissue[order(pal_tissue)]))) |> 
    ggplot(aes(Tissue, nTPM, fill = Tissue)) +
    geom_col(show.legend = F) +
    scale_fill_manual(values = pal_tissue[order(pal_tissue)]) +
    theme_hpa(angled = T) +
    ggtitle("Healthy data - HPA")
  
  ensgid <- 
    hpa_rna |>
    distinct(Gene, `Gene name`) |> 
    filter(`Gene name` == protein) |> 
    pull(Gene)
  pal_hpa <- read_tsv("data/hpa/hpa_colors.txt", col_names = F)
  
  pal_cancer <- 
    pal_hpa |> 
    filter(grepl("cancer", X1)) |> 
    select(-X2) |> 
    deframe() 
  
  pal_cancer_complete <- 
    c(pal_cancer,
      "Endometrial cancer" = "#F8BDD7",
      "Glioma" = "#FFDD00",
      "Head and neck cancer" = "#1280C4",
      "Melanoma" = "#FCCAB3",
      "Renal cancer" = "#F9A266",
      "Stomach cancer" = "#1280C4",
      "Urothelial cancer" = "#F9A266") |> 
    sort()
  
  tumor_rna  <- 
    hpa_tcga |> 
    filter(Gene == ensgid) |> 
    mutate(Cancer = case_when(Cancer == "BRCA" ~ "Breast cancer",  
                              Cancer == "CESC" ~ "Cervical cancer",  
                              Cancer == "COAD" ~ "Colorectal cancer",  
                              Cancer == "READ" ~ "Colorectal cancer",  
                              Cancer == "UCEC" ~ "Endometrial cancer",  
                              Cancer == "GBM" ~ "Glioma",  
                              Cancer == "HNSC" ~ "Head and neck cancer",  
                              Cancer == "LIHC" ~ "Liver cancer",  
                              Cancer == "LUAD" ~ "Lung cancer",  
                              Cancer == "LUSC" ~ "Lung cancer",  
                              Cancer == "SKCM" ~ "Melanoma",  
                              Cancer == "OV" ~ "Ovarian cancer",  
                              Cancer == "PAAD" ~ "Pancreatic cancer",  
                              Cancer == "PRAD" ~ "Prostate cancer",  
                              Cancer == "KICH" ~ "Renal cancer",  
                              Cancer == "KIRC" ~ "Renal cancer",  
                              Cancer == "KIRP" ~ "Renal cancer",  
                              Cancer == "STAD" ~ "Stomach cancer",  
                              Cancer == "TGCT" ~ "Testis cancer",  
                              Cancer == "THCA" ~ "Thyroid cancer",  
                              Cancer == "BLCA" ~ "Urothelial cancer"), 
           Cancer = factor(Cancer, levels = names(pal_cancer_complete))) |> 
    ggplot(aes(Cancer, FPKM, color = Cancer, fill = Cancer)) +
    geom_quasirandom(alpha = 0.6, show.legend = F, size = 0.5) +
    geom_boxplot(alpha = 0.4, color = "black", show.legend = F, outlier.color = NA, linewidth = 0.2) +
    scale_fill_manual(values = pal_cancer_complete) +
    scale_color_manual(values = pal_cancer_complete) +
    theme_hpa(angled = T) +
    ggtitle("Cancer data - TCGA")
  
  RNA_data <- 
    healthy_rna + tumor_rna +
    plot_layout(widths = c(2,1)) +
    plot_annotation(title = "Tissue RNA data") & theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                       axis.line = element_line(linewidth = 0.3),
                                                       axis.ticks = element_line(linewidth = 0.3), 
                                                       text = element_text(size = 7))
  
  # Tissue proteomics
  if(protein %in% colnames(cptac_data_paired)) {
    cptac <- 
      cptac_data_paired |> 
      select(Sample, Type, Protein_expression = protein) |>
      mutate(Disease = recode(Type,
                              Tumor = "Pancreatic cancer",
                              Normal = "Healthy"),
             Protein_expression = as.numeric(Protein_expression),
             Dataset = "CPTAC") |> 
      select(-Type) |>
      mutate(Disease = factor(Disease, levels = c(names(pal_group), names(pal_ukb)) |> unique())) |> 
      ggplot(aes(Disease,Protein_expression, color = Disease, fill = Disease)) +
      geom_quasirandom(alpha = 0.8, show.legend = F, size = 0.5) +
      geom_boxplot(alpha = 0.2, color = "black", outlier.colour = NA, show.legend = F, linewidth = 0.2) +
      theme_hpa(angled = T) +
      scale_color_manual(values = c(pal_group, pal_ukb)) +
      scale_fill_manual(values = c(pal_group, pal_ukb)) +
      xlab("") +
      ylab("Expression") +
      ggtitle("Healthy & tumor - CPTAC")
  } else {
    cptac <-  ggplot() + theme_void()
  }
 
  
  
  hpa_prot <- 
    hpa_protein |>
    mutate(Type = paste(Tissue, `Cell type`, sep = "_"),
           Level = factor(Level, levels = c("Not detected", "Low", "Medium", "High"))) |>
    filter(`Gene name` == protein) |>
    mutate(Tissue = case_when(Tissue %in% c("endometrium 1", "endometrium 2") ~ "endometrium",
                              Tissue %in% c("skin 1", "skin 2") ~ "skin",
                              Tissue %in% c("soft tissue 1", "soft tissue 2") ~ "soft tissue", 
                              Tissue %in% c("stomach 1", "stomach 2") ~ "soft tissue",
                              T ~ Tissue),
           n = case_when(Level == "Not detected" ~ 0,
                         Level == "Low" ~ 1,
                         Level == "Medium" ~ 2,
                         Level == "High" ~ 3),
           Tissue = factor(Tissue, levels = names(pal_tissue[order(pal_tissue)]))) |> 
    group_by(Tissue) |> 
    top_n(1, n) |>
    select(Tissue, Level, n) |> 
    distinct() |> 
    filter(!is.na(Tissue)) |> 
    ggplot(aes(Tissue, n, fill = Tissue)) +
    geom_col(show.legend = F) +
    #  scale_fill_manual(values = c("white", "grey80", "grey50", "grey30")) +
    scale_fill_manual(values = pal_tissue[order(pal_tissue)]) +
    theme_hpa(angled = T) +
    scale_y_continuous(breaks = c(0,1,2,3),
                       labels = c("Not detected", "Low", "Medium", "High")) +
    ggtitle("Healthy data - HPA") +
    xlab("") +
    ylab("")
  
  pathology <- 
    hpa_pathology |> 
    filter(`Gene name` == protein) |> 
    select(c(1:7)) |> 
    pivot_longer(cols = c(4:7), 
                 names_to = "Staining", 
                 values_to = "n") |> 
    group_by(Cancer) |> 
    mutate(total = sum(n),
           perc = n/total * 100,
           Staining = factor(Staining, levels = c("Not detected", "Low", "Medium", "High"))) |> 
    ggplot(aes(Cancer, perc, fill = Staining)) +
    geom_col() +
    scale_fill_manual(values = c("white", "grey80", "grey50", "grey30")) +
    theme_hpa(angled = T) +
    xlab("") +
    ylab("Patients (%)") +
    theme(legend.position = "top") +
    ggtitle("Cancer data - HPA")
  
  
  tissue_proteomics <- 
    cptac + hpa_prot + pathology +
    plot_layout(widths = c(1,4,3)) +
    plot_annotation(title = "Tissue proteomics data") & theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                              axis.line = element_line(linewidth = 0.3),
                                                              axis.ticks = element_line(linewidth = 0.3), 
                                                              text = element_text(size = 7)) 
  
  # Plasma proteomics
  plasma_hpa <- 
    cancer_data_long |>
    filter(Assay == protein) |>
    rename(GROUP = Disease) |> 
    mutate(Disease = case_when(GROUP == "PAN" ~ "Pancreatic cancer",
                               GROUP == "CAD healthy" ~ "Healthy",
                               T ~ "Other cancers")) |> 
    left_join(cancers_mapping |> 
                mutate(Cancer = ifelse(Cancer_code %in% c("AML", "CLL"), Cancer_code, Cancer),
                       Cancer = ifelse(Cancer_code %in% c("LYMPH"), "DLBCL", Cancer)), 
              by = c("GROUP" = "Cancer_code")) |> 
    mutate(Cancer = ifelse(Disease == "Healthy", "Healthy", Cancer),
           Cancer = factor(Cancer, levels = c(long_levels, "Healthy"))
    ) |> 
    select(Assay, NPX, Disease, Cancer) |> 
    ggplot(aes(Cancer, NPX, color = Disease, fill = Disease)) +
    geom_quasirandom(alpha = 0.8, show.legend = F, size = 0.5) +
    geom_boxplot(alpha = 0.2, color = "black", outlier.colour = NA, show.legend = F, linewidth = 0.2) +
    theme_hpa(angled = T) +
    scale_color_manual(values = pal_group) +
    scale_fill_manual(values = pal_group) +
    xlab("") +
    ylab("NPX") +
    ggtitle("Cancers & Healthy - HPA")
  
  
  plasma_ukb <- 
    ukb_data |> 
    filter(Assay == protein) |> 
    left_join(ukb_meta_ext, by = "id") |> 
    select(Sample = id, Protein_expression = NPX, Disease = Group) |> 
    mutate(Dataset = "UKB",
           Sample = as.character(Sample)) |> 
    mutate(Disease = factor(Disease, levels = c(names(pal_group), names(pal_ukb)) |> unique())) |> 
    ggplot(aes(Disease,Protein_expression, color = Disease, fill = Disease)) +
    geom_quasirandom(alpha = 0.8, show.legend = F, size = 0.5) +
    geom_boxplot(alpha = 0.2, color = "black", outlier.colour = NA, show.legend = F, linewidth = 0.2) +
    theme_hpa(angled = T) +
    scale_color_manual(values = c(pal_group, pal_ukb)) +
    scale_fill_manual(values = c(pal_group, pal_ukb)) +
    xlab("") +
    ylab("NPX") +
    ggtitle("UKB-PPP cohort")
  
  plasma_ukb_continuous <- 
    ukb_data |> 
    filter(Assay == protein,
           Disease != "Healthy") |>
    left_join(ukb_meta_ext, by = "id") |> 
    select(Sample = id, Protein_expression = NPX, Disease = Group, Difference) |> 
    ggplot(aes(Difference,Protein_expression, color = Disease)) +
    geom_point(show.legend = F, size = 0.5) +
    geom_smooth(inherit.aes = F, aes(Difference, Protein_expression), color = "black", linewidth = 0.5) +
    theme_hpa(angled = T) +
    scale_color_manual(values = c(pal_group, pal_ukb)) +
    xlab("Years") +
    ylab("NPX") +
    scale_x_reverse() 
  
  plasma_proteomics <- 
    plasma_hpa + plasma_ukb + plasma_ukb_continuous +
    plot_layout(widths = c(3,2,1)) +
    plot_annotation(title = "Plasma proteomics data") & theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                                              axis.line = element_line(linewidth = 0.3),
                                                              axis.ticks = element_line(linewidth = 0.3), 
                                                              text = element_text(size = 7))
  
  # Titles
  title0 <- 
    ggplot() +
    geom_text(aes(0.5, 0.5, label = protein), fontface = "bold", color = "grey20") +
    theme_void() +
    theme(plot.background = element_rect(fill = "grey80"))
  
  title1 <- 
    ggplot() +
    geom_text(aes(0.5, 0.5, label = "Overview"), fontface = "bold", color = "white") +
    theme_void() +
    theme(plot.background = element_rect(fill = "grey30"))
  
  title2 <- 
    ggplot() +
    geom_text(aes(0.5, 0.5, label = "Tissue RNA data"), fontface = "bold", color = "white") +
    theme_void() +
    theme(plot.background = element_rect(fill = "grey30"))
  
  title3 <- 
    ggplot() +
    geom_text(aes(0.5, 0.5, label = "Tissue protein data"), fontface = "bold", color = "white") +
    theme_void() +
    theme(plot.background = element_rect(fill = "grey30"))
  
  title4 <- 
    ggplot() +
    geom_text(aes(0.5, 0.5, label = "Plasma protein data"), fontface = "bold", color = "white") +
    theme_void() +
    theme(plot.background = element_rect(fill = "grey30"))
  
  # Combine
  first <- 
    general_overview | literature
  
  title0/ title1 / first / title2 / RNA_data / title3 / tissue_proteomics / title4  / plasma_proteomics  +
    plot_layout(heights = c(0.3,0.3,2.5,0.3,1.5,0.3,1.5,0.3,1.5))
  
}


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
        ggplot(aes(
          Diagnose,
          NPX,
          fill = Cohort,
          color = Cohort
        )) +
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
        ggplot(aes(
          Diagnose,
          NPX,
          fill = Cohort,
          color = Cohort
        )) +
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

plot_donut <- 
  function(proteins, type, legend = T) {
    
    
    if(type == "secretome") {
      
      donut_data <- 
        tibble(Assay = proteins) |> 
        left_join(secretome_hpa |> 
                    select(Gene, `Secretome location`), by = c("Assay" = "Gene")) |> 
        mutate(`Secretome location` = ifelse(is.na(`Secretome location`), "Not secreted", `Secretome location`),
               `Secretome location` = ifelse(`Secretome location` == "Intracellular and membrane", "Not secreted", `Secretome location`)) |> 
        count(`Secretome location`) |> 
        mutate(percentage = n / sum(n) * 100, 
               label = paste0(`Secretome location`, " (", n, ")"))  
      
      # Donut plot
      donut_data |> 
        ggplot(aes(x = 2, y = n, fill = `Secretome location`)) + 
        geom_bar(stat = "identity", width = 1, color = "white", show.legend = legend) +  
        coord_polar(theta = "y") +  
        xlim(1, 2.5) +  
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_manual(values = c(pal_secreted, "Not secreted" = "grey80")) +
        geom_text(aes(label = ifelse(n > 25, n, "")),  
                  position = position_stack(vjust = 0.5),  
                  color = "white", size = 4) +
        theme(legend.position = "right")
      
    } else if (type == "secretome-condensed") {
    
      donut_data <- 
        tibble(Assay = proteins) |> 
        left_join(secretome_hpa |> 
                    select(Gene, `Secretome location`), by = c("Assay" = "Gene")) |> 
        mutate(`Secretome location` = case_when(is.na(`Secretome location`) ~ "Not secreted", 
                                                `Secretome location` == "Intracellular and membrane" ~ "Not secreted", 
                                                `Secretome location`== "Secreted to blood" ~ "Actively secreted to blood",
                                                T ~ "Secreted to other locations")) |> 
        mutate(`Secretome location` = factor(`Secretome location`, levels = c("Actively secreted to blood", "Secreted to other locations", "Not secreted"))) |> 
        count(`Secretome location`) |> 
        mutate(percentage = n / sum(n) * 100, 
               label = paste0(`Secretome location`, " (", n, ")"))  
      
      # Donut plot
      donut_data |> 
        ggplot(aes(x = 2, y = n, fill = `Secretome location`)) + 
        geom_bar(stat = "identity", width = 1, color = "white", show.legend = legend) +  
        coord_polar(theta = "y") +  
        xlim(1, 2.5) +  
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_manual(values = pal_secretome_condensed) +
        geom_text(aes(label = ifelse(n > 25, n, "")),  
                  position = position_stack(vjust = 0.5),  
                  color = "white", size = 4) +
        theme(legend.position = "right")
    }
    else if (type == "platform") {
      
      donut_data <- 
        tibble(Assay = proteins) |> 
        left_join(olink_targets, by = "Assay") |> 
        count(Platform) |> 
        mutate(percentage = n / sum(n) * 100, 
               label = paste0(Platform, " (", n, ")"))  
      
      # Donut plot
      donut_data |> 
        ggplot(aes(x = 2, y = n, fill = Platform)) + 
        geom_bar(stat = "identity", width = 1, color = "white", show.legend = legend) +  
        coord_polar(theta = "y") +  
        xlim(1, 2.5) +  
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_manual(values = pal_platforms) +
        geom_text(aes(label = ifelse(n > 25, n, "")),  
                  position = position_stack(vjust = 0.5),  
                  color = "white", size = 4) +
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
        geom_bar(stat = "identity", width = 1, color = "white") +  
        coord_polar(theta = "y") +  
        xlim(1, 2.5) +  
        theme_void() +
        theme(legend.position = "none") +
       # scale_fill_brewer() +
#        scale_fill_manual(values = pal_platforms) +
        geom_text(aes(label = ifelse(n > 25, n, "")),  
                  position = position_stack(vjust = 0.5),  
                  color = "white", size = 4) +
        theme(legend.position = "right")
      
      
    }else if (type == "protein class") {
      
     donut_data <- 
        tibble(Assay = proteins) |> 
        left_join(secretome_hpa |> 
                    select(Gene, `Protein class`), by = c("Assay" = "Gene")) |> 
        separate_longer_delim(`Protein class`, delim = ", ") |> 
        count(`Protein class`) |> 
        mutate(percentage = n / sum(n) * 100, 
               label = paste0(`Protein class`, " (", n, ")"))  
      
      # Donut plot
      donut_data |> 
        ggplot(aes(x = 2, y = n, fill = `Protein class`)) + 
        geom_bar(stat = "identity", width = 1, color = "white", show.legend = legend) +  
        coord_polar(theta = "y") +  
        xlim(1, 2.5) +  
        theme_void() +
        theme(legend.position = "none") +
      #  scale_fill_manual(values = c(pal_secreted, "Not secreted" = "grey80")) +
        geom_text(aes(label = ifelse(n > 25, n, "")),  
                  position = position_stack(vjust = 0.5),  
                  color = "white", size = 4) +
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
        geom_bar(stat = "identity", width = 1, color = "white", show.legend = legend) +  
        coord_polar(theta = "y") +  
        xlim(1, 2.5) +  
        theme_void() +
        theme(legend.position = "none") +
        scale_fill_brewer() +
        geom_text(aes(label = ifelse(n > 25, n, "")),  
                  position = position_stack(vjust = 0.5),  
                  color = "black", size = 4) +
        theme(legend.position = "right")
      
      
    } else {
      warning("Type not recognized")
    }
   
    
    
  }

boxplot_ukb <- function(cancer, protein) {
  meta <- 
    ukb_dat |> 
    left_join(ukb_meta)  |> 
    mutate(time_to_diagnosis = Age_sample_collect - Cancer_age_diagnosis,
           Cancer_name = ifelse(Cancer_name %in% c("Colon_cancer", "Rectum_cancer"), "Colorectal_cancer", Cancer_name),
           Group = case_when(time_to_diagnosis <= -7 ~ "> 7 years before",
                             time_to_diagnosis > -7 & time_to_diagnosis <= -5 ~ "5-7 years before",
                             time_to_diagnosis > -5 & time_to_diagnosis <= -3 ~ "3-5 years before",
                             time_to_diagnosis > -3 & time_to_diagnosis <= -1 ~ "1-3 years before",
                             time_to_diagnosis > -1 & time_to_diagnosis <= 1 ~ "1 year before/after",
                             time_to_diagnosis > 1 & time_to_diagnosis <= 3 ~ "1-3 years after",
                             time_to_diagnosis > 3 ~ "> 3 years after"),
           Group = factor(Group, levels = names(pal_ukb_2[-1]))) 
  
  
  plot_dat <- 
    ukb_data |>
    left_join(meta) |> 
    filter(Assay == protein,
           Cancer_name == cancer) 
  
  plot_dat |> 
    ggplot(aes(Group, NPX, color = Group, fill = Group)) +
    geom_boxplot(show.legend = F) +
    #geom_violin(alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 95, size = 8, color = "white", show.legend = F) +
    #facet_wrap(Cancer_name~Assay, scales = "free_y", nrow = 1) +
    scale_color_manual(values = pal_ukb_2) +
    scale_fill_manual(values = pal_ukb_2) +
    theme_hpa(axis_x_title = T) +
    ggtitle(paste0(protein, " - ", cancer))
  # theme(axis.ticks = element_blank(),
  #      axis.text = element_blank()) #+
  #  ylab("Protein expression") +
  # xlab("Time to diagnosis") 
  
}


plot_importance_frequency <- function(importances, type) {
  
  dat_imp_freq <- 
    importances |> 
    group_by(Variable) |> 
    summarise(avg_importance = mean(abs(Importance))) |> 
    arrange(-avg_importance)  |> 
    left_join(importances |>
                filter(Importance > 0) |> 
                count(Variable), by = "Variable") |> 
    mutate(Type = type) 
  
  
  dat_imp_freq |> 
    ggplot(aes(avg_importance, n)) +
    geom_point(color = pal_anova[type]) +
    geom_text_repel(data = filter(dat_imp_freq, avg_importance > 0.5),
                    aes(label = Variable), size = 3) +
    theme_hpa() +
    ggtitle(type)
}

plot_top_proteins_seed <- function(protein_importance, 
                                   plot_color = "grey20",
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
    mutate(Variable = factor(Variable, levels = rev(top_proteins$Variable))) |> 
    ggplot(aes(fct_reorder(Variable, abs(Importance)), abs(Importance))) +
    geom_quasirandom(size = 0.5, color = plot_color) +
    geom_boxplot(alpha = 0.5, outlier.color = NA, fill = plot_color) +
    coord_flip() +
    xlab("")  +
    ylab("Protein importance") +
    theme_hpa() 
  
}

