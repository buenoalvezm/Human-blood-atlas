
library(RColorBrewer)

# Palettes
pal_detectability <- c("None" = "#DF7176", #FF7176
                       "Some: < 50%" = "#F0BEC1",#E8A29A
                       "Some: > 50%" = "#DFE2E1", #AFB8B4
                       "All" = "#3E6964") #174A45

pal_wellness_individuals <- colorRampPalette(brewer.pal(12, "Paired"))(96)


pal_anova <- c("Age" = "#75C8AE",    
               "Sex" = "#EB7C6A", 
               "BMI" = "#F7B84D",    
               "Disease" = "#EEE2D1")    


pal_anova_wellness <-  c("age" = "#75C8AE",    
                         "sex" = "#EB7C6A", 
                         "visit" = "#E6F0AB",    
                         "subject" = "#ABC8F6",
                         "month_of_sampling" = "#686594")    

pal_sex <- c("F" = "#8B759E", "M" = "#C3E3AF")

pal_mm <- c("Fixed effects (age & sex)" = "#83A49F", 
            "Random effects (subject)" = "#F6C9C5")

pal_heat <- colorRampPalette(c("#FFF3E0",  "#FFCC80", "#B7410E"))(100)


# Palette for disease class
getPalette3 = colorRampPalette(brewer.pal(8, "Set2"))
pal_class<-getPalette3(8)
names(pal_class)<-c("Psychiatric","Cardiovascular","Cancer","Autoimmune","Pediatric","Infection","Metabolic","Healthy") 
class_order <- c("Healthy", "Cardiovascular","Metabolic","Cancer","Psychiatric","Autoimmune","Infection","Pediatric")

rd_bu_continuous <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
pal_cor <- rd_bu_continuous[15:85]

# Levels
female_diseases <- c("Breast cancer", "Breast ductal carcinoma in situ", "Cervical cancer", "Endometrial cancer", "Ovarian cancer")
male_diseases <-  c("Prostate cancer", "Abdominal aortic aneurysm")
pediatric_diseases <- c("Pediatric CNS tumor",
                        "Pediatric bone tumor",
                        "Pediatric retinoblastoma",
                        "Pediatric lymphoma", 
                        "Pediatric neuroblastoma",
                        "Pediatric sarcoma",
                        "Pediatric kidney tumor" , 
                        "Pediatric diffuse astrocytic and oligodendro. tumor",
                        "Pediatric long COVID",
                        "Pediatric systemic inflammatory disease")

pal_ukb <- 
  c("Healthy" = "grey",
    "> 7 years before" = "#4D7654", 
    "5-7 years before" = "#748B5F", 
    "3-5 years before" = "#E3D6A0", 
    "1-3 years before" = "#C78240", 
    "1 year before/after" = "#A42F2E", 
    "1-3 years after" = "#802020", 
    "> 3 years after" = "#510402")

pal_bamse <- 
  c("4" = "#8285BD",
    "8" = "#A993BE",
    "16" = "#F698AA",
    "24" = "#F6C9A5")

pal_phase2 <-  c("#6f1926", "#E8A29A", "#f4895f",  "#f8e16f",  "#95cf92",  "#369acc",  "#9656a2",  "#cbabd1")
names(pal_phase2) <- c("BAMS_Erik Melén", "BDG2_Fredrik Edfors", "CTRL_Fredrik Edfors", "EPIL_Johan Zelano", "FIBR_Camilla Svensson", "PARD_Per Svenningsson", "PREG_Agneta Holmäng","WELL_Göran Bergström" )

pal_de <-
  c("not significant" = "#D3D3D3",
    "significant up" = "#FF7176",
    "significant down" = "#92C9DA")

getPalette3 = colorRampPalette(brewer.pal(8, "Set2"))
pal_class<-getPalette3(8)
names(pal_class)<-c("Psychiatric","Cardiovascular","Cancer","Autoimmune","Pediatric","Infection","Metabolic","Healthy") 

pal_platforms <- c("Olink Explore HT" = "#79B1A8",  # Warm muted orange  
                   "Olink Explore 3072" = "#CB9FC9",
                  "Olink Explore 1563" = "#E8B27F")  # Neutral gray  
# pal_platforms <- 
#   c("Olink Explore HT" = "#D69DCA",
#     "Olink Explore 3072" = "#B2D1B8", 
#     "Olink Explore 1463" = "#A7C7E7")

pal_binary <- 
  c("Yes" = "red",
    "No" = "grey")


pal_secretome_condensed <- c("Actively secreted to blood" = "#B30000", 
                             "Secreted to other locations" = "#F7A391", 
                             "Not secreted" = "grey80")

pal_secreted <- c("Secreted to blood" = "#B30000", 
                  "Secreted in brain" = "#FFDD00", 
                  "Secreted to digestive system" = "#1280C4", 
                  "Secreted in male reproductive system" = "#95D4F5", 
                  "Secreted in female reproductive system" = "#F8BDD7", 
                  "Secreted to extracellular matrix"  = "#7F6A9C", 
                  "Secreted in other tissues" = "#FFD480", 
                  "Secreted - unknown location" = "#A1A8AA", 
                  "Intracellular and membrane" = "#F9A266", 
                  "Immunoglobulin genes" = "#6BA592",
                  "Unknown" = "grey80")

pal_ucan <- c("Lung cancer" = "#ADC74F",
              "Colorectal cancer" = "#B89B74", 
              "Breast cancer" = "#E8A29A",
              "Ovarian cancer" = "#603479",
              "Prostate cancer" = "#E7662B")

pal_ucan_replication <- c("Replicated from Olink Explore 1463" = "#E6998D",
                          "New from Olink Explore HT" = "#5F4A8B",
                          "Other" = "grey")

pal_ucan_replication <- c("Replicated from Olink Explore 1463" = "#F4A261",  # Warm muted orange  
                          "New from Olink Explore HT" = "#2A9D8F",  # Teal/green-blue  
                          "Other" = "#A0A0A0")  # Neutral gray  

pal_platform <- c("Olink Explore 1463" = "#F4A261",  # Warm muted orange  
                  "Olink Explore HT" = "#2A9D8F")  # Neutral gray  


pal_controls<- c("All other diseases" = "#B39BC8",
                 "Class" = "#D9B382",
                 "Healthy"= "#C75D4D")
# Themes
theme_hpa <- 
  function(angled = F, axis_x = T, axis_x_title = F, axis_y = T, facet_title = T) {
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
    
    if(axis_x_title == T) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
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

