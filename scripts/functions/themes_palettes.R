
library(RColorBrewer)

# Palettes
pal_detectability <- c("None" = "#FF7176",
                       "Some: < 50%" = "#E8A29A",
                       "Some: > 50%" = "#AFB8B4",
                       "All" = "#174A45")

pal_anova_wellness <-  c("age" = "#75C8AE",    
                         "sex" = "#EB7C6A", 
                         "visit" = "#E6F0AB",    
                         "subject" = "#ABC8F6",
                         "month_of_sampling" = "#686594")    

pal_sex <- c("F" = "#8B759E", "M" = "#A3D0B4")


pal_ukb_2 <- 
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

pal_platforms <- 
  c("HT" = "#D69DCA",
    "3K" = "#D6EDDA", 
    "1.5K" = "#A7C7E7")

pal_binary <- 
  c("Yes" = "red",
    "No" = "grey")

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

pal_ucan <- c("LUNG" = "#ADC74F",
              "CRC" = "#B89B74", 
              "BRC" = "#E8A29A",
              "OVC" = "#603479",
              "PRC" = "#E7662B")

# Themes
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

