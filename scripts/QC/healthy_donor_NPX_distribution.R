#### Title: Batch effect exploration
#### Author: María Bueno Álvez
#### Description: script visualizing NPX distribution across plates for two healthy donors to assess batch effects
#### Last edited : 14/05/2025

source("scripts/functions/functions_utility.R")
source("scripts/functions/functions_visualization.R")
source("scripts/functions/themes_palettes.R")

# Look at control samples - to explore batch effects (raw data)
data_NPX <-
  readr::read_delim(
    "../Human-disease-blood-atlas/data/UD-2950_B1-B4_NPX_2023-03-29/UD-2950_B1-B4_NPX_2023-03-29.csv",
    delim = ";"
  )
manifest <- import_df("data/samples_2025-01-30.xlsx")

# Extract healthy donor data
data_controls <-
  data_NPX |>
  filter(SampleID_Original %in% c("DA_Patient5", "DA_Patient1")) |>
  select(
    SampleID,
    SampleID_Original,
    OlinkID,
    Assay,
    NPX,
    Assay_Warning,
    QC_Warning,
    BatchID,
    PlateID
  )

# Plot plate NPX distribution
order_plates <-
  data_controls |>
  distinct(PlateID, BatchID) |>
  mutate(PlateID = gsub("Run", "", PlateID),
         PlateID = as.numeric(PlateID)) |>
  arrange(BatchID, PlateID)

data_controls |>
  mutate(
    SampleID_Original = case_when(
      SampleID_Original == "DA_Patient1" ~ "Healthy donor 1",
      SampleID_Original == "DA_Patient5" ~ "Healthy donor 2"
    )
  ) |>
  group_by(SampleID_Original, BatchID, Assay) |>
  summarise(NPX = mean(NPX)) |>
  ggplot(aes(x = BatchID, y = NPX)) +
  geom_line(aes(group = interaction(Assay, SampleID_Original)),
            color = "grey",
            alpha = 0.4) +
  geom_quasirandom(aes(fill = BatchID, color = BatchID),
                   dodge.width = 0.7,
                   alpha = 0.3) +
  geom_boxplot(color = "white",
               outlier.color = NA,
               alpha = 0.6) +
  facet_wrap( ~ SampleID_Original, ncol = 1) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_hpa(angled = TRUE)

# Save results
ggsave(
  savepath_results("Manuscript-figures", "rancloud_patient1_patient5.pdf"),
  h = 10 ,
  w = 7
)
