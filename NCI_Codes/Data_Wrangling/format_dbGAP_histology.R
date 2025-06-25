# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(writexl)  # Optional: only if you want to export to Excel

# Path to your Excel file
input_file <- "dbGAP_Excel_818_Sherlock_Subjects.xlsx"

# Read the Excel file
df <- read_excel(input_file)

# Create a long-format data frame with one row per sample (Tumor or Normal)
formatted_df <- df %>%
  pivot_longer(
    cols = c(Tumor_Barcode, Normal_Barcode),
    names_to = "Sample_Type",
    values_to = "SAMPLE_ID"
  ) %>%
  mutate(
    HISTOLOGICAL_TYPE = "Adenocarcinomas",
    BODY_SITE = "Lung",
    ANALYTE_TYPE = "DNA",
    IS_TUMOR = ifelse(Sample_Type == "Tumor_Barcode", "Tumor", "Normal")
  ) %>%
  select(SAMPLE_ID, HISTOLOGICAL_TYPE, BODY_SITE, ANALYTE_TYPE, IS_TUMOR) %>%
  arrange(SAMPLE_ID)

# Optional: write to TSV
write.table(formatted_df, file = "Sherlock_Sample_Metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Or to Excel if preferred
# write_xlsx(formatted_df, path = "Sherlock_Sample_Metadata.xlsx")

# Print preview
print(head(formatted_df))
