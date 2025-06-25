# Load required libraries
library(readxl)
library(dplyr)
library(stringr)
library(writexl)  # Optional

# Read the Excel file
input_file <- "dbGAP_Excel_818_Sherlock_Subjects.xlsx"
df <- read_excel(input_file)

# Create the formatted metadata table
subject_df <- df %>%
  mutate(
    SUBJECT_ID = str_remove(Tumor_Barcode, "-T01"),
    SEX = Sex,
    AFFECTION_STATUS = NA,        # Blank
    AGE = NA,                     # Blank
    RACE = NA,                    # Blank
    STAGE = NA,                   # Blank
    GRADE = NA,                   # Blank
    Passive_Smoking = `Smoking Status`,  # Maps directly
    Death = NA,                   # Blank
    Survival_Months = NA,         # Blank
    Histology = NA                # Blank
  ) %>%
  select(SUBJECT_ID, SEX, AFFECTION_STATUS, AGE, RACE, STAGE, GRADE, Passive_Smoking, Death, Survival_Months, Histology)

# Write to TSV
write.table(subject_df, file = "Sherlock_Formatted_Subject_Metadata.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Optional: Write to Excel
# write_xlsx(subject_df, path = "Sherlock_Formatted_Subject_Metadata.xlsx")

# Preview
print(head(subject_df))
