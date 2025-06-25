library(readxl)
library(dplyr)
library(writexl)
library(stringr)

# Read the Excel files
manifest <- read_excel("TCGA_2024_Manifest_Late_Permissions.xlsx")
clinical <- read_excel("TGCA_2024_clinical.xlsx")

# Standardize case ID column names if needed
manifest <- manifest %>%
  rename(cases.case_id = any_of(c("cases_case_id", "case_id", "Case ID", "CaseID")))

clinical <- clinical %>%
  rename(cases.case_id = any_of(c("cases_case_id", "case_id", "Case ID", "CaseID")))

# Clean up case_id formatting (just in case)
manifest <- manifest %>%
  mutate(cases.case_id = str_trim(as.character(cases.case_id)))

clinical <- clinical %>%
  mutate(cases.case_id = str_trim(as.character(cases.case_id)))

# Get gender and race only, ensure one row per case
clinical_subset <- clinical %>%
  select(cases.case_id, demographic.gender, demographic.race) %>%
  distinct(cases.case_id, .keep_all = TRUE)

# Perform a left join (adds columns without duplicating rows)
manifest_updated <- manifest %>%
  left_join(clinical_subset, by = "cases.case_id")

# Write the updated manifest with added gender and race columns
write_xlsx(manifest_updated, "TCGA_2024_Manifest_Late_Permissions_with_Demographics.xlsx")
