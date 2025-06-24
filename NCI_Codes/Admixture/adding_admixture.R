# Load required libraries
library(readxl)
library(dplyr)
library(writexl)

# Read the Admixture file
admixture_df <- read_excel("New_Sherlock_Admixture.xlsx")

# Read the Manifest file
manifest_df <- read_excel("New_Sherlock_Samples_Manifest.xlsx")

# Merge based only on Sherlock_PID
merged_df <- manifest_df %>%
  left_join(admixture_df %>% 
              select(`Sherlock PID`, `Super Population`) %>% 
              distinct(),  # Ensure no duplicates in the join key
            by = c("Sherlock_PID" = "Sherlock PID"))

# Save the merged result to a new Excel file
write_xlsx(merged_df, "Updated_Sherlock_Samples_Manifest_with_SuperPopulation.xlsx")
