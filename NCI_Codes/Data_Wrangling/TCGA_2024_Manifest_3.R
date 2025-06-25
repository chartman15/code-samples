# Load necessary libraries
library(readxl)
library(dplyr)

# Load data from Excel files
manifest_data <- read_excel("TCGA_2024_Manifest.xlsx")
aliquots_data <- read_excel("TCGA_2024_aliquot.xlsx")

# Extract the core CRAM ID from the file path
manifest_data$CRAM_ID <- gsub(".*\\/([^/]+)(\\.WholeGenome\\.RP-1657\\.cram|_wgs_gdc_realn\\.cram)$", "\\1", manifest_data$CRAM_File_Path)

# First attempt: join by aliquots.submitter_id
merged_data <- manifest_data %>%
  left_join(aliquots_data, by = c("CRAM_ID" = "aliquots.submitter_id"))

# Second attempt: fill in missing case IDs by joining on aliquots.aliquot_id
missing_case_ids <- is.na(merged_data$`cases.case_id`)

if (any(missing_case_ids)) {
  backup_join <- manifest_data %>%
    filter(CRAM_ID %in% aliquots_data$`aliquots.aliquot_id`) %>%
    left_join(aliquots_data, by = c("CRAM_ID" = "aliquots.aliquot_id"))
  
  merged_data[missing_case_ids, "cases.case_id"] <- backup_join$`cases.case_id`[match(
    merged_data$CRAM_ID[missing_case_ids], backup_join$CRAM_ID)]
}

# Filter the final data to desired columns
final_data <- merged_data %>%
  select(Folder_Name, CRAM_File_Path, CRAM_ID, `cases.case_id`)

# Write the final data to a text file (tab-delimited)
write.table(final_data, "Merged_TCGADataset.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Optional: view in console
print(final_data)
