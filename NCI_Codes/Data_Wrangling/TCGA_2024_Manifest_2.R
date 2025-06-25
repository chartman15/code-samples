# Load necessary libraries
library(readxl)
library(dplyr)

# Load data from Excel files using read_excel() from the readxl package
manifest_data <- read_excel("TCGA_2024_Manifest.xlsx")
aliquots_data <- read_excel("TCGA_2024_aliquot.xlsx")

# Extract the core CRAM ID by removing suffixes like .WholeGenome.RP-1657.cram or _wgs_gdc_realn.cram
manifest_data$CRAM_ID <- gsub("(/.*?/)([^/]+)(\\.WholeGenome\\.RP-1657\\.cram|_wgs_gdc_realn\\.cram)$", "\\2", manifest_data$CRAM_File_Path)

# Merge the manifest data with aliquots data based on the CRAM_ID matching the aliquots.submitter_id
merged_data <- manifest_data %>%
  left_join(aliquots_data, by = c("CRAM_ID" = "aliquots.submitter_id"))

# If the join doesn't work directly, try matching with "aliquots.aliquot_id"
merged_data <- merged_data %>%
  left_join(aliquots_data, by = c("CRAM_ID" = "aliquots.aliquot_id"), suffix = c("", "_alt"))

# Filter the final data to only include the desired columns
final_data <- merged_data %>%
  select(Folder_Name, CRAM_File_Path, CRAM_ID, cases.case_id)

# Write the final data to a text file (tab-delimited)
write.table(final_data, "Merged_TCGADataset.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# View the final data (optional)
print(final_data)
