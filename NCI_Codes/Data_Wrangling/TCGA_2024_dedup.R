# Load required packages
library(readxl)
library(dplyr)
library(writexl)

# Define input and output file paths
input_file <- "TCGA_2024_Manifest_05_06_2025.xlsx"
output_excel <- "TCGA_2024_Manifest_05_06_2025_deduplicated.xlsx"
output_txt <- "TCGA_2024_Manifest_05_06_2025_deduplicated.txt"

# Read the Excel file
manifest_data <- read_excel(input_file)

# Remove exact duplicate rows
deduplicated_data <- manifest_data %>%
  distinct()

# Save as Excel
write_xlsx(deduplicated_data, output_excel)

# Save as tab-delimited text
write.table(deduplicated_data, output_txt, sep = "\t", row.names = FALSE, quote = FALSE)

# Optional: message to confirm
cat("Deduplication complete. Files written:\n- ", output_excel, "\n- ", output_txt, "\n")
