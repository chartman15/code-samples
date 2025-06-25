# Load libraries
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(tools) 

# Load input files
manifest <- read_excel("TCGA_2024_Manifest.xlsx")
aliquots <- read_excel("TCGA_2024_aliquot.xlsx")

# Extract possible IDs from CRAM path and Folder Name
extract_possible_ids <- function(text) {
  uuid <- str_extract_all(text, "[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}")[[1]]
  full_barcode <- str_extract_all(text, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9A-Z]{4}-[0-9A-Z]{3}-[A-Z0-9]{4}-[0-9]{2}")[[1]]
  sample_barcode <- str_extract_all(text, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]")[[1]]
  case_barcode <- str_extract_all(text, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")[[1]]
  all_ids <- unique(c(uuid, full_barcode, sample_barcode, case_barcode))
  toupper(trimws(all_ids))
}

# Combine Folder_Name and CRAM path for ID search
manifest_ids <- manifest %>%
  rowwise() %>%
  mutate(
    Combined_Text = paste(CRAM_File_Path, Folder_Name, sep = " "),
    Matched_IDs = list(extract_possible_ids(Combined_Text)),
    CRAM_File_Name = basename(CRAM_File_Path)
  ) %>%
  unnest(Matched_IDs)

# Normalize aliquots
aliquots <- aliquots %>%
  mutate(row_id = row_number())

aliquots_long <- aliquots %>%
  pivot_longer(cols = -row_id, names_to = "Matching_Column", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(value = toupper(trimws(as.character(value))))

# DEBUG: Print number of IDs that matched
matched_ids <- intersect(
  tolower(manifest_ids$Matched_IDs),
  tolower(aliquots_long$value)
)

unmatched_ids <- setdiff(
  tolower(manifest_ids$Matched_IDs),
  tolower(aliquots_long$value)
)

cat("Matched IDs:", length(matched_ids), "\n")
cat("Unmatched IDs (sample):", paste(head(unmatched_ids, 5), collapse = ", "), "\n")

# Case-insensitive join
matched <- manifest_ids %>%
  mutate(Matched_IDs_lower = tolower(Matched_IDs)) %>%
  left_join(aliquots_long %>%
              mutate(value_lower = tolower(value)),
            by = c("Matched_IDs_lower" = "value_lower")) %>%
  left_join(aliquots, by = "row_id")

# Select final output columns
output_columns <- c("Folder_Name", "CRAM_File_Path", "Matched_IDs", "CRAM_File_Name", "Matching_Column")
if ("cases.case_id" %in% colnames(matched)) {
  output_columns <- c(output_columns, "cases.case_id")
}

final <- matched %>%
  select(all_of(output_columns)) %>%
  rename(Matched_ID = Matched_IDs, Matched_CRAM_Name = CRAM_File_Name) %>%
  distinct()

# Output matched results
write.table(final, "C:/Users/hartmancas/Desktop/TCGA_2024_Manifest_final.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Optional: Save unmatched CRAMs for troubleshooting
unmatched_crams <- manifest_ids %>%
  filter(!(tolower(Matched_IDs) %in% tolower(matched_ids))) %>%
  distinct(CRAM_File_Path, Matched_IDs)

write.table(unmatched_crams, "C:/Users/hartmancas/Desktop/unmatched_cram_ids.txt", sep = "\t", quote = FALSE, row.names = FALSE)
