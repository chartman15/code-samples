library(readxl)
library(dplyr)
library(writexl)

# Load Excel files
meta <- read_excel("Sherlock_Sample_&_Subject_Metadata.xlsx")
wgs <- read_excel("WGS_Analyzed_samples.xlsx")

# Standardize column names to lowercase
colnames(meta) <- tolower(colnames(meta))
colnames(wgs) <- tolower(colnames(wgs))

# Join using sample_id (from meta) == barcode (from wgs)
meta_joined <- meta %>%
  left_join(wgs %>% select(barcode, source_material), by = c("sample_id" = "barcode"))

# Check some matches: print a few rows where source_material was found
cat("Samples with source_material info found:\n")
print(meta_joined %>%
        filter(!is.na(source_material)) %>%
        select(sample_id, body_site, source_material) %>%
        head(10))

# Fill body_site from source_material if missing
meta_final <- meta_joined %>%
  mutate(
    body_site = ifelse(is.na(body_site) | body_site == "", source_material, body_site)
  ) %>%
  select(-source_material)

# Write the updated metadata file
write_xlsx(meta_final, "Sherlock_Sample_Metadata_BodySite_Filled.xlsx")
