# Load required libraries
library(readxl)
library(dplyr)

# ==== 1. Set file paths ====
ped_file <- "somalier_manifest_mapped.ped"
manifest_file <- "TCGA_2024_Manifest_Late_Permissions.xlsx"
output_ped <- "combined_somalier_manifest.ped"

# ==== 2. Read the original PED file ====
ped <- read.delim(ped_file, header = TRUE, stringsAsFactors = FALSE)

# Convert all columns to character (to prevent type mismatch)
ped <- ped %>%
  mutate(across(everything(), as.character))

# ==== 3. Read the manifest Excel file ====
manifest <- read_excel(manifest_file)

# ==== 4. Prepare the manifest data in PED format ====
manifest <- manifest %>%
  mutate(
    sex_code = case_when(
      tolower(demographic.gender) == "male" ~ "1",
      tolower(demographic.gender) == "female" ~ "2",
      TRUE ~ "0"  # unknown or missing sex
    )
  )

# Create new PED-style data from manifest
new_ped <- manifest %>%
  transmute(
    fam_id = as.character(ifelse(is.na(Study), "TCGA", Study)),
    sample_id = as.character(Sample_ID_RG_Tag),
    paternal_id = "0",
    maternal_id = "0",
    sex = sex_code,
    phenotype = "-9"
  )

# ==== 5. Filter to avoid duplicates ====
new_ped_filtered <- new_ped %>%
  filter(!sample_id %in% ped$sample_id)

# ==== 6. Combine original PED with new entries ====
combined_ped <- bind_rows(ped, new_ped_filtered)

# ==== 7. Write to a new PED file ====
write.table(
  combined_ped,
  file = output_ped,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Combined PED file written to:", output_ped, "\n")
