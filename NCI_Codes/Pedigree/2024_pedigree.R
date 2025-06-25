# Load required libraries
library(readr)
library(readxl)
library(dplyr)

# ---- 1. Read input files ----

# Read Somalier output (skip header comment lines)
somalier_data <- read_tsv("somalier.samples.tsv", col_names = TRUE)

# Read Excel manifest (contains CRAM_ID and Sex info)
manifest <- read_excel("Manifest_410_WGS_Lung_Year_2024.xlsx")

# Read sample ID map: CRAM_ID <tab> Somalier_ID
id_map <- read_tsv("sample_id_map.tsv", col_names = c("CRAM_ID", "Somalier_ID"))

# ---- 2. Merge all data ----

# Link Somalier sample_id to CRAM_ID
somalier_with_cram <- somalier_data %>%
  left_join(id_map, by = c("sample_id" = "Somalier_ID"))

# Merge with manifest to get Sex info
merged <- somalier_with_cram %>%
  left_join(manifest, by = "CRAM_ID")

# ---- 3. Map Sex to PED codes ----

# PED sex values: 1 = male, 2 = female, -9 = unknown
merged <- merged %>%
  mutate(
    sex_ped = case_when(
      Sex == "Male" ~ 1,
      Sex == "Female" ~ 2,
      TRUE ~ -9
    )
  )

# ---- 4. Create PED dataframe ----

ped <- merged %>%
  select(sample_id, sex_ped) %>%
  mutate(
    family_id = sample_id,
    paternal_id = -9,
    maternal_id = -9,
    phenotype = -9
  ) %>%
  select(family_id, sample_id, paternal_id, maternal_id, sex_ped, phenotype) %>%
  rename(sex = sex_ped)

# ---- 5. Write PED file with correct header ----

# Write header
writeLines("#family_id\tsample_id\tpaternal_id\tmaternal_id\tsex\tphentotype", "somalier_manifest_mapped.ped")

# Append data without column names or row names
write.table(
  ped,
  "somalier_manifest_mapped.ped",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)
