# Load required libraries
library(readr)
library(readxl)
library(dplyr)

# ---- 1. Read input files ----

# Read Somalier output
somalier_data <- read_tsv("somalier.samples.tsv", col_names = TRUE)

# Read Excel manifest (contains CRAM_ID and Sex info)
manifest <- read_excel("Manifest_430_WGS_Lung_Year_2022.xlsx")

# ---- 2. Prepare Somalier sample IDs to match CRAM_ID ----

# Add `.cram` extension to sample_id so it matches manifest$CRAM_ID
somalier_data <- somalier_data %>%
  mutate(CRAM_ID = paste0(sample_id, ".cram"))

# ---- 3. Merge with manifest to get Sex info ----

merged <- somalier_data %>%
  left_join(manifest, by = "CRAM_ID")

# ---- 4. Map Sex to PED codes ----

# PED sex values: 1 = male, 2 = female, -9 = unknown
merged <- merged %>%
  mutate(
    sex_ped = case_when(
      Sex == "Male" ~ 1,
      Sex == "Female" ~ 2,
      TRUE ~ -9
    )
  )

# ---- 5. Create PED dataframe ----

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

# ---- 6. Write PED file with correct header ----

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
