library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(purrr)

# Load manifest
manifest <- read_excel("New_TCGA_Manifest_Permissions.xlsx")
colnames(manifest) <- trimws(colnames(manifest))

# Use Sample_ID_RG_Tag instead of Sample_ID
# Group tumor/normal samples by case ID using Sample_ID_RG_Tag
grouped <- manifest %>%
  group_by(`cases.case_id`) %>%
  summarize(
    tumor_ids = list(Sample_ID_RG_Tag[`samples.tissue_type` == "Tumor"]),
    normal_ids = list(Sample_ID_RG_Tag[`samples.tissue_type` == "Normal"]),
    .groups = "drop"
  )

# Generate all tumor-normal pairs per case
expand_pairs <- function(tumors, normals) {
  if (length(tumors) >= 1 && length(normals) >= 1) {
    expand.grid(tumor = tumors, normal = normals, stringsAsFactors = FALSE)
  } else {
    NULL
  }
}

# Apply pairing logic
paired <- grouped %>%
  mutate(pairs = map2(tumor_ids, normal_ids, expand_pairs)) %>%
  select(pairs) %>%
  filter(!map_lgl(pairs, is.null)) %>%
  unnest(pairs)

# Collect all Sample_ID_RG_Tag values used in pairs
paired_ids <- unique(c(paired$tumor, paired$normal))

# Identify ungrouped samples
ungrouped <- manifest %>%
  filter(!(Sample_ID_RG_Tag %in% paired_ids)) %>%
  select(Sample_ID_RG_Tag)

# Write output files
write_delim(paired, "tumor_normal_groups.txt", delim = ",", col_names = FALSE)
write_delim(ungrouped, "ungrouped_samples.txt", delim = "\n", col_names = FALSE)
