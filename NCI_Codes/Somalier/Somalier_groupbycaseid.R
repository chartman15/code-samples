library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(purrr)

# Load manifest
manifest <- read_excel("New_TCGA_Manifest_Permissions.xlsx")
colnames(manifest) <- trimws(colnames(manifest))  # Remove any whitespace in column names

# Group all Sample_ID_RG_Tag values by case
grouped_all <- manifest %>%
  group_by(`cases.case_id`) %>%
  summarize(
    sample_ids = list(Sample_ID_RG_Tag),
    .groups = "drop"
  )

# Create all-vs-all sample pairs within the same case
expand_all_pairs <- function(ids) {
  if (length(ids) >= 2) {
    # Get all combinations excluding self-pairs
    pairs <- expand.grid(a = ids, b = ids, stringsAsFactors = FALSE) %>%
      filter(a != b)
    
    # Remove duplicate pairs by sorting each row and keeping unique combinations
    unique(t(apply(pairs, 1, sort))) %>%
      as.data.frame() %>%
      setNames(c("sample1", "sample2"))
  } else {
    NULL
  }
}

# Apply pairing logic to all grouped samples
paired_all <- grouped_all %>%
  mutate(pairs = map(sample_ids, expand_all_pairs)) %>%
  select(pairs) %>%
  filter(!map_lgl(pairs, is.null)) %>%
  unnest(pairs)

# Collect all Sample_ID_RG_Tag values used in pairs
paired_ids <- unique(c(paired_all$sample1, paired_all$sample2))

# Identify ungrouped samples
ungrouped <- manifest %>%
  filter(!(Sample_ID_RG_Tag %in% paired_ids)) %>%
  select(Sample_ID_RG_Tag)

# Write output files
write_delim(paired_all, "somalier_same_case_pairs.txt", delim = ",", col_names = FALSE)
write_delim(ungrouped, "ungrouped_samples.txt", delim = "\n", col_names = FALSE)
