library(dplyr)
library(readxl)
library(readr)
library(tidyr)

# Load manifest
manifest <- read_excel("WGS_ID_unanalyzed_batch1_11202024.xlsx")
colnames(manifest) <- trimws(colnames(manifest))

# Remove .cram extension from CRAM_ID
manifest <- manifest %>%
  mutate(CRAM_ID_clean = sub("\\.cram$", "", CRAM_ID))

# Group by Sex and aggregate CRAM_IDs into list
sex_group <- manifest %>%
  group_by(Sex) %>%
  summarize(CRAM_IDs = list(CRAM_ID_clean), .groups = "drop")

# Create separate files for Male and Female groups
sex_group %>%
  filter(Sex == "Female") %>%
  mutate(line = paste(c("Female", unlist(CRAM_IDs)), collapse = ",")) %>%
  pull(line) %>%
  write_lines("female_groups.txt")

sex_group %>%
  filter(Sex == "Male") %>%
  mutate(line = paste(c("Male", unlist(CRAM_IDs)), collapse = ",")) %>%
  pull(line) %>%
  write_lines("male_groups.txt")
