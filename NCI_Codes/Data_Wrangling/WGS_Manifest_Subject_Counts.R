library(readxl)
library(dplyr)
library(tidyr)
library(readr)

# Read Excel file
df <- read_excel("WGS_NewSherlock_samples.xlsx")

# Standardize tissue_attribute values
df <- df %>%
  select(Subject, Study, tissue_attribute) %>%
  mutate(tissue_attribute = trimws(tolower(tissue_attribute))) %>%
  mutate(tissue_attribute = case_when(
    tissue_attribute == "tumor" ~ "Tumor",
    tissue_attribute == "normal" ~ "Normal",
    tissue_attribute == "peritumoral" ~ "Peritumoral",
    TRUE ~ NA_character_
  )) %>%
  distinct()

# Determine tissue combinations per subject
subject_summary <- df %>%
  group_by(Study, Subject) %>%
  summarize(
    has_tumor = any(tissue_attribute == "Tumor"),
    has_normal = any(tissue_attribute == "Normal"),
    has_peritumoral = any(tissue_attribute == "Peritumoral"),
    .groups = "drop"
  ) %>%
  mutate(
    category = case_when(
      has_tumor & has_normal & has_peritumoral ~ "Multi-class",
      has_tumor & has_normal ~ "Tumor-Normal Pair",
      has_tumor & has_peritumoral ~ "Tumor-Peritumoral Pair",
      has_normal & has_peritumoral ~ "Normal-Peritumoral Pair",
      has_tumor & !has_normal & !has_peritumoral ~ "Unmatched Tumor",
      has_normal & !has_tumor & !has_peritumoral ~ "Unmatched Normal",
      has_peritumoral & !has_tumor & !has_normal ~ "Unmatched Peritumoral",
      TRUE ~ "Unknown"
    )
  )

# Count subjects by category per study
study_summary <- subject_summary %>%
  count(Study, category) %>%
  pivot_wider(
    names_from = category,
    values_from = n,
    values_fill = 0
  )

# Save summary
write.table(
  study_summary,
  file = "study_summary.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Summary saved to 'study_summary.txt'\n")

# ⬇️ Save unmatched subjects (per study)
unmatched_subjects <- subject_summary %>%
  filter(category %in% c("Unmatched Tumor", "Unmatched Normal", "Unmatched Peritumoral")) %>%
  arrange(Study, category, Subject)

write_tsv(unmatched_subjects, "unmatched_subjects.tsv")
cat("Unmatched subject list saved to 'unmatched_subjects.tsv'\n")

# ⬇️ Identify subjects that appear in more than one study
multi_study_subjects <- df %>%
  distinct(Subject, Study) %>%
  group_by(Subject) %>%
  filter(n() > 1) %>%
  arrange(Subject)

write_tsv(multi_study_subjects, "multi_study_subjects.tsv")
cat("Subjects found in multiple studies saved to 'multi_study_subjects.tsv'\n")
