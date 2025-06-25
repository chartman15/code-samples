# Load required libraries
library(readxl)
library(dplyr)

# Read Excel files
df_main <- read_excel("WGS_Analyzed_samples.xlsx")
df_race <- read_excel("WGS_Analyzed_samples_race.xlsx")

# Merge on 'subject_id'
df_merged <- df_main %>%
  left_join(df_race %>% select(subject_id, RACE_DERIVED), by = "subject_id")

# Write to a tab-delimited text file with blanks instead of NA
write.table(df_merged, file = "WGS_Analyzed_samples_with_race.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

