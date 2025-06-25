library(dplyr)
library(readxl)
library(stringr)
library(writexl)

# Load the data from the Excel file, skipping the first row (which contains data, not headers)
df <- read_excel("WGS_ID_unanalyzed_batch1_11202024.xlsx", skip = 1)

# Check the column names to confirm they are now correctly set
colnames(df)

# Create a new data frame with the desired columns
df_transformed <- df %>%
  mutate(
    CRAM_ID = basename(cram_file),  # Extract file name from the full file path
    File_ID = cram_file,            # Full file path
    Attribute = `Tissue Attribute`,  # Tissue Attribute column
    Source_Material = `Source Material`,  # Source Material column
    Sherlock_PID = `Sherlock PID`,  # Sherlock PID column
    Site_Study = `Site/Study`      # Site/Study column
  ) %>%
  select(CRAM_ID, File_ID, Attribute, Source_Material, Sherlock_PID, Site_Study)  # Reorder columns

# Save the transformed data as an Excel file
write_xlsx(df_transformed, "transformed_data.xlsx")

