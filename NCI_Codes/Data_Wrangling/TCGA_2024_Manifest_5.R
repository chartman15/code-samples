# Load necessary libraries
library(readxl)
library(dplyr)

# Load merged TCGA dataset and sample metadata
merged_data <- read_excel("Merged_TCGADataset.xlsx")
sample_metadata <- read_excel("TGCA_2024_sample.xlsx")

# Perform a left join to bring in the sample metadata
final_merged <- merged_data %>%
  left_join(sample_metadata %>% 
              select(samples.sample_id, 
                     samples.tissue_type, 
                     samples.specimen_type, 
                     samples.sample_type),
            by = "samples.sample_id")

# Export the updated merged dataset to a tab-delimited .txt file
write.table(final_merged, "Merged_TCGADataset_with_Sample_Metadata.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Optional: print a preview
print(head(final_merged))
