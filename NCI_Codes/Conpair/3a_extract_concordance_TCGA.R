# Load required library
library(stringr)

# Function to extract concordance value from a file
extract.concordance <- function(filename.input) {
  input <- read.delim(filename.input, header = FALSE, stringsAsFactors = FALSE, sep = " ")
  con <- input[1, 2]  # Assumes the concordance value is in the second column of the first row
  return(con)
}

# Define input directory path
dir.input <- "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/output"

# List all sample directories matching "Sample" followed by 3 or more digits
subjects <- list.files(dir.input, pattern = "^Sample\\d{3,}$", full.names = TRUE)

# Construct paths to each expected concordance file
filenames.input <- paste0(subjects, "/", basename(subjects), "_homozygous_concordance.txt")

# Check which files actually exist
idx <- which(file.exists(filenames.input))

# Subset to existing files
subjects <- basename(subjects)[idx]
filenames.input <- filenames.input[idx]

# Name the vector for clarity in output
names(filenames.input) <- subjects

# Extract concordance values
out <- sapply(filenames.input, extract.concordance)

# Create a data frame for output
output <- data.frame(sample = names(out), concordance = out)

# Write results to a tab-delimited text file
write.table(output,
            "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/Conpair_tumour_normal_concordance.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
