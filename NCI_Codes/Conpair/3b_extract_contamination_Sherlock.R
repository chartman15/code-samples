library(stringr)

# Function to extract contamination values from the file
extract.contamination <- function(filename.input) {
  # Read the entire file
  input <- readLines(filename.input)
  
  # Use regular expressions to extract the contamination values
  normal_contamination <- as.numeric(sub("Normal sample contamination level: ([0-9.]+)%", "\\1", input[grep("Normal sample contamination level", input)]))
  tumor_contamination <- as.numeric(sub("Tumor sample contamination level: ([0-9.]+)%", "\\1", input[grep("Tumor sample contamination level", input)]))
  
  # Return both values as a named vector (Normal, Tumor)
  return(c(Normal = normal_contamination, Tumor = tumor_contamination))
}

# Set the directory input path
dir.input <- "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/output/Sherlock"

# Get the list of sample directories (Sample001, Sample002, etc.)
subjects <- list.files(dir.input, pattern = "^Sample\\d{3}$", full.names = TRUE)

# Create the corresponding file paths for contamination files
filenames.input <- paste0(subjects, "/", basename(subjects), "_cross_contamination.txt")

# Check which contamination files exist
idx <- which(file.exists(filenames.input))

# Filter the subjects and filenames based on existing files
subjects <- basename(subjects)[idx]
filenames.input <- filenames.input[idx]

# Create a named vector for the contamination values
names(filenames.input) <- subjects

# Extract contamination values for each subject using sapply
contamination_values <- t(sapply(filenames.input, extract.contamination))

# Create a data frame with sample names and corresponding contamination values
contamination_data <- data.frame(sample = rownames(contamination_values), contamination_values)

# Calculate summary statistics for Normal and Tumor contamination values
mean_normal <- mean(contamination_values[, "Normal"], na.rm = TRUE)
median_normal <- median(contamination_values[, "Normal"], na.rm = TRUE)
sd_normal <- sd(contamination_values[, "Normal"], na.rm = TRUE)
mean_tumor <- mean(contamination_values[, "Tumor"], na.rm = TRUE)
median_tumor <- median(contamination_values[, "Tumor"], na.rm = TRUE)
sd_tumor <- sd(contamination_values[, "Tumor"], na.rm = TRUE)

# Create a summary table for contamination values
summary_stats <- data.frame(
  Normal_Mean = mean_normal,
  Normal_Median = median_normal,
  Normal_SD = sd_normal,
  Tumor_Mean = mean_tumor,
  Tumor_Median = median_tumor,
  Tumor_SD = sd_tumor
)

# Print summary statistics to console
print(summary_stats)

# Write the contamination data to a file
write.table(contamination_data, "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/contamination_data.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# Write the summary statistics to a separate file
write.table(summary_stats, "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/contamination_summary.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

# Print a message confirming completion
cat("Contamination data and summary statistics have been successfully written.\n")
