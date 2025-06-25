library(stringr)

# Function to extract contamination values from a file
extract.contamination <- function(filename.input) {
  input <- readLines(filename.input)
  
  # Extract normal contamination level
  normal_line <- grep("Normal sample contamination level", input, value = TRUE)
  normal_contamination <- if (length(normal_line) > 0) {
    as.numeric(sub(".*Normal sample contamination level:\\s*([0-9.]+)%.*", "\\1", normal_line))
  } else {
    NA
  }
  
  # Extract tumor contamination level
  tumor_line <- grep("Tumor sample contamination level", input, value = TRUE)
  tumor_contamination <- if (length(tumor_line) > 0) {
    as.numeric(sub(".*Tumor sample contamination level:\\s*([0-9.]+)%.*", "\\1", tumor_line))
  } else {
    NA
  }
  
  return(c(Normal = normal_contamination, Tumor = tumor_contamination))
}

# Set the input directory path
dir.input <- "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/New_TCGA_Batch_2024/output"

# Get list of sample directories (e.g., Sample001, Sample1000, etc.)
subjects <- list.files(dir.input, pattern = "^Sample\\d{3,}$", full.names = TRUE)

# Construct paths to contamination files
filenames.input <- paste0(subjects, "/", basename(subjects), "_cross_contamination.txt")

# Filter to existing files
idx <- which(file.exists(filenames.input))
subjects <- basename(subjects)[idx]
filenames.input <- filenames.input[idx]
names(filenames.input) <- subjects

# Extract contamination values
contamination_values <- t(sapply(filenames.input, extract.contamination))

# Convert to data frame
contamination_data <- data.frame(sample = rownames(contamination_values), contamination_values, row.names = NULL)

# Write contamination data
write.table(contamination_data,
            file = file.path(dir.input, "contamination_data.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Check for missing contamination values
missing_rows <- which(is.na(contamination_values[, "Normal"]) | is.na(contamination_values[, "Tumor"]))
if (length(missing_rows) > 0) {
  cat("Warning: Missing contamination values in the following samples:\n")
  print(rownames(contamination_values)[missing_rows])
}

# Summary statistics (ignoring missing values)
summary_stats <- data.frame(
  Normal_Mean   = mean(contamination_values[, "Normal"], na.rm = TRUE),
  Normal_Median = median(contamination_values[, "Normal"], na.rm = TRUE),
  Normal_SD     = sd(contamination_values[, "Normal"], na.rm = TRUE),
  Tumor_Mean    = mean(contamination_values[, "Tumor"], na.rm = TRUE),
  Tumor_Median  = median(contamination_values[, "Tumor"], na.rm = TRUE),
  Tumor_SD      = sd(contamination_values[, "Tumor"], na.rm = TRUE)
)

# Write summary statistics
write.table(summary_stats,
            file = file.path(dir.input, "contamination_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Completion message
cat("Contamination data and summary statistics have been successfully written.\n")
