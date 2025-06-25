# library(stringr)
# 
# extract.concordance <- function(
#     filename.input
# ) {
#     input <- read.delim(filename.input, header=F, stringsAsFactors=F, sep=" ")
#     con <- input[1,2]
#     return(con)
# }
# #dir.input <- "/data/hoangph/Conpair/output/Sherlock"
# dir.input <- "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/output/Sherlock" 
# subjects <- grep("-T", list.files(dir.input), value=T)
# filenames.input <- paste0(dir.input, "/", subjects, "/", subjects, "_homozygous_concordance.txt")
# idx <- which(file.exists(filenames.input))
# subjects <- subjects[idx]
# filenames.input <- filenames.input[idx]
# names(filenames.input) <- subjects
# out <- sapply(subjects, function(subject) extract.concordance(filenames.input[subject])) 
# output <- data.frame(sample=subjects,concordance=out)
# write.table(output,"/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/Conpair_tumour_normal_concordance.txt",col.names=T,row.names=F,quote=F,sep='\t')


library(stringr)

# Define the function to extract concordance from a file
extract.concordance <- function(filename.input) {
  input <- read.delim(filename.input, header = FALSE, stringsAsFactors = FALSE, sep = " ")
  con <- input[1, 2]  # Assuming the concordance is always in the second column of the first row
  return(con)
}

# Set the directory input path
dir.input <- "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/output/Sherlock"

# Get the list of sample directories (Sample001, Sample002, etc.)
subjects <- list.files(dir.input, pattern = "^Sample\\d{3}$", full.names = TRUE)

# Create the corresponding file paths for concordance files
filenames.input <- paste0(subjects, "/", basename(subjects), "_homozygous_concordance.txt")

# Check which concordance files exist
idx <- which(file.exists(filenames.input))

# Filter the subjects and filenames based on existing files
subjects <- basename(subjects)[idx]
filenames.input <- filenames.input[idx]

# Create a named vector for the concordance values
names(filenames.input) <- subjects

# Extract concordance values for each subject using sapply
out <- sapply(filenames.input, extract.concordance)

# Create a data frame with sample names and corresponding concordance values
output <- data.frame(sample = names(out), concordance = out)

# Write the output to a file
write.table(output, "/data/Sherlock_Lung/CalebHartman/Phuc_Conpair/Conpair_tumour_normal_concordance.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
