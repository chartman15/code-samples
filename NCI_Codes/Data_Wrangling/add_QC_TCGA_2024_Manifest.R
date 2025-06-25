library(readr)
library(readxl)
library(dplyr)
library(openxlsx)

# Load files
manifest <- read_csv("New_TCGA_Manifest_Permissions.csv")
qc <- read_excel("TCGA_2024_QC_Concordance_Contamination.xlsx")

# Extract just the filename from CRAM path
manifest$filename <- basename(manifest$CRAM_File_Path)
qc$tumor_filename <- basename(qc$filename_bam_tumour)
qc$normal_filename <- basename(qc$filename_bam_normal)

# Match based on filename only (sample-level identifier)
match_tumor <- match(manifest$filename, qc$tumor_filename)
match_normal <- match(manifest$filename, qc$normal_filename)

# Prefer tumor match; fall back on normal
match_best <- ifelse(!is.na(match_tumor), match_tumor, match_normal)

# Add notes to manifest
manifest <- manifest %>%
  mutate(
    Notes_contamination = qc$Notes_contamination[match_best],
    Notes_concordance = qc$Notes_concordance[match_best]
  )

# Save output files
write.table(
  manifest,
  file = "Merged_Manifest_with_QC.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  eol = "\n"
)

write.xlsx(manifest, "Merged_Manifest_with_QC.xlsx", overwrite = TRUE)

cat("✅ Matching complete. Notes columns added based on filename.\n")
