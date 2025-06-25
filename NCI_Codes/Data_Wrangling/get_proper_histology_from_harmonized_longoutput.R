library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

# Load Excel files
dbgap <- read_excel("dbGAP_Excel_818_Sherlock_Subjects.xlsx")
wgs <- read_excel("WGS_Analyzed_samples.xlsx")
harmonized <- read_excel("sherlock_2025March_Harmonized.xlsx")

# Standardize column names to lowercase
colnames(dbgap) <- tolower(colnames(dbgap))
colnames(wgs) <- tolower(colnames(wgs))
colnames(harmonized) <- tolower(colnames(harmonized))

id_cols_to_use <- c("subject", "subject_id", "sherlock_pid")

get_histology_by_barcode <- function(barcode) {
  if (is.na(barcode) || barcode == "") return(NA)
  wgs_row <- wgs %>% filter(barcode == !!barcode)
  if (nrow(wgs_row) == 0) return(NA)
  possible_ids <- unique(na.omit(unlist(wgs_row %>% select(any_of(id_cols_to_use)))))
  if (length(possible_ids) == 0) return(NA)
  matched <- harmonized %>%
    filter(
      subject_id %in% possible_ids |
        sherlock_pid %in% possible_ids
    )
  if (nrow(matched) == 0) return(NA)
  unique_histologies <- unique(matched$histology_composite)
  if (length(unique_histologies) > 1) {
    warning(paste0("Multiple histologies found for barcode ", barcode, ": ", paste(unique_histologies, collapse = ", ")))
  }
  return(unique_histologies[1])
}

# Add tumor and normal histology columns
dbgap_with_histology <- dbgap %>%
  mutate(
    tumor_histology = sapply(tumor_barcode, get_histology_by_barcode),
    normal_histology = sapply(normal_barcode, get_histology_by_barcode)
  )

# Convert to long format: stack tumor and normal histologies with sample_type indicator
dbgap_long_histology <- dbgap_with_histology %>%
  pivot_longer(
    cols = c(tumor_histology, normal_histology),
    names_to = "sample_type",
    values_to = "histology"
  ) %>%
  select(everything())  # optional: reorder if needed

# Write the long format table to a new Excel file
write_xlsx(dbgap_long_histology, "dbgap_with_histology_long_format.xlsx")
