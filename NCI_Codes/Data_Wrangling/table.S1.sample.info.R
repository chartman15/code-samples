library(readxl)
library(dplyr)
library(writexl)

# Step 1: Load Excel files
table_s1 <- read_excel("table.S1.sample.info.xlsx")
wgs_current <- read_excel("All_WGS_Current.xlsx")
sherlock <- read_excel("sherlock_2025March.xlsx")

# Step 2: Normalize ID columns for consistent matching
wgs_current <- wgs_current %>%
  mutate(
    Subject_ID = toupper(Subject),
    Sherlock_PID_upper = toupper(sherlock_pid),
    Sample_ID = Barcode
  )

sherlock <- sherlock %>%
  mutate(
    SUBJECT_ID_upper = toupper(SUBJECT_ID),
    Sherlock_PID_upper = toupper(Sherlock_PID)
  )

# Step 3: Join on Subject_ID (Subject ↔ SUBJECT_ID)
join1 <- wgs_current %>%
  left_join(
    sherlock %>% select(SUBJECT_ID_upper, CURRENT_LAST_RES),
    by = c("Subject_ID" = "SUBJECT_ID_upper")
  ) %>%
  rename(CURRENT_LAST_RES_subj = CURRENT_LAST_RES)

# Step 4: Join on Sherlock_PID (sherlock_pid ↔ Sherlock_PID)
join2 <- join1 %>%
  left_join(
    sherlock %>% select(Sherlock_PID_upper, CURRENT_LAST_RES),
    by = "Sherlock_PID_upper"
  ) %>%
  rename(CURRENT_LAST_RES_pid = CURRENT_LAST_RES)

# Step 5: Prioritize CURRENT_LAST_RES from Subject first, then fallback to PID
wgs_enriched <- join2 %>%
  mutate(
    CURRENT_LAST_RES = coalesce(CURRENT_LAST_RES_subj, CURRENT_LAST_RES_pid)
  ) %>%
  select(Sample_ID, Subject_ID, Study, CURRENT_LAST_RES) %>%
  rename(WGS_data_source_fill = Study)

# Step 6: Deduplicate — ensure only 1 row per Sample_ID
wgs_enriched_dedup <- wgs_enriched %>%
  group_by(Sample_ID) %>%
  slice(1) %>%  # pick the first match per Sample_ID
  ungroup()

# Step 7: Join with table_s1
merged_final <- table_s1 %>%
  left_join(wgs_enriched_dedup, by = "Sample_ID")

# Step 8: Fill in missing fields only if NA
table_s1_filled <- merged_final %>%
  mutate(
    WGS_data_source = coalesce(WGS_data_source, WGS_data_source_fill),
    Country_state = coalesce(Country_state, CURRENT_LAST_RES)
  ) %>%
  select(all_of(names(table_s1)))  # Keep original structure

# Step 9: Write to Excel
write_xlsx(table_s1_filled, "table.S1.sample.info.updated.xlsx")
