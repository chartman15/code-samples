# Paper: Quantitative and qualitative mutational impact of ionizing radiation on normal cells (PMID 38359788)

# Data: Supp Table 3-4.  SBS Constitution. 
#       Supp Table 7-8.  ID Constitution. 
#       Supp Table 5-6.  SBS Activities.
#       Supp Table 9-10. ID Activities. 


library(tidyverse)
library(readxl)

setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/KoNA-PanCancer-CG2024")

# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')

### Signature activity data--------------------------------------------------------------------------------------------------------------------------------------------


# denovo SBS96 - Mouse
TableData <- read_xlsx(path = 'Supp_Table_5.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,`mSBS-A`:`mSBS-D`)

Exposure_refdata_tmp1 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'KoNA-PanCancer-CG2024', Dataset = 'WGS', Cancer_Type="PanCancer", Organ="Various",Signature_set_name = 'Denovo_SBS96_Mouse') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# denovo sbs96 - Human
TableData <- read_xlsx(path = 'Supp_Table_6.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,`hSBS-A`:`hSBS-G`)

Exposure_refdata_tmp2 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'KoNA-PanCancer-CG2024', Dataset = 'WGS', Cancer_Type="PanCancer", Organ="Various",Signature_set_name = 'Denovo_SBS96_Human') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# denovo ID83 - Mouse
TableData <- read_xlsx(path = 'Supp_Table_9.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,`mID-A`:`mID-D`)

Exposure_refdata_tmp3 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'KoNA-PanCancer-CG2024', Dataset = 'WGS', Cancer_Type="PanCancer", Organ="Various",Signature_set_name = 'Denovo_ID83_Mouse') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# denovo ID83 - Human 
TableData <- read_xlsx(path = 'Supp_Table_10.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,`hID-A`:`hID-D`)

Exposure_refdata_tmp4 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'KoNA-PanCancer-CG2024', Dataset = 'WGS', Cancer_Type="PanCancer", Organ="Various",Signature_set_name = 'Denovo_ID83_Human') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# Combine all activity data 

Exposure_refdata_tmp_final <- bind_rows(
  Exposure_refdata_tmp1,
  Exposure_refdata_tmp2,
  Exposure_refdata_tmp3,
  Exposure_refdata_tmp4,
) %>% 
  arrange(Study,Dataset,Cancer_Type,Organ,Signature_set_name,Signature_name,Sample)


### mutational signature data----------------------------------------------------------------------------------------------------------------------------------------------- 

# denovo SBS96 - Mouse
TableData <- read_xlsx(path = 'Supp_Table_3.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
colnames(TableData)[2] <- "MutationContext"

signature_refsets_tmp1 <- TableData %>% 
  pivot_longer(cols = -c(MutationType,MutationContext), names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS96',Signature_set_name = "Denovo_SBS96_Mouse", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, MutationContext, Contribution)

# denovo SBS96 - Human 
TableData <- read_xlsx(path = 'Supp_Table_4.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
colnames(TableData)[2] <- "MutationContext"

signature_refsets_tmp2 <- TableData %>% 
  pivot_longer(cols = -c(MutationType,MutationContext), names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS96',Signature_set_name = "Denovo_SBS96_Human", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, MutationContext, Contribution)

# denovo ID83 - Mouse
TableData <- read_xlsx(path = 'Supp_Table_7.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
colnames(TableData)[2] <- "MutationContext"

signature_refsets_tmp3 <- TableData %>% 
  pivot_longer(cols = -c(MutationType,MutationContext), names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "Denovo_ID83_Mouse", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, MutationContext, Contribution)

# denovo ID83 - Human 
TableData <- read_xlsx(path = 'Supp_Table_8.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
colnames(TableData)[2] <- "MutationContext"

signature_refsets_tmp4 <- TableData %>% 
  pivot_longer(cols = -c(MutationType,MutationContext), names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "Denovo_ID83_Human", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, MutationContext, Contribution)

# combined all signature data

signature_refsets_tmp_final <- bind_rows(
  signature_refsets_tmp1,
  signature_refsets_tmp2,
  signature_refsets_tmp3,
  signature_refsets_tmp4
) %>% 
  arrange(Source,Profile,Signature_set_name,Dataset,Signature_name,MutationType,MutationContext)

### Mutational profile matrix data----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Exposure_refdata_tmp_final %>%
  select(Study,Cancer_Type,Sample,Dataset,Signature_set_name,Signature_name,Exposure)

seqmatrix_refdata_tmp_final <- Exposure_refdata_tmp_final %>%
  left_join(signature_refsets_tmp_final %>% select(Profile,Signature_set_name,Signature_name,MutationType,MutationContext,Contribution)) %>%
  group_by(Study,Cancer_Type,Sample,Dataset,Profile,MutationType,MutationContext,Signature_set_name) %>%
  summarise(Mutations = sum(Exposure*Contribution)) %>%
  ungroup() %>%
  filter(str_starts(Signature_set_name,"Denovo_")) %>%  # only select the estimation based on denovo results
  select(-Signature_set_name) %>%
  arrange(Study,Cancer_Type,Sample,Dataset,Profile,MutationType,MutationContext)


# Save major R objects-----------------------------------------------------------

exposures_refdata <- Exposure_refdata_tmp_final
save(exposures_refdata, file='KoNA-PanCancer-CG2024_WGS_exposure_refdata.RData') # study name, dataset and exposure_refdata.RData

signatures_refsets <- signature_refsets_tmp_final
save(signatures_refsets, file='KoNA-PanCancer-CG2024_WGS_signature_refsets.RData') # study name, dataset and signature_refsets.RData

seqmatrix_refdata <- seqmatrix_refdata_tmp_final
save(seqmatrix_refdata, file='KoNA-PanCancer-CG2024_WGS_seqmatrix_refdata.RData') # study name, dataset and seqmatrix_refdata.RData


# NEXT, Collect sample level study data -----------------------------------



