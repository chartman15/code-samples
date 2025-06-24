
# Paper: Accurate and sensitive mutational signature analysis with MuSiCal (PMID 38361034)

# Data: Supp Table 1. Catalog of ID signatures. 
#       Supp Table 2. Catalog of SBS signatures. 
#       Supp Table 3. Assignment (exposure) of ID signatures in PCAWG samples. 
#       Supp Table 4. Assignment (exposure) of SBS signatures in PCAWG samples.

library(tidyverse)
library(readxl)

setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/PCAWG-PanCancer-NG2024")

# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')


### Supp Tables 3-4: Signature activity data--------------------------------------------------------------------------------------------------------------------------------------------

# denovo IDxxx
TableData <- read_xlsx(path = 'Supplementary_Table_3.xlsx',col_names = T)
colnames(TableData)[1] <- "Sample"
colnames(TableData) [27] <- "Organ"

sdata <- TableData %>% select(Sample, Organ, ID1:ID26)

Exposure_refdata_tmp1 <- sdata %>% 
  pivot_longer(cols = -c(Sample,Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'PCAWG-PanCancer-NG2024', Dataset = 'WGS', Cancer_Type="PanCancer",Signature_set_name = 'Denovo_IDxxx') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# denovo SBSxxx
TableData <- read_xlsx(path = 'Supplementary_Table_4.xlsx',col_names = T)
colnames(TableData)[1] <- "Sample"
colnames(TableData) [86] <- "Organ"

sdata <- TableData %>% select(Sample, Organ, SBS1:SBS100)

Exposure_refdata_tmp2 <- sdata %>% 
  pivot_longer(cols = -c(Sample,Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'PCAWG-PanCancer-NG2024', Dataset = 'WGS', Cancer_Type="PanCancer",Signature_set_name = 'Denovo_SBSxxx') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# combine all 
Exposure_refdata_tmp_final <- bind_rows(
  Exposure_refdata_tmp1,
  Exposure_refdata_tmp2) %>% 
  arrange(Study,Dataset,Cancer_Type,Organ,Signature_set_name,Signature_name,Sample)


### Supp Tables 1-2: Mutational signature data-------------------------------------------------------------------------------------------------------------------------------------------

# denovo IDxxx
TableData <- read_xlsx(path = 'Supplementary_Table_1.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"


signature_refsets_tmp1 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='IDxxx',Signature_set_name = "Denovo_IDxxx", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

# denovo SBSxxx
TableData <- read_xlsx(path = 'Supplementary_Table_2.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"


signature_refsets_tmp2 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBSxxx',Signature_set_name = "Denovo_SBSxxx", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

# combine all
signature_refsets_tmp_final <- bind_rows(
  signature_refsets_tmp1,
  signature_refsets_tmp2) %>% 
  arrange(Source,Profile,Signature_set_name,Dataset,Signature_name,MutationType)


### Mutational profile matrix data----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Exposure_refdata_tmp_final %>%
  select(Study,Cancer_Type,Sample,Dataset,Signature_set_name,Signature_name,Exposure)

seqmatrix_refdata_tmp_final <- Exposure_refdata_tmp_final %>%
  left_join(signature_refsets_tmp_final %>% select(Profile,Signature_set_name,Signature_name,MutationType,Contribution)) %>%
  group_by(Study,Cancer_Type,Sample,Dataset,Profile,MutationType,Signature_set_name) %>%
  summarise(Mutations = sum(Exposure*Contribution)) %>%
  ungroup() %>%
  filter(str_starts(Signature_set_name,"Denovo_")) %>%  # only select the estimation based on denovo results
  select(-Signature_set_name) %>%
  arrange(Study,Cancer_Type,Sample,Dataset,Profile,MutationType)


# Save major R objects-----------------------------------------------------------

exposures_refdata <- Exposure_refdata_tmp_final
save(exposures_refdata, file='PCAWG-PanCancer-NG2024_WGS_exposure_refdata.RData') # study name, dataset and exposure_refdata.RData

signatures_refsets <- signature_refsets_tmp_final
save(signatures_refsets, file='PCAWG-PanCancer-NG2024_WGS_signature_refsets.RData') # study name, dataset and signature_refsets.RData

seqmatrix_refdata <- seqmatrix_refdata_tmp_final
save(seqmatrix_refdata, file='PCAWG-PanCancer-NG2024_WGS_seqmatrix_refdata.RData') # study name, dataset and seqmatrix_refdata.RData


