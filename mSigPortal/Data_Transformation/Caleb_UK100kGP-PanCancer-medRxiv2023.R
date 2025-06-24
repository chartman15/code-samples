
# Paper: Comprehensive repertoire of the chromosomal alteration and mutational signatures across 16 cancer types from 10,983 cancer patients (medRxiv)

# Data: Supp Table 4a-4e: Signature activities. 
#       Supp Table 5a-5e: Mutational Signature Data (Constitution).
#       Supp Table 3: Description of signatures extracted in this study.
#       Supp Table 11: Software and data used in this study and SigProfilerExtractor parameters used for signature extraction in each mutation type.

library(tidyverse)
library(readxl)

setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/UK100kGP-PanCancer-medRxiv2023")

# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')


### Supp Table 4a-4e: Signature activity data--------------------------------------------------------------------------------------------------------------------------------------------

# decomposed SBS288 - COSMIC
TableData <- read_xlsx(path = 'Supplementary_Table_4a.xlsx', col_names = T)
colnames(TableData)[1] <- "Sample"
colnames(TableData)[2] <- "Organ"

sdata <- TableData %>% select(Sample,Organ,SBS1:SBS94)

Exposure_refdata_tmp1 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_SBS288') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed SBS96 - Degasperi et al. 2022 (SBS105/141/129)
sdata <- TableData %>% select(Sample,Organ,SBS95:SBS97)

Exposure_refdata_tmp2 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Cancer_Reference_Signatures_2022_GRCh37_SBS96') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed DBS78 - COSMIC
TableData <- read_xlsx(path = 'Supplementary_Table_4b.xlsx', col_names = T)
colnames(TableData)[1] <- "Sample"
colnames(TableData)[2] <- "Organ"

sdata <- TableData %>% select(Sample, Organ, DBS1:DBS11)

Exposure_refdata_tmp3 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_DBS78') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed DBS78 - Degasperi et al. 2022 (DBS28/17/27/12/26/29)
sdata <- TableData %>% select(Sample, Organ, DBS12:DBS17)

Exposure_refdata_tmp4 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Cancer_Reference_Signatures_2022_GRCh37_DBS78') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed DBS78 - Novel or Not linearly independent 
sdata <- TableData %>% select(Sample, Organ, DBS18:DBS19)

Exposure_refdata_tmp5 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Novel_DBS78') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed ID83 - COSMIC
TableData <- read_xlsx(path = 'Supplementary_Table_4c.xlsx', col_names = T)
colnames(TableData)[1] <- "Sample"
colnames(TableData)[2] <- "Organ"

sdata <- TableData %>% select(Sample, Organ, ID1:ID18)

Exposure_refdata_tmp6 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_ID83') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed ID83 - Novel
sdata <- TableData %>% select(Sample, Organ, ID19:ID22)

Exposure_refdata_tmp7 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Novel_ID83') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed CNV48 - COSMIC
TableData <- read_xlsx(path = 'Supplementary_Table_4d.xlsx', col_names = T)
colnames(TableData)[1] <- "Sample"
colnames(TableData)[2] <- "Organ"

sdata <- TableData %>% select(Sample, Organ, CN1:CN24)

Exposure_refdata_tmp8 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh37_CN48') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed CNV48 - Novel
sdata <- TableData %>% select(Sample, Organ, CN25)

Exposure_refdata_tmp9 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Novel_CN48') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed SV32 - Nik-Zainal et al. 2016 (RS1/2/3/4/5/6)
TableData <- read_xlsx(path = 'Supplementary_Table_4e.xlsx', col_names = T)
colnames(TableData)[1] <- "Sample"
colnames(TableData)[2] <- "Organ"

sdata <- TableData %>% select(Sample, Organ, SV1:SV6)

Exposure_refdata_tmp10 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Organ-specific_Cancer_Signatures_GRCh37_RS32') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed SV32 - Degasperi et al. 2020 (RefSigR7/10/1/8/6b)
sdata <- TableData %>% select(Sample, Organ, SV7, SV8, SV10, SV11, SV13)

Exposure_refdata_tmp11 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Cancer_Reference_Signatures_GRCh37_RS32') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed SV32 - Novel or Not linearly independent 
sdata <- TableData %>% select(Sample, Organ, SV9, SV12)

Exposure_refdata_tmp12 <- sdata %>% 
  pivot_longer(cols = c(-Sample,-Organ),names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100kGP-PanCancer-medRxiv2023', Dataset = 'WGS', Cancer_Type="PanCancer", Signature_set_name = 'Novel_RS32') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# combine all 
Exposure_refdata_tmp_final <- bind_rows(
  Exposure_refdata_tmp1,
  Exposure_refdata_tmp2,
  Exposure_refdata_tmp3,
  Exposure_refdata_tmp4,
  Exposure_refdata_tmp5,
  Exposure_refdata_tmp6,
  Exposure_refdata_tmp7,
  Exposure_refdata_tmp8,
  Exposure_refdata_tmp9,
  Exposure_refdata_tmp10,
  Exposure_refdata_tmp11,
  Exposure_refdata_tmp12
) %>% 
  arrange(Study,Dataset,Cancer_Type,Organ,Signature_set_name,Signature_name,Sample)


### Supp Table 5a-5e: mutational signature data----------------------------------------------------------------------------------------------------------------------------------------------- 

# decomposed SBS
TableData <- read_xlsx(path = 'Supplementary_Table_5a.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,SBS1:SBS94)

signature_refsets_tmp1 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS288',Signature_set_name = "COSMIC_v3.2_Signatures_GRCh38_SBS288", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

TableData <- read_xlsx(path = 'Supplementary_Table_5a.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,SBS95:SBS97)

signature_refsets_tmp2 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS96',Signature_set_name = "Cancer_Reference_Signatures_2022_GRCh37_SBS96", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

signature_refsets_tmpa <- bind_rows(signature_refsets_tmp1,signature_refsets_tmp2)


# decomposed DBS
TableData <- read_xlsx(path = 'Supplementary_Table_5b.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,DBS1:DBS11)

signature_refsets_tmp3 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "COSMIC_v3.3_Signatures_GRCh38_DBS78", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

TableData <- read_xlsx(path = 'Supplementary_Table_5b.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,DBS12:DBS17)

signature_refsets_tmp4 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "Cancer_Reference_Signatures_2022_GRCh37_DBS78", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

TableData <- read_xlsx(path = 'Supplementary_Table_5b.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,DBS18:DBS19)

signature_refsets_tmp5 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "Novel_DBS78", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

signature_refsets_tmpb <- bind_rows(signature_refsets_tmp3,signature_refsets_tmp4,signature_refsets_tmp5)


# decomposed ID
TableData <- read_xlsx(path = 'Supplementary_Table_5c.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,ID1:ID18)

signature_refsets_tmp6 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "COSMIC_v3.3_Signatures_GRCh38_ID83", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

TableData <- read_xlsx(path = 'Supplementary_Table_5c.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,ID19:ID22)

signature_refsets_tmp7 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "Novel_ID83", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

signature_refsets_tmpc <- bind_rows(signature_refsets_tmp6,signature_refsets_tmp7)


# decomposed CNV
TableData <- read_xlsx(path = 'Supplementary_Table_5d.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,CN1:CN24)

signature_refsets_tmp8 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='CN48',Signature_set_name = "COSMIC_v3.3_Signatures_GRCh37_CN48", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

TableData <- read_xlsx(path = 'Supplementary_Table_5d.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,CN25)

signature_refsets_tmp9 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='CN48',Signature_set_name = "Novel_CN48", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

signature_refsets_tmpd <- bind_rows(signature_refsets_tmp8,signature_refsets_tmp9)


# decomposed SV
TableData <- read_xlsx(path = 'Supplementary_Table_5e.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType,SV1:SV6)

signature_refsets_tmp10 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='RS32',Signature_set_name = "Organ-specific_Cancer_Signatures_GRCh37_RS32", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

TableData <- read_xlsx(path = 'Supplementary_Table_5e.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType, SV7, SV8, SV10, SV11, SV13)

signature_refsets_tmp11 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='RS32',Signature_set_name = "Cancer_Reference_Signatures_GRCh37_RS32", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

TableData <- read_xlsx(path = 'Supplementary_Table_5e.xlsx',col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"
sdata <- TableData %>% select(MutationType, SV9, SV12)

signature_refsets_tmp12 <- sdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='RS32',Signature_set_name = "Novel_RS32", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

signature_refsets_tmpe <- bind_rows(signature_refsets_tmp10,signature_refsets_tmp11,signature_refsets_tmp12)

# combine all 
signature_refsets_tmp_final <- bind_rows(
  signature_refsets_tmpa,
  signature_refsets_tmpb,
  signature_refsets_tmpc,
  signature_refsets_tmpd,
  signature_refsets_tmpe
) %>% 
  arrange(Source,Profile,Signature_set_name,Dataset,Signature_name,MutationType)



# Save major R objects-----------------------------------------------------------

exposures_refdata <- Exposure_refdata_tmp_final
save(exposures_refdata, file='UK100kGP-PanCancer-medRxiv2023_WGS_exposure_refdata.RData') # study name, dataset and exposure_refdata.RData

signatures_refsets <- signature_refsets_tmp_final
save(signatures_refsets, file='UK100kGP-PanCancer-medRxiv2023_WGS_signature_refsets.RData') # study name, dataset and signature_refsets.RData


