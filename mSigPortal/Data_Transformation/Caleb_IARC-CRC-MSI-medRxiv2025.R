# Paper: Geographic and age-related variations in mutational processes in colorectal cancer (medRxiv) 


# Data: Supp Note Table 2-6. MSI Denovo Constitution data.   
#       Supp Note Table 8-12. MSI Decomposed Activities data. 

# Note: MSI = microsatellite instability
# Note: CosName = COSMIC naming scheme; OtherName = other naming scheme  

library(tidyverse)
library(readxl)

#setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/IARC-CRC-medRxiv2025")
#set_wd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')

# Signature activity data--------------------------------------------------------------------------------------------------------------------------------------------

## SBS activities ##-------------------------------------

# CosName decomposed MSI SBS-288

TableData <- read_xlsx(path = 'Supp_Note_Table_8_Activities_MSI_Decomposed_SBS.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,SBS1:SBS93)

Exposure_refdata_tmp1 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_SBS288_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# OtherName decomposed MSI SBS-288 

TableData <- read_xlsx(path = 'Supp_Note_Table_8_Activities_MSI_Decomposed_SBS.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,SBS_I_MSI:SBS_O_MSI)

Exposure_refdata_tmp2 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_SBS288_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure) 


## ID activities ##-----------------------------------------------

# CosName decomposed MSI ID-83

TableData <- read_xlsx(path = 'Supp_Note_Table_9_Activities_MSI_Decomposed_ID.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,ID1:ID2)

Exposure_refdata_tmp3 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_ID83_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# OtherName decomposed MSI ID-83 - Not present  

## DBS activities ##---------------------------------------------------------- 

# CosName decomposed MSI DBS-78

TableData <- read_xlsx(path = 'Supp_Note_Table_10_Activities_MSI_Decomposed_DBS.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,DBS2:DBS18)

Exposure_refdata_tmp4 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_DBS78_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# OtherName decomposed MSI DBS-78

TableData <- read_xlsx(path = 'Supp_Note_Table_10_Activities_MSI_Decomposed_DBS.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,DBS_B_MSI)

Exposure_refdata_tmp5 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_DBS78_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

## CN activities ##-------------------------------------------------------------

# CosName decomposed MSI CN-68

TableData <- read_xlsx(path = 'Supp_Note_Table_11_Activities_MSI_Decomposed_CN.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,CN1:CN25)

Exposure_refdata_tmp6 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_CN68_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# OtherName decomposed MSI CN-68 - Not present 

## SV activities ##------------------------------------------------------------------

# CosName decomposed MSI SV-38 

TableData <- read_xlsx(path = 'Supp_Note_Table_12_Activities_MSI_Decomposed_SV.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,SV3:SV7)

Exposure_refdata_tmp7 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_SV38_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# OtherName decomposed MSI SV-38 

TableData <- read_xlsx(path = 'Supp_Note_Table_12_Activities_MSI_Decomposed_SV.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

sdata <- TableData %>% select(Sample,SV_D)

Exposure_refdata_tmp8 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-CRC-medRxiv2025', Dataset = 'WGS', Cancer_Type="Colon-CRC", Organ="Colon",Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_SV38_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


Exposure_refdata_tmp_final <- bind_rows(
  Exposure_refdata_tmp1,
  Exposure_refdata_tmp2,
  Exposure_refdata_tmp3,
  Exposure_refdata_tmp4,
  Exposure_refdata_tmp5,
  Exposure_refdata_tmp6,
  Exposure_refdata_tmp7,
  Exposure_refdata_tmp8
) %>% 
  arrange(Study,Dataset,Cancer_Type,Organ,Signature_set_name,Signature_name,Sample)



# Mutational signature data----------------------------------------------------------------------------------------------------------------------------------------------- 

## SBS profile ##----------------------------------------------

# denovo MSI SBS-288
TableData <- read_xlsx(path = 'Supp_Note_Table_2_MutProfile_MSI_DeNovo_SBS.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp1 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS288',Signature_set_name = "GRCh38_SBS288_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp1$Contribution <- as.numeric(as.character(signature_refsets_tmp1$Contribution))


## ID profile ##------------------------------------------------------------

# denovo MSI ID-83
TableData <- read_xlsx(path = 'Supp_Note_Table_3_MutProfile_MSI_DeNovo_ID.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp2 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "GRCh38_ID83_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp2$Contribution <- as.numeric(as.character(signature_refsets_tmp2$Contribution))

## DBS profile ##---------------------------------------------------------

# denovo MSI DBS-78
TableData <- read_xlsx(path = 'Supp_Note_Table_4_MutProfile_MSI_DeNovo_DBS.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp3 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "GRCh38_DBS78_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp3$Contribution <- as.numeric(as.character(signature_refsets_tmp3$Contribution))

## CN profile ##-------------------------------------------------

# denovo MSI CN-68
TableData <- read_xlsx(path = 'Supp_Note_Table_5_MutProfile_MSI_DeNovo_CN.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp4 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='CN68',Signature_set_name = "GRCh38_CN68_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp4$Contribution <- as.numeric(as.character(signature_refsets_tmp4$Contribution))

## SV profile ##----------------------------------------------------

# denovo MSI SV-38
TableData <- read_xlsx(path = 'Supp_Note_Table_6_MutProfile_MSI_DeNovo_SV.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp5 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SV38',Signature_set_name = "GRCh38_SV38_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp5$Contribution <- as.numeric(as.character(signature_refsets_tmp5$Contribution)) 



## SBS decomposed profile ##-----------------------------------------------------------------

# decomposed MSI signature SBS288, extract BOTH decomposed COSMIC signatures and new decomposed signatures
# new signature names from "COSMICv3.4_MSI_Decomposition"

sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh38_SBS288_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp6 <- bind_rows(
  signature_refsets %>%
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_SBS288') %>% filter(Signature_name %in% sigtmp), # v3.4 changed to v3.2 & GRCh38 changed to GRCh37
  
  signature_refsets_tmp1 %>%
    mutate(Signature_name = if_else(Signature_name == 'SBS_B_MSI','SBS44', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>%
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_SBS288_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp6$Signature_name)))


## ID decomposed profile ##--------------------------------------------------------------------

sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh38_ID83_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp7 <- bind_rows(
  signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh37_ID83') %>% filter(Signature_name %in% sigtmp), # GRCh38 changed to GRCh37
  
  signature_refsets_tmp2 %>%
    mutate(Signature_name = if_else(Signature_name == 'ID_C_MSI','ID2', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_ID83_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp7$Signature_name))) 


## DBS decomposed profile ##-------------------------------------------------------------------------

sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh38_DBS78_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp8 <- bind_rows(
  signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh38_DBS78') %>% filter(Signature_name %in% sigtmp), 
  
  signature_refsets_tmp3 %>%
    mutate(Signature_name = if_else(Signature_name == 'DBS_D_MSI','DBS8', Signature_name)) %>%
    mutate(Signature_name = if_else(Signature_name == 'DBS_F_MSI','DBS14', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>% 
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_DBS78_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp8$Signature_name))) 


## CN decomposed profile ##-----------------------------------------------------------------------------

sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh38_CN68_Denovo') %>% pull(Signature_name) %>% unique()

# signature_refsets_tmp9 <- bind_rows(
#   signature_refsets %>%
#     filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh37_CN48') %>% filter(Signature_name %in% sigtmp), # No CN68 in signature_refsets


signature_refsets_tmp9 <- signature_refsets_tmp4 %>%
  mutate(Signature_name = if_else(Signature_name == 'CN_A_MSI','CN1', Signature_name)) %>%
  mutate(Signature_name = if_else(Signature_name == 'CN_C_MSI','CN2', Signature_name)) %>%
  
  filter(Signature_name %in% sigtmp) %>%
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_CN68_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp9$Signature_name)))


## SV decomposed profile ##---------------------------------------------------------------------------

sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh38_SV38_Denovo') %>% pull(Signature_name) %>% unique()

# signature_refsets_tmp10 <- bind_rows(
#   signature_refsets %>%
#     filter(Signature_set_name == 'COSMIC_v3.4_Signatures_GRCh38_SV32') %>% filter(Signature_name %in% sigtmp), # No SV38 in signature_refsets 

signature_refsets_tmp10 <- signature_refsets_tmp5 %>%
  mutate(Signature_name = if_else(Signature_name == 'SV_B_MSI','SV_D', Signature_name)) %>%
  
  filter(Signature_name %in% sigtmp) %>%
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.4_Signatures_GRCh38_SV38_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp10$Signature_name)))

# combined all signature data

signature_refsets_tmp_final <- bind_rows(
  signature_refsets_tmp1,
  signature_refsets_tmp2,
  signature_refsets_tmp3,
  signature_refsets_tmp4,
  signature_refsets_tmp5,
  signature_refsets_tmp6,
  signature_refsets_tmp7,
  signature_refsets_tmp8,
  signature_refsets_tmp9,
  signature_refsets_tmp10
  
) %>%
  arrange(Source,Profile,Signature_set_name,Dataset,Signature_name,MutationType)

# Mutational profile matrix data----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Exposure_refdata_tmp_final %>%
  select(Study,Cancer_Type,Sample,Dataset,Signature_set_name,Signature_name,Exposure)

seqmatrix_refdata_tmp_final <- Exposure_refdata_tmp_final %>%
  left_join(signature_refsets_tmp_final %>% select(Profile,Signature_set_name,Signature_name,MutationType,Contribution)) %>%
  group_by(Study,Cancer_Type,Sample,Dataset,Profile,MutationType,Signature_set_name) %>%
  summarise(Mutations = sum(Exposure*Contribution)) %>%
  ungroup() %>%
  filter(str_starts(Signature_set_name,"COSMIC")) %>% 
  select(-Signature_set_name) %>%
  arrange(Study,Cancer_Type,Sample,Dataset,Profile,MutationType)


# Save major R objects-----------------------------------------------------------

exposures_refdata <- Exposure_refdata_tmp_final
save(exposures_refdata, file='IARC-CRC-MSI-medRxiv2025_WGS_exposure_refdata.RData') # study name, dataset and exposure_refdata.RData

signatures_refsets <- signature_refsets_tmp_final
save(signatures_refsets, file='IARC-CRC-MSI-medRxiv2025_WGS_signature_refsets.RData') # study name, dataset and signature_refsets.RData

seqmatrix_refdata <- seqmatrix_refdata_tmp_final %>% mutate(Mutations = as.integer(round(Mutations)))
save(seqmatrix_refdata, file='IARC-CRC-MSI-medRxiv2025_WGS_seqmatrix_refdata.RData') # study name, dataset and seqmatrix_refdata.RData

