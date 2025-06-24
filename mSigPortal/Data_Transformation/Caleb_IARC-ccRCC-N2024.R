
# Paper: Geographic variation of mutagenic exposures in kidney cancer genomes (PMID 38693263)


# Data: Supp Table 2-4. Denovo Constitution data. 
#       Supp Table 5.   Decomposed activities data. 
#       Supp Table 6.   Denovo activities data. 
#       Supp Table 7.   COSMIC decomposition info. 

library(tidyverse)
library(readxl)

#setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/IARC-ccRCC-N2024")
#set_wd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')


### Signature activity data--------------------------------------------------------------------------------------------------------------------------------------------

TableData <- read_xlsx(path = 'Supplementary_Table_6.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

# denovo SBS1536
sdata <- TableData %>% select(Sample,SBS_A:SBS_M)

Exposure_refdata_tmp1 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-ccRCC-N2024', Dataset = 'WGS', Cancer_Type="Kidney-ccRCC", Organ="Kidney",Signature_set_name = 'GRCh38_SBS1536_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# denovo DBS78
sdata <- TableData %>% select(Sample, DBS_A:DBS_D)

Exposure_refdata_tmp2 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-ccRCC-N2024', Dataset = 'WGS', Cancer_Type="Kidney-ccRCC", Organ="Kidney",Signature_set_name = 'GRCh38_DBS78_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# denovo ID83
sdata <- TableData %>% select(Sample, ID_A:ID_G)

Exposure_refdata_tmp3 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-ccRCC-N2024', Dataset = 'WGS', Cancer_Type="Kidney-ccRCC", Organ="Kidney",Signature_set_name = 'GRCh38_ID83_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


TableData <- read_xlsx(path = 'Supplementary_Table_5.xlsx',skip = 2,col_names = T)
colnames(TableData)[1] <- "Sample"

# decomposed SBS96 (1536 collapsed into 96 classification for COSMIC decomposition)
sdata <- TableData %>% select(Sample,SBS1:SBS40c)

Exposure_refdata_tmp4 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-ccRCC-N2024', Dataset = 'WGS', Cancer_Type="Kidney-ccRCC", Organ="Kidney",Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_SBS96_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# decomposed DBS78
sdata <- TableData %>% select(Sample,DBS2:DBS_D)

Exposure_refdata_tmp5 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-ccRCC-N2024', Dataset = 'WGS', Cancer_Type="Kidney-ccRCC", Organ="Kidney",Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_DBS78_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# decomposed ID83
sdata <- TableData %>% select(Sample,ID1:ID_C)

Exposure_refdata_tmp6 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-ccRCC-N2024', Dataset = 'WGS', Cancer_Type="Kidney-ccRCC", Organ="Kidney",Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_ID83_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


Exposure_refdata_tmp_final <- bind_rows(
  Exposure_refdata_tmp1,
  Exposure_refdata_tmp2,
  Exposure_refdata_tmp3,
  Exposure_refdata_tmp4,
  Exposure_refdata_tmp5,
  Exposure_refdata_tmp6
) %>% 
  arrange(Study,Dataset,Cancer_Type,Organ,Signature_set_name,Signature_name,Sample)


### mutational signature data----------------------------------------------------------------------------------------------------------------------------------------------- 

# denovo SBS1536
TableData <- read_xlsx(path = 'Supplementary_Table_2.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp1 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS1536',Signature_set_name = "GRCh38_SBS1536_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp1$Contribution <- as.numeric(as.character(signature_refsets_tmp1$Contribution))


# denovo DBS78
TableData <- read_xlsx(path = 'Supplementary_Table_3.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp2 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "GRCh38_DBS78_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

signature_refsets_tmp2$Contribution <- as.numeric(as.character(signature_refsets_tmp2$Contribution))


# denovo ID83
TableData <- read_xlsx(path = 'Supplementary_Table_4.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp3 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "GRCh38_ID83_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

signature_refsets_tmp3$Contribution <- as.numeric(as.character(signature_refsets_tmp3$Contribution))


# decomposed signature SBS96, extract BOTH decomposed COSMIC signatures and new decomposed signatures 
# new signature names from Supp Table 7
sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.3_Signatures_GRCh38_SBS96_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp4 <- bind_rows(
    signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.3_Signatures_GRCh38_SBS96') %>% filter(Signature_name %in% sigtmp),
  
    signature_refsets_tmp1 %>%
    mutate(Signature_name = if_else(Signature_name == 'SBS_A','SBS40b', Signature_name)) %>%
    mutate(Signature_name = if_else(Signature_name == 'SBS_B','SBS40a', Signature_name)) %>%
    mutate(Signature_name = if_else(Signature_name == 'SBS_F','SBS40c', Signature_name)) %>%
    mutate(Signature_name = if_else(Signature_name == 'SBS_I','SBS22b', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_SBS96_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp4$Signature_name))) # signature SBS22a was renamed from SBS22 but is not considered a new decomposed signature (not added to this list) 


# decomposed signature DBS78, extract BOTH decomposed COSMIC signatures and new decomposed signatures 
# new signature names from Supp Table 7
sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.3_Signatures_GRCh38_DBS78_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp5 <- bind_rows(
    signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.3_Signatures_GRCh38_DBS78') %>% filter(Signature_name %in% sigtmp),
  
    signature_refsets_tmp2 %>% 
    #mutate(Signature_name = if_else(Signature_name == 'DBS_D','DBS20', Signature_name)) %>% 
    
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_DBS78_Denovo' )

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp5$Signature_name))) 


# decomposed signature ID83, extract BOTH decomposed COSMIC signatures and new decomposed signatures 
# new signature names from Supp Table 7
sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.3_Signatures_GRCh38_ID83_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp6 <- bind_rows(
    signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.3_Signatures_GRCh37_ID83') %>% # COSMIC_v3.3_Signatures_GRCh38_ID83 doesn't existed. Replace with GRCh37.
    filter(Signature_name %in% sigtmp),
  
    signature_refsets_tmp3 %>% 
    #mutate(Signature_name = if_else(Signature_name == 'ID_C','ID23', Signature_name)) %>% 
    mutate(Signature_name = if_else(Signature_name == 'ID_F','ID1', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.3_Signatures_GRCh38_ID83_Denovo' )

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp6$Signature_name)))


# combined all signature data

signature_refsets_tmp_final <- bind_rows(
  signature_refsets_tmp1,
  signature_refsets_tmp2,
  signature_refsets_tmp3,
  signature_refsets_tmp4,
  signature_refsets_tmp5,
  signature_refsets_tmp6
) %>%
  arrange(Source,Profile,Signature_set_name,Dataset,Signature_name,MutationType)


### Mutational profile matrix data----------------------------------------------------------------------------------------------------------------------------------------------------------------------

Exposure_refdata_tmp_final %>%
  select(Study,Cancer_Type,Sample,Dataset,Signature_set_name,Signature_name,Exposure)

seqmatrix_refdata_tmp_final <- Exposure_refdata_tmp_final %>%
  left_join(signature_refsets_tmp_final %>% select(Profile,Signature_set_name,Signature_name,MutationType,Contribution)) %>%
  group_by(Study,Cancer_Type,Sample,Dataset,Profile,MutationType,Signature_set_name) %>%
  summarise(Mutations = sum(Exposure*Contribution)) %>%
  ungroup() %>%
  filter(str_starts(Signature_set_name,"GRCh")) %>%  # only select the estimation based on denovo results
  select(-Signature_set_name) %>%
  arrange(Study,Cancer_Type,Sample,Dataset,Profile,MutationType)


# Save major R objects-----------------------------------------------------------

exposures_refdata <- Exposure_refdata_tmp_final
save(exposures_refdata, file='IARC-ccRCC-N2024_WGS_exposure_refdata.RData') # study name, dataset and exposure_refdata.RData

signatures_refsets <- signature_refsets_tmp_final
save(signatures_refsets, file='IARC-ccRCC-N2024_WGS_signature_refsets.RData') # study name, dataset and signature_refsets.RData

seqmatrix_refdata <- seqmatrix_refdata_tmp_final %>% mutate(Mutations = as.integer(round(Mutations)))
save(seqmatrix_refdata, file='IARC-ccRCC-N2024_WGS_seqmatrix_refdata.RData') # study name, dataset and seqmatrix_refdata.RData


# NEXT, Collect sample level study data -----------------------------------

