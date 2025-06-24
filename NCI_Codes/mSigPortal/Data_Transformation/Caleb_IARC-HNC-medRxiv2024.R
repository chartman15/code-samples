# Paper: The Complexity of Tobacco Smoke-Induced Mutagenesis in Head and Neck Cancer (medRxiv)

# Data: Supp Table 4-6. Denovo Constitution data (SBS, DBS, ID).  
#       Supp Table 7.   COSMIC decomposition info (SBS, DBS, ID). 
#       Supp Table 8.   Denovo Activities data (SBS, DBS, ID).
#       Supp Table 9.   Decomposed Activities data (SBS, DBS, ID). 
#       Supp Table 17.  Denovo Constitution data (CNV).
#       SN Table 4.     COSMIC decomposition info (CNV).
#       SN Table 5.     Denovo Activities data (CNV).
#       SN Table 6.     Decomposed Activities data (CNV). 


library(tidyverse)
library(readxl)

setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/IARC-HNC-medRxiv2024")

# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')


### Signature activity data--------------------------------------------------------------------------------------------------------------------------------------------

TableData <- read_xlsx(path = 'Supp_Table_8.xlsx',skip = 1,col_names = T)
colnames(TableData)[1] <- "Sample"

# denovo SBS1536
sdata <- TableData %>% select(Sample,SBS_A:SBS_O)

Exposure_refdata_tmp1 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'Denovo_SBS1536') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# denovo DBS78
sdata <- TableData %>% select(Sample,DBS_A:DBS_D)

Exposure_refdata_tmp2 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'Denovo_DBS78') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# denovo ID83
sdata <- TableData %>% select(Sample,ID_A:ID_G)

Exposure_refdata_tmp3 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'Denovo_ID83') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

TableData <- read_xlsx(path = 'SN_Table_5.xlsx',skip = 1,col_names = T)
colnames(TableData)[1] <- "Sample"

# denovo CNV48
sdata <- TableData %>% select(Sample,CN1:CN_G)

Exposure_refdata_tmp4 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'Denovo_CN48') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

TableData <- read_xlsx(path = 'Supp_Table_9.xlsx',skip = 1,col_names = T)
colnames(TableData)[1] <- "Sample"

# decomposed SBS96
sdata <- TableData %>% select(Sample,SBS1:SBS_L)

Exposure_refdata_tmp5 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_SBS96_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed DBS78
sdata <- TableData %>% select(Sample,DBS1:DBS_D)

Exposure_refdata_tmp6 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_DBS78_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decomposed ID83 
sdata <- TableData %>% select(Sample,ID1:ID14)

Exposure_refdata_tmp7 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_ID83_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

TableData <- read_xlsx(path = 'SN_Table_6.xlsx',skip = 1,col_names = T)
colnames(TableData)[1] <- "Sample"

# decomposed CNV48
sdata <- TableData %>% select(Sample,CN1:CN_G)

Exposure_refdata_tmp8 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'IARC-HNC-medRxiv2024', Dataset = 'WGS', Cancer_Type="HNC", Organ="Head_Neck",Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_CN48_Denovo') %>% 
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


### mutational signature data----------------------------------------------------------------------------------------------------------------------------------------------- 

# denovo SBS1536
TableData <- read_xlsx(path = 'Supp_Table_4.xlsx',skip = 1,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp1 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS1536',Signature_set_name = "Denovo_SBS1536", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp1$Contribution <- as.numeric(as.character(signature_refsets_tmp1$Contribution))

# denovo DBS78
TableData <- read_xlsx(path = 'Supp_Table_5.xlsx',skip = 1,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp2 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "Denovo_DBS78", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp2$Contribution <- as.numeric(as.character(signature_refsets_tmp2$Contribution))

# denovo ID83
TableData <- read_xlsx(path = 'Supp_Table_6.xlsx',skip = 1,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp3 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "Denovo_ID83", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp3$Contribution <- as.numeric(as.character(signature_refsets_tmp3$Contribution))

# denovo CNV48
TableData <- read_xlsx(path = 'Supp_Table_17.xlsx',skip = 1,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp4 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='CN48',Signature_set_name = "Denovo_CN48", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution) 

signature_refsets_tmp4$Contribution <- as.numeric(as.character(signature_refsets_tmp4$Contribution))


# decomposed signature SBS96, extract BOTH decomposed COSMIC signatures and new decomposed signatures 
# new signature names from Supp Table 7
sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_SBS96_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp5 <- bind_rows(
    signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_SBS96') %>% filter(Signature_name %in% sigtmp),
    signature_refsets_tmp1 %>%
      
    mutate(Signature_name = if_else(Signature_name == 'SBS1536I','SBS_I', Signature_name)) %>%
    mutate(Signature_name = if_else(Signature_name == 'SBS1536L','SBS_L', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_SBS96_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp5$Signature_name)))


# decomposed signature DBS78, extract BOTH decomposed COSMIC signatures and new decomposed signatures 
# new signature names from Supp Table 7
sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_DBS78_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp6 <- bind_rows(
    signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_DBS78') %>% filter(Signature_name %in% sigtmp),
    signature_refsets_tmp2 %>%
    
    mutate(Signature_name = if_else(Signature_name == 'DBS78D','DBS_D', Signature_name)) %>%
    #mutate(Signature_name = if_else(Signature_name == 'DBS78B','DBS1', Signature_name)) %>%
    #mutate(Signature_name = if_else(Signature_name == 'DBS78C','DBS4', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_DBS78_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp6$Signature_name)))


# decomposed signature ID83, extract BOTH decomposed COSMIC signatures and new decomposed signatures 
# new signature names from Supp Table 7
sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_ID83_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp7 <- bind_rows(
    signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_ID83') %>% filter(Signature_name %in% sigtmp), # change to GRCh37
    signature_refsets_tmp3 %>%
    
    #mutate(Signature_name = if_else(Signature_name == 'ID83A','ID3', Signature_name)) %>%
    
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_ID83_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp7$Signature_name)))


# decomposed signature CNV48, extract BOTH decomposed COSMIC signatures and new decomposed signatures 
# new signature names from SN Table 4
sigtmp <- Exposure_refdata_tmp_final %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_CN48_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp8 <- bind_rows(
    signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_CN48') %>% filter(Signature_name %in% sigtmp), 
    signature_refsets_tmp4 %>%
      
    filter(Signature_name %in% sigtmp)) %>% 
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_CN48_Denovo')

# check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp8$Signature_name))) # only collects CN_G because signature names are the same for decomposed and denovo CNVs 

# combined all signature data

signature_refsets_tmp_final <- bind_rows(
  signature_refsets_tmp1,
  signature_refsets_tmp2,
  signature_refsets_tmp3,
  signature_refsets_tmp4,
  signature_refsets_tmp5,
  signature_refsets_tmp6,
  signature_refsets_tmp7,
  signature_refsets_tmp8
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
  filter(str_starts(Signature_set_name,"Denovo_")) %>%  # only select the estimation based on denovo results
  select(-Signature_set_name) %>%
  arrange(Study,Cancer_Type,Sample,Dataset,Profile,MutationType)


# Save major R objects-----------------------------------------------------------

exposures_refdata <- Exposure_refdata_tmp_final
save(exposures_refdata, file='IARC-HNC-medRxiv2024_WGS_exposure_refdata.RData') # study name, dataset and exposure_refdata.RData

signatures_refsets <- signature_refsets_tmp_final
save(signatures_refsets, file='IARC-HNC-medRxiv2024_WGS_signature_refsets.RData') # study name, dataset and signature_refsets.RData

seqmatrix_refdata <- seqmatrix_refdata_tmp_final
save(seqmatrix_refdata, file='IARC-HNC-medRxiv2024_WGS_seqmatrix_refdata.RData') # study name, dataset and seqmatrix_refdata.RData


# NEXT, Collect sample level study data -----------------------------------


