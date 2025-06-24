

library(tidyverse)
library(readxl)

setwd("C:/Users/hartmancas/Desktop/mSigPortal_Data_Cleaning")

# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')


# collecting signature activity data --------------------------------------

tdata <- read_xlsx(path = 'Supplementary_Table_13.xlsx',skip = 4,col_names = T)
colnames(tdata)[1] <- "Sample"

# denovo SBS96
sdata <- tdata %>% select(Sample,SBS96A:SBS96AA)
exposure_refdata_tmp1 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100k-CRC-N2023', Dataset = 'WGS', Cancer_Type="Colorectal_Carcinoma", Organ="Colon",Signature_set_name = 'Denovo_SBS96') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# denovo DBS78
sdata <- tdata %>% select(Sample,DBS78A:DBS78H)
exposure_refdata_tmp2 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100k-CRC-N2023', Dataset = 'WGS', Cancer_Type="Colorectal_Carcinoma", Organ="Colon",Signature_set_name = 'Denovo_DBS78') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# denovo ID83
sdata <- tdata %>% select(Sample,ID83A:ID83K)
exposure_refdata_tmp3 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100k-CRC-N2023', Dataset = 'WGS', Cancer_Type="Colorectal_Carcinoma", Organ="Colon",Signature_set_name = 'Denovo_ID83') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decompsite SBS96
sdata <- tdata %>% select(Sample,SBS1:`SBS-CRC2`)
exposure_refdata_tmp4 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100k-CRC-N2023', Dataset = 'WGS', Cancer_Type="Colorectal_Carcinoma", Organ="Colon",Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_SBS96_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


# decompsite DBS78
sdata <- tdata %>% select(Sample,DBS2:`DBS-CRC5`)
exposure_refdata_tmp5 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100k-CRC-N2023', Dataset = 'WGS', Cancer_Type="Colorectal_Carcinoma", Organ="Colon",Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_DBS78_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)

# decompsite ID83
sdata <- tdata %>% select(Sample,ID1:`ID-CRC2`)
exposure_refdata_tmp6 <- sdata %>% 
  pivot_longer(cols = -Sample,names_to = "Signature_name",values_to = 'Exposure') %>% 
  mutate(Study = 'UK100k-CRC-N2023', Dataset = 'WGS', Cancer_Type="Colorectal_Carcinoma", Organ="Colon",Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_ID83_Denovo') %>% 
  select(Study,Dataset,Cancer_Type,Organ,Sample,Signature_set_name,Signature_name,Exposure)


exposure_refdata_tmp <- bind_rows(
  exposure_refdata_tmp1,
  exposure_refdata_tmp2,
  exposure_refdata_tmp3,
  exposure_refdata_tmp4,
  exposure_refdata_tmp5,
  exposure_refdata_tmp6
) %>% 
  arrange(Study,Dataset,Cancer_Type,Organ,Signature_set_name,Signature_name,Sample)


# collecting mutational signature data ---------------------------------------------------

#signature denovo SBS96
tdata <- read_xlsx(path = 'Supplementary_Table_11.xlsx',skip = 3,col_names = T,range = "A4:AB100",.name_repair = 'check_unique')
colnames(tdata)[1] <- "MutationType"

signature_refsets_tmp1 <- tdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS96',Signature_set_name = "Denovo_SBS96", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source,Profile,Signature_set_name,Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

#signature denovo DBS78
tdata <- read_xlsx(path = 'Supplementary_Table_11.xlsx',skip = 3,col_names = T,range = "BD4:BL82",.name_repair = 'check_unique')
colnames(tdata)[1] <- "MutationType"

signature_refsets_tmp2 <- tdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "Denovo_DBS78", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source,Profile,Signature_set_name,Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)


#signature denovo ID83
tdata <- read_xlsx(path = 'Supplementary_Table_11.xlsx',skip = 3,col_names = T,range = "BU4:CF87",.name_repair = 'check_unique')
colnames(tdata)[1] <- "MutationType"

signature_refsets_tmp3 <- tdata %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "Denovo_ID83", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source,Profile,Signature_set_name,Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)


# decomposition signature SBS96, only extract if not collected in signature_refsets
# new signature name from the Supplementary Table 12

sigtmp <- exposure_refdata_tmp %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_SBS96_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp4 <- bind_rows(
  signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_SBS96') %>% 
    filter(Signature_name %in% sigtmp),
  
  signature_refsets_tmp1 %>% 
    mutate(Signature_name = if_else(Signature_name == 'SBS96J',"SBS-CRC1",Signature_name)) %>% 
    mutate(Signature_name = if_else(Signature_name == 'SBS96AA',"SBS-CRC2",Signature_name)) %>% 
    filter(Signature_name %in% sigtmp)
) %>% 
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_SBS96_Denovo' )

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp4$Signature_name)))


# decomposition signature DBS78, only extract if not collected in signature_refsets
# new signature name from the Supplementary Table 12

sigtmp <- exposure_refdata_tmp %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_DBS78_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp5 <- bind_rows(
  signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_DBS78') %>% 
    filter(Signature_name %in% sigtmp),
  
  signature_refsets_tmp2 %>% 
    mutate(Signature_name = if_else(Signature_name == 'DBS78A',"DBS-CRC1", Signature_name)) %>% 
    mutate(Signature_name = if_else(Signature_name == 'DBS78B',"DBS-CRC2", Signature_name)) %>% 
    mutate(Signature_name = if_else(Signature_name == 'DBS78C',"DBS-CRC3", Signature_name)) %>% 
    mutate(Signature_name = if_else(Signature_name == 'DBS78F',"DBS-CRC4", Signature_name)) %>% 
    mutate(Signature_name = if_else(Signature_name == 'DBS78G',"DBS-CRC5", Signature_name)) %>% 
    filter(Signature_name %in% sigtmp)
) %>% 
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_DBS78_Denovo' )

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp5$Signature_name)))


# decomposition signature ID83, only extract if not collected in signature_refsets
# new signature name from the Supplementary Table 12

sigtmp <- exposure_refdata_tmp %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh38_ID83_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp6 <- bind_rows(
  signature_refsets %>% 
    filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_ID83') %>% # COSMIC_v3.2_Signatures_GRCh38_ID83 doesn't existed. replace with GRCh37
    filter(Signature_name %in% sigtmp),
  
  signature_refsets_tmp3 %>% 
    mutate(Signature_name = if_else(Signature_name == 'ID83D',"ID-CRC1", Signature_name)) %>% 
    mutate(Signature_name = if_else(Signature_name == 'ID83K',"ID-CRC2", Signature_name)) %>% 
    filter(Signature_name %in% sigtmp)
) %>% 
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh38_ID83_Denovo' )

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp6$Signature_name)))


# combined all signature data


signature_refsets_tmp <- bind_rows(
  signature_refsets_tmp1,
  signature_refsets_tmp2,
  signature_refsets_tmp3,
  signature_refsets_tmp4,
  signature_refsets_tmp5,
  signature_refsets_tmp6
) %>% 
  arrange(Source,Profile,Signature_set_name,Dataset,Signature_name,MutationType)



# Collect mutational profile matrix data ----------------------------------

exposure_refdata_tmp %>%
  select(Study,Cancer_Type,Sample,Dataset,Signature_set_name,Signature_name,Exposure)

seqmatrix_refdata_tmp <- exposure_refdata_tmp %>%
  left_join(signature_refsets_tmp %>% select(Profile,Signature_set_name,Signature_name,MutationType,Contribution)) %>%
  group_by(Study,Cancer_Type,Sample,Dataset,Profile,MutationType,Signature_set_name) %>%
  summarise(Mutations = sum(Exposure*Contribution)) %>%
  ungroup() %>%
  filter(str_starts(Signature_set_name,"Denovo_")) %>%  # only select the estimation based on denovo result
  select(-Signature_set_name) %>%
  arrange(Study,Cancer_Type,Sample,Dataset,Profile,MutationType)




# Save major R object -----------------------------------------------------------

seqmatrix_refdata <- seqmatrix_refdata_tmp
save(seqmatrix_refdata, file='UK100k-CRC-N2023_WGS_seqmatrix_refdata.RData') # study name, dataset and seqmatrix_refdata.Rdata

exposure_refdata <- exposure_refdata_tmp
save(exposure_refdata, file='UK100k-CRC-N2023_WGS_exposure_refdata.RData') # study name, dataset and exposure_refdata.Rdata

signature_refsets <-  signature_refsets_tmp
save(exposure_refdata, file='UK100k-CRC-N2023_WGS_signature_refsets.RData') # study name, dataset and siganture_refsets.RData



# NEXT, Collect sample level study data -----------------------------------





