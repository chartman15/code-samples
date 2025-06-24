# Paper: Deep whole-genome analysis of 494 hepatocellular carcinomas (PMID 38355797)

# Data: Supp Table 3a-3c: Signature contribution. 
#       Extended Data Fig. 2: Profiles of all mutational signatures in CLCA.


library(tidyverse)
library(readxl)

setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/CLCA-HCC-N2024")

# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')

### Supp Table 3a-3c: mutational signature data----------------------------------------------------------------------------------------------------------------------------------------------- 

# decomposed SBS96
TableData <- read_xlsx(path = 'Supplementary_Table_3a.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp1 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name', values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='SBS96',Signature_set_name = "COSMIC_v3.2_Signatures_GRCh37_SBS96_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

sigtmp <- signature_refsets_tmp1 %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_SBS96_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp1 <- bind_rows(signature_refsets %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_SBS96') %>%
    filter(Signature_name %in% sigtmp), signature_refsets_tmp1) %>%
    mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh37_SBS96_Denovo')

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp1$Signature_name)))



# decomposed DBS78
TableData <- read_xlsx(path = 'Supplementary_Table_3b.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp2 <- TableData %>% 
  pivot_longer(cols = -MutationType, names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='DBS78',Signature_set_name = "COSMIC_v3.2_Signatures_GRCh37_DBS78_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

sigtmp <- signature_refsets_tmp2 %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_DBS78_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp2 <- bind_rows(signature_refsets %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_DBS78') %>%
  filter(Signature_name %in% sigtmp), signature_refsets_tmp2) %>%
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh37_DBS78_Denovo')

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp2$Signature_name)))



# decomposed ID83
TableData <- read_xlsx(path = 'Supplementary_Table_3c.xlsx',skip = 2,col_names = T,.name_repair = 'check_unique')
colnames(TableData)[1] <- "MutationType"

signature_refsets_tmp3 <- TableData %>% 
  pivot_longer(cols = -MutationType,names_to = 'Signature_name',values_to = 'Contribution') %>% 
  mutate(Source='Study_signatures',Profile='ID83',Signature_set_name = "COSMIC_v3.2_Signatures_GRCh37_ID83_Denovo", Dataset = "WGS", Strand_info = "N", Strand = NA_character_) %>% 
  select(Source, Profile, Signature_set_name, Dataset, Strand_info, Strand, Signature_name, MutationType, Contribution)

sigtmp <- signature_refsets_tmp3 %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_ID83_Denovo') %>% pull(Signature_name) %>% unique()

signature_refsets_tmp3 <- bind_rows(signature_refsets %>% filter(Signature_set_name == 'COSMIC_v3.2_Signatures_GRCh37_ID83') %>%
                                      filter(Signature_name %in% sigtmp), signature_refsets_tmp3) %>%
  mutate(Source = 'Study_signatures',Signature_set_name = 'COSMIC_v3.2_Signatures_GRCh37_ID83_Denovo')

#check all signatures included
unique(sort(sigtmp) == sort(unique(signature_refsets_tmp3$Signature_name)))



# combined all signature data

signature_refsets_tmp_final <- bind_rows(
  signature_refsets_tmp1,
  signature_refsets_tmp2,
  signature_refsets_tmp3
) %>% 
  arrange(Source,Profile,Signature_set_name,Dataset,Signature_name,MutationType)



# Save major R objects-----------------------------------------------------------

signatures_refsets <- signature_refsets_tmp_final
save(signatures_refsets, file='CLCA-HCC-N2024_WGS_signature_refsets.RData') # study name, dataset and signature_refsets.RData

