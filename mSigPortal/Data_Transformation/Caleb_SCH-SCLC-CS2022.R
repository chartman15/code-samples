
# Paper: Molecular subtyping of small-cell lung cancer based on mutational signatures with different genomic features and therapeutic strategies (PMID 36178064)

# Data: Supp Table 2. Proportions of somatic mutational signatures and the cluster subtype of each patient.

library(tidyverse)
library(readxl)

setwd("C:/Users/hartmancas/Desktop/mSigPortal/mSigPortal_Data_Cleaning/SCH-SCLC-CS2022")

# load collected mutational signatures, including different COSMIC signatures
load('signature_refsets.RData')


### mutational signature data----------------------------------------------------------------------------------------------------------------------------------------------- 

# decomposed signature SBS96, extract BOTH decomposed COSMIC signatures and new decomposed signatures 

