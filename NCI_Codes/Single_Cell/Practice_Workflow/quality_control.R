# 8/31/2023
# Caleb Hartman 
# HBC Single-cell RNA-seq workshop

#-------------------------------------------------------------------------------------
# ***Single-cell RNA-seq analysis - Quality control (QC) setup***

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

sessionInfo()

# Reading in a single sample

# How to read in 10X data for a single sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")

# Turn count matrix into a Seurat object (output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts, min.features = 100) 

# The 'min.features' argument specifies the minimum number of genes that need to be detected per cell. This argument will filter out poor quality cells that likely just have random barcodes encapsulated without any cell present.

# Explore the metadata
head(ctrl@meta.data)

# orig.ident: this often contains the sample identity if known, but will default to “SeuratProject”

# nCount_RNA: number of UMIs per cell

# nFeature_RNA: number of genes detected per cell

# Reading in multiple samples 

# In practice, you will likely have several samples that you will need to read in data for, and that can get tedious and error-prone if you do it one at a time. So, to make the data import into R more efficient we can use a for loop, which will iterate over a series of commands for each of the inputs given and create seurat objects for each of our samples.

# In R, the for loop has the following structure/syntax:
#
#   
#   for (variable in input)
#   {
#     command1
#     command2
#     command3
#   }

# Create a Seurat object for each sample
for (file in c("ctrl_raw_feature_bc_matrix", "stim_raw_feature_bc_matrix"))
{
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, min.features = 100, 
  project = file)
  assign(file, seurat_obj)
}

# Check the metadata in the new Seurat objects
head(ctrl_raw_feature_bc_matrix@meta.data)
head(stim_raw_feature_bc_matrix@meta.data)

# Create a merged Seurat object; do this so we can run QC for both sample groups together 
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, y = stim_raw_feature_bc_matrix, add.cell.id = c("ctrl", "stim"))

# Because the same cell IDs can be used for different samples, we add a sample-specific prefix to each of our cell IDs using the 'add.cell.id' argument.

# Seurat now has functionality to merge many samples together. You can do this quite easily by adding all sample objects to the y argument in a vector format. An example is provided below:

# merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, y = c(stim1_raw_feature_bc_matrix, stim2_raw_feature_bc_matrix, stim3_raw_feature_bc_matrix), add.cell.id = c("ctrl", "stim1", "stim2", "stim3"))

head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)

# ------------------------------------------------------------------------------------
# ***Single-cell RNA-seq analysis - Quality control (QC) analysis***

view(merged_seurat@meta.data)

# In order to create the appropriate plots for the quality control analysis, we need to calculate some additional metrics. These include:
  
# 1. number of genes detected per UMI: this metric gives us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data) (complexity/novelty score)

# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)


# 2. mitochondrial ratio: this metric will give us a percentage of cell reads originating from the mitochondrial genes

# For each cell, the function takes the sum of counts across all genes (features) belonging to the “Mt-“ set, and then divides by the count sum for all genes (features). This value is multiplied by 100 to obtain a percentage value. For our analysis, rather than using a percentage value we would prefer to work with the ratio value. As such, we will reverse that last step performed by the function by taking the output value and dividing by 100.

# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# NOTE: The pattern provided (“^MT-“) works for human gene names. You may need to adjust the pattern argument depending on your organism of interest. Additionally, if you weren’t using gene names as the gene ID then this function wouldn’t work as we have used it above as the pattern will not suffice. Since there are caveats to using this function, it is advisable to manually compute this metric. If you are interested, we have code available to compute this metric on your own.

# add cell IDs and condition info to metadata of merged_seurat object 

# Create metadata dataframe
metadata <- merged_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"

# Rename columns
metadata <- metadata %>%
dplyr::rename(seq_folder = orig.ident, nUMI = nCount_RNA, nGene = nFeature_RNA)

view(metadata)

# final metadata table will have rows that correspond to each cell, and columns with information about those cells

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 250)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score) (ratio of nGenes over nUMI)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs (upper right quadrant of the plot). Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot. With this plot we also evaluate the slope of the line, and any scatter of data points in the bottom right hand quadrant of the plot. These cells have a high number of UMIs but only a few number of genes. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).

# Mitochondrial read fractions are only high in particularly low count cells with few detected genes (darker colored data points). This could be indicative of damaged/dying cells whose cytoplasmic mRNA has leaked out through a broken membrane, and thus, only mRNA located in the mitochondria is still conserved. We can see from the plot, that these cells are filtered out by our count and gene number thresholds.

# Filtering 

# Cell filtering 

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, subset = (nUMI >= 500) & (nGene >= 250) &(log10GenesPerUMI > 0.80) & (mitoRatio < 0.20))

# Gene filtering 

# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene whether or not there are more than zero counts per cell
nonzero <- counts > 0

# Now, we will perform some filtering by prevalence. If a gene is only expressed in a handful of cells, it is not particularly meaningful as it still brings down the averages for all other cells it is not expressed in. For our data we choose to keep only genes which are expressed in 10 or more cells. By using this filter, genes which have zero counts in all cells will effectively be removed.

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

#-------------------------------------------------------------------------------------
# Exercises

# Extract the new metadata from the filtered Seurat object using the code provided below:
  
# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

# 1. Perform all of the same QC plots using the filtered data.

# Visualize the number of cell counts per sample
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells - clean metadata")

# Visualize the number UMIs/transcripts per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 250)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score) (ratio of nGenes over nUMI)
metadata_clean %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

# 2. Report the number of cells left for each sample, and comment on whether the number of cells removed is high or low. Can you give reasons why this number is still not ~12K (which is how many cells were loaded for the experiment)?

# Answer: a little under 15k cells for ctrl and stim. Number cells removed fairly low. In an ideal world, you would expect the number of unique cellular barcodes to correpsond to the number of cells you loaded. However, this is not the case as capture rates of cells are only a proportion of what is loaded. For example, the inDrops cell capture efficiency is higher (70-80%) compared to 10X (used here) which is between 50-60%.

# The cell numbers can also vary by protocol, producing cell numbers that are much higher than what we loaded. For example, during the inDrops protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. Similarly, with the 10X protocol there is a chance of obtaining only a barcoded bead in the emulsion droplet (GEM) and no actual cell. Both of these, in addition to the presence of dying cells can lead to a higher number of cellular barcodes than cells.

# 3. After filtering for nGene per cell, you should still observe a small shoulder to the right of the main peak. What might this shoulder represent?

# Answer: maybe a distinct cell population with more transcriptome diversity (more genes)
  
# 4. When plotting the nGene against nUMI do you observe any data points in the bottom right quadrant of the plot? What can you say about these cells that have been removed?

# Answer: The cells that were removed were those with high nUMI but low numbers of genes detected. These cells had many captured transcripts but represent only a small number of genes. These low complexity cells could represent a specific cell type (i.e. red blood cells which lack a typical transcriptome), or could be due to some other strange artifact or contamination.
#-------------------------------------------------------------------------------------

# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")

# END 