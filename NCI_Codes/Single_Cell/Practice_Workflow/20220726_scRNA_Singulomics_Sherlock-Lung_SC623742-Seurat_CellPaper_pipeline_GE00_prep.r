library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(clustree)
library(gridExtra)
library(ggridges)
library(tidyverse)

rm(list=ls())
setwd("/data/Sherlock_Lung/zhaow5/scRNA_test/SC623742/")



piano.data <- Read10X("./filtered_feature_bc_matrix")
piano.raw <- CreateSeuratObject(piano.data,project="SC623742/")
aggr.file <- read_csv("Singulomics_Lung_SC623742_aggregation.csv")
sample.idx <- gsub("\\S+-([1|2])","\\1",rownames(piano.raw@meta.data))
orig_id <- as.character(1:nrow(aggr.file))
new_id <- aggr.file$sample_id
piano.raw@meta.data[,'sample'] <- as.character(plyr::mapvalues(sample.idx, from=orig_id,to=new_id))
table(piano.raw@meta.data$sample)


----------------------------------------------------------------------------------------------------------



piano.raw[["percent.mt"]] <- PercentageFeatureSet(piano.raw, pattern="^MT-")

ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = piano.raw@assays$RNA@data), value = TRUE)
percent.ribo <- Matrix::colSums(piano.raw@assays$RNA@counts[ribo.genes, ])/Matrix::colSums(piano.raw@assays$RNA@data)
piano.raw <- AddMetaData(object=piano.raw, metadata=percent.ribo, col.name = 'percent.ribo')
piano.raw$log10GenesPerUMI <- log10(piano.raw$nFeature_RNA) / log10(piano.raw$nCount_RNA)

setwd("/data/Sherlock_Lung/zhaow5/scRNA_test/SC623742//General")


#pdf("QC_violin_plot.pdf",height=6,width=15)
#VlnPlot(piano.raw,features=c("nCount_RNA","nFeature_RNA","percent.mt"))
#dev.off()

-----------------------------------------------------------------------------------------------------------------


piano <- subset(piano.raw,subset=nFeature_RNA>500 & nCount_RNA>1000 & percent.mt<5 & log10GenesPerUMI>0.8)

table(piano.raw$sample)
table(piano$sample)
rawcnt <- table(piano.raw$sample)
filter.cnt <- table(piano$sample)
out <- cbind(rawcnt,filter.cnt)
write.table(out,"cell_counts_before_after_filter.txt",sep="\t")


-----------------------------------------------------------------------------------------------------------------


piano <- NormalizeData(piano,normalization.method="LogNormalize")

load("/data/zhaow5/GENOMEDIR/scRNAseq/cycle.rda")

piano <- CellCycleScoring(piano, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

piano <- ScaleData(piano,features=rownames(piano))
piano <- FindVariableFeatures(piano,dispersion.cutoff=c(0.5,Inf),mean.cutoff=c(0.1,Inf) )
piano <- RunPCA(piano)

pdf(paste("PCA_cell_cycle.pdf",sep=""),width=15,height=6)
DimPlot(piano,reduction = "pca",group.by= "Phase",split.by = "Phase")
dev.off()

------------------------------------------------------------------------------------------------------------------

print(piano[["pca"]],dims=1:5,nfeatures=10)
pdf("GE01_PCA_analysis.pdf",width=10,height=10)
VizDimLoadings(piano, dims=1:4,reduction="pca")
DimPlot(piano,reduction="pca",group.by="sample")
DimHeatmap(piano, dims=1:9,cells=500, balanced=TRUE)
ElbowPlot(piano,reduction="pca",ndims=30)
dev.off()

my.cols <- rev(brewer.pal(11, "RdYlBu"))
tissue.type.genes <- c("EPCAM","CLDN5","COL1A2","PTPRC")
pdf(paste("GE01_PCA_tissue_genes.pdf",sep=""),width=10,height=10)
FeaturePlot(piano, reduction="pca",features=tissue.type.genes,cols=my.cols,ncol=2,min.cutoff="q10",max.cutoff="q90",pt.size=0.2)
dev.off()

---------------------------------------------------------------------------------------------------------------------

### first round of clustering (group cells by epithelial/immune/stromal markers)
n.pcs <- 25
#res.used <- 1.8
seed.use=123
res.used <- seq(0.1,1.5,by=0.2)

for(i in res.used){
  piano <- FindNeighbors(piano,dims=1:n.pcs)
  piano <- FindClusters(piano,resolution=i)
}

clus.tree.out <- clustree(piano) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
pdf(paste("GE01_General_cluster_comparison_",n.pcs,"pcs_seed",seed.use,".pdf",sep=""),width=14,height=12)
print(clus.tree.out)
dev.off()


n.pcs <- 25
#res.used <- 1.8
seed.use=123
res.used=0.5

piano <- FindNeighbors(piano,dims=1:n.pcs)
piano <- FindClusters(piano,resolution=res.used)
piano <- RunUMAP(piano, dims=1:n.pcs, seed.use=seed.use)
pdf(paste("GE01_General_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_UMAP_graphic_cluster.pdf",sep=""))
DimPlot(piano, reduction="umap",label=TRUE,pt.size=0.1)
dev.off()

pdf(paste("GE01_General_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_UMAP_graphic_sample.pdf",sep=""))
DimPlot(piano, reduction="umap",label=TRUE,pt.size=0.1,group.by="sample")
dev.off()

piano <- RunTSNE(piano, dims=1:n.pcs, seed.use=seed.use)
pdf(paste("GE01_General_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_TSNE_graphic_cluster.pdf",sep=""))
DimPlot(piano, reduction="tsne",label=TRUE,pt.size=0.1)
dev.off()

pdf(paste("GE01_General_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_TSNE_graphic_sample.pdf",sep=""))
DimPlot(piano, reduction="tsne",label=TRUE,pt.size=0.1,group.by="sample")
dev.off()

tissue.type.genes <- c("EPCAM","CLDN5","COL1A2","PTPRC")
pdf(paste("GE01_tSNE_tissue_genes.pdf",sep=""),width=10,height=10)
FeaturePlot(piano, reduction="tsne",features=tissue.type.genes,cols=my.cols,ncol=2,min.cutoff="q10",max.cutoff="q90",pt.size=0.2)
dev.off()

---------------------------------------------------------------------------------------------------------------------------------------

#### integrate samples
split_seurat <- SplitObject(piano, split.by = "sample")
#split_seurat <- split_seurat[c("ctrl", "stim")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]],variable.features.n = 20000, vars.to.regress = c("percent.mt"))
}



# Select the most variable features to use for integration


nfeatures=15000
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = nfeatures) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

saveRDS(seurat_integrated, paste("GE00_integrated_seurat_n",nfeatures,".rds",sep=""))
saveRDS(piano, "GE00_aggr_seurat.rds")

