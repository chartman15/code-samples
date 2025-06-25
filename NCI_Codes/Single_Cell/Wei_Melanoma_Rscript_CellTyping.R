###Calculate and plot  expression of cell markers
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(clustree)
library(gridExtra)
library(ggridges)
library(DoubletFinder)
library(tidyverse)
require(ggridges)
require(ggplot2)

n.pcs <- 30
res.used=0.9
seed.use=123

########################################################################################
cell.genes <- read_tsv("/data/zhaow5/Melanoma/Genelist/Melanoma_cell_type_markers.txt")

genes_to_check <- cell.genes %>% filter(Type=="Primary") %>% pull(gene)
#genes_to_check <- cell.genes %>% select(gene)
genes_to_check <-  as.character(unique(genes_to_check))
genes_to_check <-genes_to_check[genes_to_check %in%rownames(tiss.merge@assays$SCT)]

########################################################################################


pdf(paste("GE01_General_SCT_Harmony_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_primary_marker_DotPlot.pdf", sep=""),  width = 12, height = 20)
DotPlot(tiss.merge, features=genes_to_check) + coord_flip()
dev.off()


#######################################################################################
setwd("/data/Sherlock_Lung/zhaow5/scRNA_test/melanoma/Patient_2253/01_General")
load("Patient_2253_general_cluster.RData")



#cell.genes <- cell.genes[genes_to_check,]

#umap.coor <- Embeddings(tiss.merge,reduction = "umap")


cell.genes <- read_tsv("/data/zhaow5/Melanoma/Genelist/Melanoma_cell_type_markers.txt")

genes_to_check <- cell.genes %>% filter(Type=="Primary")  %>% filter(gene%in%rownames(tiss.merge@assays$SCT))
#genes_to_check <- cell.genes %>% select(gene)

data_temp <- as.data.frame(as.matrix(GetAssayData(object=tiss.merge[["SCT"]], slot="data")))
tsne.coor <- Embeddings(tiss.merge,reduction="tsne")
umap.coor <- Embeddings(tiss.merge,reduction="umap")

ggplot.list.umap <- list()
ggplot.list.tsne <- list()

ggplot.list.2 <- list()
for(i in 1:nrow(genes_to_check)){
  gene.name <- genes_to_check$gene[i]
  cell.type <- genes_to_check$cell[i]
  gene.exp <- data_temp[gene.name,]
  clusters <- tiss.merge@meta.data$seurat_clusters
  
  df.tsne <- as.data.frame(cbind(tsne.coor, as.data.frame(t(gene.exp)),clusters))
  colnames(df.tsne) <- c(colnames(tsne.coor),"gene.exp","clusters")
  df.umap <- as.data.frame(cbind(umap.coor, as.data.frame(t(gene.exp)),clusters))
  colnames(df.umap) <- c(colnames(umap.coor),"gene.exp","clusters")
  
  ggplot.list.umap[[i]] <- ggplot(df.umap, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour=gene.exp), size=0.7) +
    scale_colour_gradientn(colours=c("lightgray", "blue"),limits=quantile(df.umap$gene.exp,probs=c(0.05,0.95)),oob=scales::squish) +
    labs(title=paste(gene.name,"(",cell.type,")",sep=""))+
    theme(plot.title=element_text(hjust=0.5))
  
  ggplot.list.tsne[[i]] <- ggplot(df.tsne, aes(tSNE_1, tSNE_2)) +
    #  ggplot.list[[i]] <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour=gene.exp), size=0.7) +
    scale_colour_gradientn(colours=c("lightgray", "blue"),limits=quantile(df.tsne$gene.exp,probs=c(0.05,0.95)),oob=scales::squish) +
    labs(title=paste(gene.name,"(",cell.type,")",sep=""))+
    theme(plot.title=element_text(hjust=0.5))
  
  ggplot.list.2[[i]] <-ggplot(df.umap) +
    geom_boxplot(aes(x=clusters,y=gene.exp))  + 
    ggtitle(paste(gene.name,"(",cell.type,")",sep="")) + ylab("gene expression (log)")+
    theme(plot.title=element_text(hjust=0.5))
}

n<- length(ggplot.list.tsne)
ncol <- floor(sqrt(n))
pdf(paste("GE01_General_SCT_Harmony_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_Boxplot_primary_markers.pdf", sep=""), width = 20, height = 20)
do.call("grid.arrange", c(ggplot.list.2, ncol=ncol))
dev.off()

#pdf(paste("GE01_General_SCT_Harmony_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_TSNE_4Sets_markers.pdf", sep=""),width=30,height=30)
#do.call("grid.arrange", c(ggplot.list.tsne, ncol=ncol))
#dev.off()

#pdf(paste("GE01_General_SCT_Harmony_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_UMAP_4Sets_markers.pdf", sep=""),width=20,height=20)
#do.call("grid.arrange", c(ggplot.list.umap, ncol=ncol))
#dev.off()


pdf(paste("GE01_General_SCT_Harmony_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_UMAP_primary_markers.pdf", sep=""),width=30,height=30)
FeaturePlot(tiss.merge, reduction="umap",features=genes_to_check$gene,ncol=floor(sqrt(nrow(genes_to_check))),min.cutoff="q10",max.cutoff="q90",pt.size=0.1)
dev.off()

pdf(paste("GE01_General_SCT_Harmony_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_TSNE_primary_markers.pdf", sep=""),width=30,height=30)
FeaturePlot(tiss.merge, reduction="tsne",features=genes_to_check$gene,ncol=floor(sqrt(nrow(genes_to_check))),min.cutoff="q10",max.cutoff="q90",pt.size=0.1)
dev.off()





###Calculate and plot average expression of broad cell type specific gene sets 

require(ggridges)
require(ggplot2)
# # # 
cell.genes <- read.table("/data/zhaow5/Melanoma/Genelist/Melanoma_cell_type_markers.txt", header=T,sep="\t")
cell.genes <- cell.genes[as.character(cell.genes$Type)=="Primary",]
cell.types  <- as.character(unique(cell.genes$cell))
tsne.coor <- Embeddings(tiss.merge,reduction = "tsne")
## 
ggplot.list <- list()
ggplot.list.2 <- list()
# 
# rm(temp)
data_temp <- as.data.frame(as.matrix(GetAssayData(object = tiss.merge, slot = "data")))
for(i in 1:length(unique(cell.types))){
  genes <- as.character(cell.genes$gene[which(cell.genes$cell==cell.types[i])])
  genes <- intersect(genes,rownames(data_temp))
  gene.exp <- colMeans(as.matrix(data_temp[genes,]))[row.names(tsne.coor)]
  #  gene.exp <- colSums(as.matrix(data_temp[genes,]))[row.names(tsne.coor)]
  clusters <- tiss@meta.data$RNA_snn_res.1.1
  # Make ggplot friendly 
  temp <- as.data.frame(cbind(tsne.coor, as.data.frame(gene.exp), as.data.frame(clusters)))
  # Plot with ggplot 
  # ggplot.list[[i]] <- ggplot(temp, aes(UMAP_1, UMAP_2)) + 
  ggplot.list[[i]]<- ggplot(temp, aes(tSNE_1, tSNE_2)) + 
    geom_point(aes(colour = gene.exp),size=1) + 
    scale_colour_gradient(low = "grey95", high = "red") + 
    labs(title = cell.types[i], subtitle = paste(genes, collapse = ", "))
  # Boxplot per cluster 
  ggplot.list.2[[i]] <- ggplot(temp, aes(x = clusters, y = gene.exp)) + 
    geom_boxplot() + 
    ggtitle(cell.types[i]) + ylab("Average gene expression (log)")
}
# Plot all 
require(grid)
require(gridExtra)
require(gbm)
n <- length(ggplot.list)
nCol <- floor(sqrt(n))
# Expression on UMAP


pdf(paste("General_cluster_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_TSNE_with_avg_expression_of_cell_markers.pdf", sep=""),20,20)
do.call("grid.arrange", c(ggplot.list, ncol=nCol))
dev.off()

# Expression per cluster boxplots 

pdf(paste("General_cluster_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_boxplots_with_avg_expression_of_cell_markers.pdf", sep=""),20,20)
do.call("grid.arrange", c(ggplot.list.2, ncol=nCol))
dev.off()


#pdf(paste("Nonmixed_cluster_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_UMAP_with_average_expression_of_cell_markers.pdf", sep=""),15,15)
#do.call("grid.arrange", c(ggplot.list, ncol=nCol))
#dev.off()
# Expression per cluster boxplots 
#pdf(paste("Nonmixed_cluster_",n.pcs,"pcs_res",res.used,"_seed",seed.use,"_boxplots_with_average_expression_of_cell_markers.pdf", sep=""),15,15)
#do.call("grid.arrange", c(ggplot.list.2, ncol=nCol))
#dev.off()
#########################
###  general cluster annotation
tiss.merge@meta.data$GE01_seurat_cluster_id <- tiss.merge@meta.data$seurat_clusters
tiss.merge@meta.data$GE01_main_cluster <- tiss.merge@meta.data$seurat_clusters
cluster.id <- sort(as.numeric(unique(as.character(tiss.merge@meta.data$GE01_seurat_cluster_id))))
free_annotation <- c("fibroblasts, smooth muscle", #0
                     "melanocytes", #1
                     "keratinocytes", #2
                     "endothelial", #3
                     "mixed", #4
                     "mixed", #5
                     "melanocytes", #6
                     "fibroblasts, smooth muscle", #7
                     "mixed", #8
                     "keratinocytes", #9
                     "fibroblasts, smooth muscle", #10
                     "melanocytes", #11
                     "keratinocytes", #12
                     "melanocytes", #13
                     "mixed", #14
                     "mixed", #15
                     "mixed", #16
                     "keratinocytes", #17
                     "endothelial", #18
                     "fibroblasts, smooth muscle", #19
                     "keratinocytes", #20
                     "mixed", #21
                     "melanocytes", #22
                     "melanocytes", #23
                     "fibroblasts, smooth muscle",#24
                     "melanocytes", #25
                     "melanocytes", #26
                     "mixed", #27
                     "melanocytes", #28
                     "neuron cells", #29
                     "keratinocytes" #30
) 

tiss.merge@meta.data$GE01_main_cluster <- as.character(plyr::mapvalues(tiss.merge@meta.data$seurat_clusters, from=cluster.id, to=free_annotation))

table(tiss.merge@meta.data$GE01_main_cluster)
save(tiss.merge,file="Patient_2253_general_cluster.RData")

##############################
#