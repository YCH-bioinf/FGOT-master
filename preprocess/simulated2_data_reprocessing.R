#Load the required libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

getwd()
setwd(".")

######################################################################
# Load data
batch1_counts <- read.table(file = "simulated_data/simulation2_batch1.txt")
batch2_counts <- read.table(file = "simulated_data/simulation2_batch2.txt")

label1 <- read.table(file = './simulated_data/simulation2_label1.txt')
label2 <- read.table(file = './simulated_data/simulation2_label2.txt')

######################################################################
# batch1 analysis
# Create Seurat object
batch1 <- CreateSeuratObject(counts = batch1_counts)

# RNA analysis
batch1 <- NormalizeData(batch1)
x1 = batch1@assays$RNA@data
batch1 <- FindVariableFeatures(batch1, selection.method = "vst", nfeatures = 2000)
all.genes_rna <- rownames(batch1)
batch1 <- ScaleData(batch1, features = all.genes_rna)
batch1 <-RunPCA(batch1, features = VariableFeatures(object = batch1))
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# ElbowPlot(batch1)
batch1 <- FindNeighbors(batch1, dims = 1:15)
batch1 <- FindClusters(batch1, resolution = 0.6)
batch1 <- RunUMAP(batch1, dims = 1:15, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

Idents(batch1) <- "seurat_clusters"
DimPlot(batch1, reduction = "umap.rna",label = TRUE) # According to seurat_clusters
Idents(batch1) <- label1
DimPlot(batch1, reduction = "umap.rna",label = TRUE) # According to true label

# markers
markers <- FindAllMarkers(batch1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers[c('cluster','gene')], file = "./simulated_feature_selected_data/Batch1_marker_genes.txt", sep = "\t")

batch1_top50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# pdf('E17.5_group1_markers_heatmap.pdf')
DoHeatmap(batch1, features = batch1_top50$gene)


######################################################################
# batch2 analysis
# Create Seurat object
batch2 <- CreateSeuratObject(counts = batch2_counts)

# RNA analysis
batch2 <- NormalizeData(batch2)
x2 = batch2@assays$RNA@data
batch2 <- FindVariableFeatures(batch2, selection.method = "vst", nfeatures = 2000)
all.genes_rna <- rownames(batch2)
batch2 <- ScaleData(batch2, features = all.genes_rna)
batch2 <-RunPCA(batch2, features = VariableFeatures(object = batch2))
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# ElbowPlot(batch1)
batch2 <- FindNeighbors(batch2, dims = 1:15)
batch2 <- FindClusters(batch2, resolution = 0.1)
batch2 <- RunUMAP(batch2, dims = 1:15, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

Idents(batch2) <- label2
DimPlot(batch2, reduction = "umap.rna",label = TRUE)

# markers
markers <- FindAllMarkers(batch2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers[c('cluster','gene')], file = "./simulated_feature_selected_data/Batch2_marker_genes.txt", sep = "\t")

batch2_top50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# pdf('E17.5_group1_markers_heatmap.pdf')
DoHeatmap(batch2, features = batch2_top50$gene)

batch1_top50_gene = batch1_top50$gene
batch2_top50_gene = batch2_top50$gene
fea50 = intersect(batch1_top50_gene, batch2_top50_gene)

###########################################################

x1_50 = x1[fea50,]
write.table(x1_50, file = "./simulated_feature_selected_data/Batch1_feature50_selected_data.txt", sep = "\t")

write.table(batch1_top50[c('cluster','gene')], file = "./simulated_feature_selected_data/Batch1_top50_marker_genes.txt", sep = "\t")


x2_50 = x2[fea50,]
write.table(x2_50, file = "./simulated_feature_selected_data/Batch2_feature50_selected_data.txt", sep = "\t")

write.table(batch2_top50[c('cluster','gene')], file = "./simulated_feature_selected_data/Batch2_top50_marker_genes.txt", sep = "\t")


