#Load the required libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)

getwd()
setwd("Z:/yangchenghui/my_ych_project_simulation1/simulated_data1/")

######################################################################
# Load data
batch1_counts <- read.table(file = "simulation1_batch1.txt")
batch2_counts <- read.table(file = "simulation1_batch2.txt")

label1 <- read.table(file = "simulation1_label1.txt")
label2 <- read.table(file = "simulation1_label2.txt")
label1 <- label1$x
label2 <- label2$x

######################################################################
# batch1 analysis
# Create Seurat object
batch1 <- CreateSeuratObject(counts = batch1_counts)

# RNA analysis
batch1 <- NormalizeData(batch1)
x1 <- GetAssayData(batch1, slot = "data")
#x1 = batch1@assays$RNA@data
batch1 <- FindVariableFeatures(batch1, selection.method = "vst", nfeatures = 2000)
all.genes_rna <- rownames(batch1)
batch1 <- ScaleData(batch1, features = all.genes_rna)
batch1 <-RunPCA(batch1, features = VariableFeatures(object = batch1))
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# ElbowPlot(batch1)
batch1 <- FindNeighbors(batch1, dims = 1:15)
batch1 <- FindClusters(batch1, resolution = 0.5)
batch1 <- RunUMAP(batch1, dims = 1:15, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

Idents(batch1) <- label1
DimPlot(batch1, reduction = "umap.rna",label = TRUE)

# markers
markers <- FindAllMarkers(batch1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers[c('cluster','gene')], file = "Z:/yangchenghui/my_ych_project_simulation1/feature_selected_simulated1_data/Batch1_marker_genes.txt", sep = "\t")

batch1_top50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# pdf('E17.5_group1_markers_heatmap.pdf')
DoHeatmap(batch1, features = batch1_top50$gene)


######################################################################
# batch2 analysis
# Create Seurat object
batch2 <- CreateSeuratObject(counts = batch2_counts)

# RNA analysis
batch2 <- NormalizeData(batch2)

x2 <- GetAssayData(batch2, slot = "data")
#x2 = batch2@assays$RNA@data
batch2 <- FindVariableFeatures(batch2, selection.method = "vst", nfeatures = 2000)
all.genes_rna <- rownames(batch2)
batch2 <- ScaleData(batch2, features = all.genes_rna)
batch2 <-RunPCA(batch2, features = VariableFeatures(object = batch2))
# print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# ElbowPlot(batch1)
batch2 <- FindNeighbors(batch2, dims = 1:15)
batch2 <- FindClusters(batch2, resolution = 0.5)
batch2 <- RunUMAP(batch2, dims = 1:15, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

Idents(batch2) <- label2
DimPlot(batch2, reduction = "umap.rna",label = TRUE)

# markers
markers <- FindAllMarkers(batch2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers[c('cluster','gene')], file = "Z:/yangchenghui/my_ych_project_simulation1/feature_selected_simulated1_data/Batch2_marker_genes.txt", sep = "\t")

batch2_top50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
# pdf('E17.5_group1_markers_heatmap.pdf')
DoHeatmap(batch2, features = batch2_top50$gene)


fea = intersect(batch1_top50$gene, batch2_top50$gene)

###########################################################

x11 = x1[fea,]
write.table(x11, file = "Z:/yangchenghui/my_ych_project_simulation1/feature_selected_simulated1_data/Batch1_feature50_selected_data.txt", sep = "\t")

df <- batch1_top50[c('cluster','gene')]
selected_rows <- df %>% filter(gene %in% fea)
write.table(selected_rows, file = "Z:/yangchenghui/my_ych_project_simulation1/feature_selected_simulated1_data/Batch1_top50_marker_genes.txt", sep = "\t")

x22 = x2[fea,]
write.table(x22, file = "Z:/yangchenghui/my_ych_project_simulation1/feature_selected_simulated1_data/Batch2_feature50_selected_data.txt", sep = "\t")

