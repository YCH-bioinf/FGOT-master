# aggregated pbmc rna
rm(list=ls())
setwd("Z:/yangchenghui/BMMC_Zhang")
library(Seurat)
library(ggplot2)
library(patchwork)

rna_counts <- Read10X(data.dir = './rna/')
label <- read.csv2('RNA_celltype.txt')
label <- label$celltype
metas <- data.frame(label = label)
rownames(metas) <- colnames(rna_counts)

# Create Seurat object
combined <- CreateSeuratObject(counts = rna_counts, meta.data = metas,min.cells = 3, min.features = 200) 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 &nCount_RNA < 20000 & percent.mt < 10)

combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
#all.genes <- rownames(combined)
#combined <- ScaleData(combined, features = all.genes)
combined <- ScaleData(combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 1)
combined <- RunUMAP(combined, n.neighbors = 30, dims = 1:30)
DimPlot(combined, reduction = "umap")

DimPlot(combined, reduction = "umap", group.by = 'label',cols = c('B'='#F68282','CD4'='#31C53F','CD8'='#faf4cf','CLP'='#B95FBB','CMP'='#D4D915',
                                                                  'Ery'='#4B4BF7','GMP'='#ff9a36','HSC'='#2FF18B','LMPP'='#aeadb3','MEP'='#E6C122',
                                                                  'Mono'='#CCB1F1','MPP'='#25aff5','NK'='#A4DFF2'))
DimPlot(combined, reduction = "umap", group.by = "label",
        cols = c('B'='#F68282','CD4'='#31C53F','Ery'='#4B4BF7','CD8'='#faf4cf','NK'='#A4DFF2','CLP'='#B95FBB','Mono'='#CCB1F1','CMP'='#2FF18B',
                 'GMP'='#2FF18B','HSC'='#2FF18B','LMPP'='#2FF18B','MEP'='#2FF18B',
                 'MPP'='#2FF18B'))
saveRDS(combined,file = "BMMC_RNA_Seurat.rds")
FeaturePlot(combined,features = c('nFeature_RNA','nCount_RNA','percent.mt'))


Idents(combined) <- combined$label
fc <- FoldChange(combined, ident.1 = "CD8", ident.2 = "CD4")
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
head(fc)
FeaturePlot(combined,features = c('CAVIN2','GYPA','IFIT1B','SNCA','TMPRSS9'))

table(combined$label)
# remove clusters with its number smaller than 20
Idents(combined) <- combined$label
combined_new <- subset(combined, idents = c("B","CD4","CD8","CLP","Ery","Mono","NK"))

combined_new <- NormalizeData(object = combined_new, normalization.method = "LogNormalize", scale.factor = 10000)
combined_new <- FindVariableFeatures(combined_new, selection.method = "vst", nfeatures = 3000)
#all.genes <- rownames(combined)
#combined <- ScaleData(combined, features = all.genes)
combined_new <- ScaleData(combined_new)
combined_new <- RunPCA(combined_new, features = VariableFeatures(object = combined_new))
combined_new <- FindNeighbors(combined_new, dims = 1:30)
combined_new <- FindClusters(combined_new, resolution = 1)
combined_new <- RunUMAP(combined_new, n.neighbors = 30, dims = 1:30)

DimPlot(combined_new, reduction = "umap", group.by = 'label',cols = c('B'='#F68282','CD4'='#31C53F','CD8'='#faf4cf','CLP'='#B95FBB',
                                                                  'Ery'='#4B4BF7',
                                                                  'Mono'='#CCB1F1','NK'='#A4DFF2'))

saveRDS(combined_new,file = "BMMC_RNA_Seurat_after_filter.rds")
rna_data <- combined_new@assays$RNA@layers$data

dim(rna_data)
hvg <- VariableFeatures(object = combined_new)
#id <- match(hvg, rownames(combined_new@assays$RNA$counts))
#identical(rownames(combined_new@assays$RNA$counts)[id], hvg)
#rna_data_s <- rna_data[id,]
rownames(rna_data) <- rownames(combined_new@assays$RNA$counts)
colnames(rna_data) <- colnames(combined_new@assays$RNA$counts)

write.table(rna_data, file = "BMMC_RNA_data.txt",quote=FALSE,sep="\t")
write.table(hvg, file = "BMMC_RNA_hvg.txt",quote=FALSE,sep="\t")

metas_sub <- combined_new$label
write.table(metas_sub, file = "BMMC_RNA_label.txt",quote=FALSE,sep="\t")

snn_smooth <-  combined_new[['RNA_snn']] # used for compute cost
write.table(snn_smooth, file = "BMMC_rna_snn.txt", quote = FALSE,sep = "\t")
