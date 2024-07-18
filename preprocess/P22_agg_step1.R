rm(list=ls())
getwd()
setwd("D:/ych小组会文件/data/P22mousebrain")
library(devtools)
library(Seurat)
library(Signac)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
set.seed(1234)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
#library(DiffBind)

P22_mousebrain <- readRDS("D:/ych小组会文件/data/P22mousebrain/P22mousebrain_spatial_RNA_ATAC.rds")
image <- P22_mousebrain@images
#save(image, file = "D:/ych小组会文件/studyingR/P22_project/P22_image.RData")

atac_counts <- P22_mousebrain@assays[["peaks"]]@counts
metas <- P22_mousebrain@meta.data[["ATAC_clusters"]]
metas <- as.data.frame(metas)
rownames(metas) <- colnames(atac_counts)

combined <- CreateSeuratObject(counts = atac_counts, meta.data = metas, assay = "ATAC") 
#combined[["ATAC"]] <- P22_mousebrain@assays[["peaks"]]
combined@images <- P22_mousebrain@images

# ATAC analysis
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = 'lsi', dims = 2:20)
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:20)
combined <- FindClusters(object = combined, verbose = FALSE, algorithm = 3, resolution = 0.6)
DimPlot(object = combined, label = TRUE) + ggtitle("ATAC UMAP")
SpatialDimPlot(combined, label = TRUE, label.size = 3) + ggtitle("Spatial ATAC")

saveRDS(combined, file = "D:/ych小组会文件/studyingR/P22_project/P22_atac_seurat.rds")


source("D:/ych小组会文件/studyingR/P22_project/codes/Aggregate_data.R")
source("D:/ych小组会文件/studyingR/P22_project/codes/generate_aggregated_data.R")
source("D:/ych小组会文件/studyingR/P22_project/codes/estimateSizeFactorsForMatrix.R")

snn_matrix <- read.table(file = "Z:/yangchenghui/my_ych_project_P22brain/output/P22atac_snn.txt",sep = ",",header = TRUE, row.names = 1)

start_time <- Sys.time()
new_data <- Aggregate_data(combined, k_neigh = 10, snn_matrix = snn_matrix) ### 1.8min
# paired data need to aggregate on preprocess raw data
end_time <- Sys.time()
end_time - start_time

rownames(new_data$atac) <- rownames(atac_counts)
colnames(new_data$atac) <- colnames(atac_counts)
save(new_data, file = "D:/ych小组会文件/studyingR/P22_project/P22_agg_atac_k10.RData")

#################################################### aggregate RNA

rna_counts <- P22_mousebrain@assays[["SCT"]]@counts
metas <- P22_mousebrain@meta.data[["RNA_clusters"]]
metas <- as.data.frame(metas)
rownames(metas) <- colnames(rna_counts)

combined <- CreateSeuratObject(counts = rna_counts, meta.data = metas, assay = "RNA") 
combined@images <- P22_mousebrain@images

combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(object = combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
combined <- RunUMAP(combined, dims = 2:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, verbose = FALSE, algorithm = 3, resolution = 0.7)
DimPlot(object = combined, label = TRUE) + ggtitle("RNA UMAP")
SpatialDimPlot(combined, label = TRUE, label.size = 3) + ggtitle("Spatial RNA")

saveRDS(combined, file = "D:/ych小组会文件/studyingR/P22_project/P22_rna_seurat.rds")


source("D:/ych小组会文件/studyingR/P22_project/codes/Aggregate_data.R")
source("D:/ych小组会文件/studyingR/P22_project/codes/generate_aggregated_data.R")
source("D:/ych小组会文件/studyingR/P22_project/codes/estimateSizeFactorsForMatrix.R")

snn_matrix <- read.table(file = "Z:/yangchenghui/my_ych_project_P22brain/output/P22rna_snn.txt",sep = ",",header = TRUE, row.names = 1)

start_time <- Sys.time()
new_data <- Aggregate_data(combined, k_neigh = 10, snn_matrix = snn_matrix) ### 27.34771 secs
# paired data need to aggregate on preprocess raw data
end_time <- Sys.time()
end_time - start_time

rownames(new_data$rna) <- rownames(rna_counts)
colnames(new_data$rna) <- colnames(rna_counts)
save(new_data, file = "D:/ych小组会文件/studyingR/P22_project/P22_agg_rna_k10.RData")


################################ Agg_ATAC 可视化效果
load("D:/ych小组会文件/studyingR/P22_project/P22_agg_atac_k10.RData")
atac <- new_data$atac
agg_sample <- new_data[["cell_sample"]]
write.table(agg_sample, file = "Z:/yangchenghui/my_ych_project_P22brain_new/agg_P22_data/P22_agg_way.txt", sep = "\t")

combined <- CreateSeuratObject(counts = atac, assay = "ATAC") 
#combined[["ATAC"]] <- P22_mousebrain@assays[["peaks"]]
combined@images <- P22_mousebrain@images
