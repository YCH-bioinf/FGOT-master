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

load("D:/ych小组会文件/studyingR/P22_project/P22_image.RData")

load("D:/ych小组会文件/studyingR/P22_project/P22_agg_rna_k10.RData")
rna <- new_data$rna

load("D:/ych小组会文件/studyingR/P22_project/P22_agg_atac_k10.RData")
atac <- new_data$atac
peak <- rownames(atac)
peak1 <- strsplit(peak,"-")
peaks <- matrix(0, nrow = length(peak),ncol = 1)
for (i in 1:length(peaks)) {
  peaks[i] <- paste0(peak1[[i]][1],":",peak1[[i]][2],"-",peak1[[i]][3])
}
peaks <- as.vector(peaks)
rownames(atac) <- peaks

write.table(rna, file = "D:/ych小组会文件/studyingR/P22_project/P22_agg_rna.txt", quote = FALSE,sep = "\t")
write.table(atac, file = "D:/ych小组会文件/studyingR/P22_project/P22_agg_atac.txt", quote = FALSE,sep = "\t")

combined <- CreateSeuratObject(counts = rna)
chrom_assay <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = 'mm10',
  min.cells = 0
)
combined[["ATAC"]] <- chrom_assay
combined@images <- image

# RNA analysis
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(object = combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(object = combined)
combined <- RunPCA(combined, pc.genes = combined@var.genes, pcs.compute = 40, do.print = FALSE)
PCAPlot(object = combined)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, verbose = FALSE, algorithm = 3, resolution = 0.7)
DimPlot(object = combined, label = TRUE) + ggtitle("RNA UMAP")
SpatialDimPlot(combined, label = TRUE, label.size = 3) + ggtitle("Spatial RNA")


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.3)

DimPlot(object = combined, reduction = "umap.atac", label = TRUE) + ggtitle("ATAC UMAP")
SpatialDimPlot(combined, label = TRUE, label.size = 3) + ggtitle("Spatial ATAC")

DimPlot(object = combined, reduction = "wnn.umap", label = TRUE) + ggtitle("WNN UMAP")
SpatialDimPlot(combined, label = TRUE, label.size = 3) + ggtitle("Spatial WNN")

saveRDS(combined, file = "D:/ych小组会文件/studyingR/P22_project/P22_Seurat_after_agg_new.rds")
combined <- readRDS("D:/ych小组会文件/studyingR/P22_project/P22_Seurat_after_agg_new.rds")


# step2. identify putative cell clusters
source("D:/ych小组会文件/FGOT-zlh/Codes/Identify_putative_clusters.R")
future.globals.maxSize = 10000*1024^2
options(future.globals.maxSize = future.globals.maxSize)

combined <- identifyClusters(combined, graph.name = 'wsnn')

combined$putative_clusters <- Idents(combined)

# step3. identify confident cells
source("D:/ych小组会文件/FGOT-zlh/Codes/identifyConfidentCells.R")
combined <- identifyConfidentCells(combined, graph.name = 'wsnn')
confident.cells <-  Misc(combined, slot = 'confident.cells') # used for anchors
snn_smooth <- SNN <- combined[['wsnn']] # used for compute cost
putative_clusters <- combined$putative_clusters
# write out the result
# data needs: aggregated RNA and ATAC; focus genes; anchors; putative clusters; snn
write.table(confident.cells, file = "D:/ych小组会文件/studyingR/P22_project/P22_anchors_new.txt", quote = FALSE, sep = "\t")
write.table(putative_clusters, file = "D:/ych小组会文件/studyingR/P22_project/P22_putative_clusters_new.txt", quote = FALSE, sep = "\t")
write.table(snn_smooth, file = "D:/ych小组会文件/studyingR/P22_project/P22_wsnn_new.txt", quote = FALSE,sep = "\t")



