rm(list=ls())
setwd("Z:/yangchenghui/FGOT_scMultiSim_simulation/simulation_data_discrete_4types/")

# devtools::install_github("ZhangLabGT/scMultiSim@main")

library(scMultiSim)
# scmultisim_help("options")
# run_shiny()

data(GRN_params_100)
length(unique(GRN_params_100$regulated.gene)) #94
length(unique(GRN_params_100$regulator.gene)) #6
common_genes <- intersect(GRN_params_100$regulator.gene, GRN_params_100$regulated.gene)
# data(GRN_params_1139)

common_genes <- intersect(GRN_params_100$regulator.gene, GRN_params_100$regulated.gene)

# rand.seed设置随机种子，用于确保模拟结果的可重复性
# GRN指定使用的基因调控网络
# num.cells设置模拟的细胞数量
# num.cifs设置CIF（细胞身份因子）的数量
# cif.sigma控制CIF的标准差，影响细胞群体结构的清晰度
# tree选择用于模拟的分化树，可选Phyla5(),Phyla3(),Phyla1(),ape::read.tree(text = "(A:1,B:1,C:1,D:1);")
# diff.cif.fraction控制受细胞群体影响的CIF比例
# do.velocity设置是否使用动力学模型模拟RNA速度数据
# discrete.cif设置为TRUE表示模拟离散细胞群体
# region.distrib  (default: c(0.1, 0.5, 0.4))
## celltype1

results1 <- sim_true_counts(list(
  rand.seed = 1,
  GRN = GRN_params_100,
  tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
  num.cells = 300, # 细胞数量
  num.cif = 20,
  cif.sigma = 0.1,
  discrete.cif = T,
  atac.effect = 0.9
))

plot_tsne(log2(results1$counts + 1),
          results1$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(results1$atacseq_data + 1),
          results1$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True ATAC-seq Tsne')
saveRDS(results1, file = "type1_scMultiSim_result.rds")

cell_ids_pop1 <- results1[["cell_meta"]]$cell_id[results1[["cell_meta"]]$pop == 1]
cell_ids_pop1 <- paste0("type1_", cell_ids_pop1)
  
rna_count1 <- results1[["counts"]]
colnames(rna_count1) <- paste0("type1_cell", 1:ncol(rna_count1))
rownames(rna_count1) <- ifelse(grepl("^gene", rownames(rna_count1)),
                              rownames(rna_count1),
                              paste0("gene", rownames(rna_count1)))
atac_count1 <- results1[["atac_counts"]]
colnames(atac_count1) <- colnames(rna_count1)
rownames(atac_count1) <- paste0("peak", seq_len(nrow(atac_count1)))

rna_count1 <- rna_count1[, cell_ids_pop1]
atac_count1 <- atac_count1[, cell_ids_pop1]
dim(rna_count1) #110genes * cells
dim(atac_count1) #330peaks* cells

region_gene_mt1 <- results1[["region_to_gene"]]
rownames(region_gene_mt1) <- rownames(atac_count1)
colnames(region_gene_mt1) <- rownames(rna_count1)
dim(region_gene_mt1) # 330 110
length(which(region_gene_mt1 != 0, arr.ind = TRUE)) #296

meta1 <- results1[["cell_meta"]]
meta1$cell_id <- paste0("type1_", meta1$cell_id)
library(dplyr)
pop_to_celltype <- c(
  "1" = "celltype1",
  "2" = "celltype2",
  "3" = "celltype3",
  "4" = "celltype4"
)
meta1 <- meta1 %>% mutate(celltype = pop_to_celltype[pop])
meta1 <- meta1[meta1$cell_id %in% cell_ids_pop1, c("cell_id", "celltype")]

### celltype2
results2 <- sim_true_counts(list(
  rand.seed = 2,
  GRN = GRN_params_100,
  tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
  num.cells = 300, # 细胞数量
  num.cif = 20,
  cif.sigma = 0.1,
  discrete.cif = T,
  atac.effect = 0.9
))

plot_tsne(log2(results2$counts + 1),
          results2$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(results2$atacseq_data + 1),
          results2$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True ATAC-seq Tsne')
saveRDS(results2, file = "type2_scMultiSim_result.rds")

cell_ids_pop2 <- results2[["cell_meta"]]$cell_id[results2[["cell_meta"]]$pop == 2]
cell_ids_pop2 <- paste0("type2_", cell_ids_pop2)

rna_count2 <- results2[["counts"]]
colnames(rna_count2) <- paste0("type2_cell", 1:ncol(rna_count2))
rownames(rna_count2) <- ifelse(grepl("^gene", rownames(rna_count2)),
                               rownames(rna_count2),
                               paste0("gene", rownames(rna_count2)))
atac_count2 <- results2[["atac_counts"]]
colnames(atac_count2) <- colnames(rna_count2)
rownames(atac_count2) <- paste0("peak", seq_len(nrow(atac_count2)))

rna_count2 <- rna_count2[, cell_ids_pop2]
atac_count2 <- atac_count2[, cell_ids_pop2]
dim(rna_count2) #110genes * cells
dim(atac_count2) #330peaks* cells

region_gene_mt2 <- results2[["region_to_gene"]]
rownames(region_gene_mt2) <- rownames(atac_count2)
colnames(region_gene_mt2) <- rownames(rna_count2)
dim(region_gene_mt2) # 330 110
length(which(region_gene_mt2 != 0, arr.ind = TRUE)) #276

meta2 <- results2[["cell_meta"]]
meta2$cell_id <- paste0("type2_", meta2$cell_id)
library(dplyr)
pop_to_celltype <- c(
  "1" = "celltype1",
  "2" = "celltype2",
  "3" = "celltype3",
  "4" = "celltype4"
)
meta2 <- meta2 %>% mutate(celltype = pop_to_celltype[pop])
meta2 <- meta2[meta2$cell_id %in% cell_ids_pop2, c("cell_id", "celltype")]

## 合并rna_count，atac_count, meta
rna_count <- cbind(rna_count1,rna_count2)
atac_count <- cbind(atac_count1,atac_count2)
meta <- rbind(meta1,meta2)

plot_tsne(log2(rna_count + 1),
          meta$celltype, perplexity = 10,
          legend = 'celltype', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(atac_count + 1),
          meta$celltype, perplexity = 10,
          legend = 'pop', plot.name = 'True ATAC-seq Tsne')


### celltype3
results3 <- sim_true_counts(list(
  rand.seed = 3,
  GRN = GRN_params_100,
  tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
  num.cells = 300, # 细胞数量
  num.cif = 20,
  cif.sigma = 0.1,
  discrete.cif = T,
  atac.effect = 0.9
))

plot_tsne(log2(results3$counts + 1),
          results3$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(results3$atacseq_data + 1),
          results3$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True ATAC-seq Tsne')
saveRDS(results3, file = "type3_scMultiSim_result.rds")

cell_ids_pop3 <- results3[["cell_meta"]]$cell_id[results3[["cell_meta"]]$pop == 3]
cell_ids_pop3 <- paste0("type3_", cell_ids_pop3)

rna_count3 <- results3[["counts"]]
colnames(rna_count3) <- paste0("type3_cell", 1:ncol(rna_count3))
rownames(rna_count3) <- ifelse(grepl("^gene", rownames(rna_count3)),
                               rownames(rna_count3),
                               paste0("gene", rownames(rna_count3)))
atac_count3 <- results3[["atac_counts"]]
colnames(atac_count3) <- colnames(rna_count3)
rownames(atac_count3) <- paste0("peak", seq_len(nrow(atac_count3)))

rna_count3 <- rna_count3[, cell_ids_pop3]
atac_count3 <- atac_count3[, cell_ids_pop3]
dim(rna_count3) #110genes * cells
dim(atac_count3) #330peaks* cells

region_gene_mt3 <- results3[["region_to_gene"]]
rownames(region_gene_mt3) <- rownames(atac_count3)
colnames(region_gene_mt3) <- rownames(rna_count3)
dim(region_gene_mt3) # 330 110
length(which(region_gene_mt3 != 0, arr.ind = TRUE)) #286

meta3 <- results3[["cell_meta"]]
meta3$cell_id <- paste0("type3_", meta3$cell_id)
library(dplyr)
pop_to_celltype <- c(
  "1" = "celltype1",
  "2" = "celltype2",
  "3" = "celltype3",
  "4" = "celltype4"
)
meta3 <- meta3 %>% mutate(celltype = pop_to_celltype[pop])
meta3 <- meta3[meta3$cell_id %in% cell_ids_pop3, c("cell_id", "celltype")]

## 合并rna_count，atac_count, meta
rna_count <- cbind(rna_count,rna_count3)
atac_count <- cbind(atac_count,atac_count3)
meta <- rbind(meta,meta3)

plot_tsne(log2(rna_count + 1),
          meta$celltype, perplexity = 10,
          legend = 'celltype', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(atac_count + 1),
          meta$celltype, perplexity = 10,
          legend = 'pop', plot.name = 'True ATAC-seq Tsne')


### celltype4
results4 <- sim_true_counts(list(
  rand.seed = 4,
  GRN = GRN_params_100,
  tree = ape::read.tree(text = "(A:1,B:1,C:1,D:1);"),
  num.cells = 300, # 细胞数量
  num.cif = 20,
  cif.sigma = 0.1,
  discrete.cif = T,
  atac.effect = 0.9
))

plot_tsne(log2(results4$counts + 1),
          results4$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(results4$atacseq_data + 1),
          results4$cell_meta$pop, perplexity = 10,
          legend = 'pop', plot.name = 'True ATAC-seq Tsne')
saveRDS(results4, file = "type4_scMultiSim_result.rds")

cell_ids_pop4 <- results4[["cell_meta"]]$cell_id[results4[["cell_meta"]]$pop == 4]
cell_ids_pop4 <- paste0("type4_", cell_ids_pop4)

rna_count4 <- results4[["counts"]]
colnames(rna_count4) <- paste0("type4_cell", 1:ncol(rna_count4))
rownames(rna_count4) <- ifelse(grepl("^gene", rownames(rna_count4)),
                               rownames(rna_count4),
                               paste0("gene", rownames(rna_count4)))
atac_count4 <- results4[["atac_counts"]]
colnames(atac_count4) <- colnames(rna_count4)
rownames(atac_count4) <- paste0("peak", seq_len(nrow(atac_count4)))

rna_count4 <- rna_count4[, cell_ids_pop4]
atac_count4 <- atac_count4[, cell_ids_pop4]
dim(rna_count4) #110genes * cells
dim(atac_count4) #330peaks* cells

region_gene_mt4 <- results4[["region_to_gene"]]
rownames(region_gene_mt4) <- rownames(atac_count4)
colnames(region_gene_mt4) <- rownames(rna_count4)
dim(region_gene_mt4) # 330 110
length(which(region_gene_mt4 != 0, arr.ind = TRUE)) #286
identical(region_gene_mt3, region_gene_mt4) #FALSE

meta4 <- results4[["cell_meta"]]
meta4$cell_id <- paste0("type4_", meta4$cell_id)
library(dplyr)
pop_to_celltype <- c(
  "1" = "celltype1",
  "2" = "celltype2",
  "3" = "celltype3",
  "4" = "celltype4"
)
meta4 <- meta4 %>% mutate(celltype = pop_to_celltype[pop])
meta4 <- meta4[meta4$cell_id %in% cell_ids_pop4, c("cell_id", "celltype")]

## 合并rna_count，atac_count, meta
rna_count <- cbind(rna_count,rna_count4)
atac_count <- cbind(atac_count,atac_count4)
meta <- rbind(meta,meta4)

plot_tsne(log2(rna_count + 1),
          meta$celltype, perplexity = 10,
          legend = 'celltype', plot.name = 'True RNA Counts Tsne')
plot_tsne(log2(atac_count + 1),
          meta$celltype, perplexity = 10,
          legend = 'pop', plot.name = 'True ATAC-seq Tsne')

########## save results
write.table(rna_count, file = "rna_count.txt", sep = "\t",  quote = FALSE)
write.table(atac_count, file = "atac_count.txt", sep = "\t",  quote = FALSE)
write.table(meta, file = "celltype_info.txt", sep = "\t",  quote = FALSE)

write.table(region_gene_mt1, file = "peak_gene_mt1.txt", sep = "\t",  quote = FALSE)
write.table(region_gene_mt2, file = "peak_gene_mt2.txt", sep = "\t",  quote = FALSE)
write.table(region_gene_mt3, file = "peak_gene_mt3.txt", sep = "\t",  quote = FALSE)
write.table(region_gene_mt4, file = "peak_gene_mt4.txt", sep = "\t",  quote = FALSE)



####umap
library(Seurat)
library(ggplot2)

obj <- CreateSeuratObject(counts = rna_count, assay = "RNA")
obj[["RNA"]]
# Perform standard analysis of each modality independently RNA analysis
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:30)

obj <- AddMetaData(obj, metadata = meta)
DimPlot(obj, reduction = "umap", group.by = "celltype", label = TRUE) + ggtitle("RNA UMAP by Cell Type")

####获得 wnn cost
obj_atac <- CreateSeuratObject(counts = atac_count, assay = "ATAC")
obj_atac[["ATAC"]]
obj_atac <- NormalizeData(obj_atac)
obj_atac <- FindVariableFeatures(obj_atac)
obj_atac <- ScaleData(obj_atac)
obj_atac <- RunPCA(obj_atac, reduction.name = "atac_pca")
obj_atac <- RunUMAP(obj_atac, reduction = 'atac_pca', dims = 1:30)

obj_atac <- AddMetaData(obj_atac, metadata = meta)
DimPlot(obj_atac, reduction = "umap", group.by = "celltype", label = TRUE) + ggtitle("ATAC UMAP by Cell Type")

obj@reductions[["atac_pca"]] <- obj_atac@reductions[["atac_pca"]]
obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "atac_pca"), dims.list = list(1:50, 1:50))
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn",resolution = 0.5,algorithm = 3, verbose = FALSE)
length(unique(Idents(obj)))
DimPlot(obj, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

saveRDS(obj, file = "simu_wnn_seurat_obj.rds")
wnn_smooth <-  obj[['wsnn']]
write.table(wnn_smooth, file = "simu_wnn.txt", quote = FALSE,sep = "\t")

