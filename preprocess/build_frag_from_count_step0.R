##################### convert data matrix to fragments file
rm(list=ls())
setwd("/Users/lihuazhang/Documents/projects/FGOT/Datasets/BMMC")
library(data.table)

library(Seurat)
library(ggplot2)
library(patchwork)


require(magrittr)
require(readr)
require(Matrix)
require(tidyr)
require(dplyr)

# peak-bc matrix
mex_dir_path <- "./atac/"

mtx_path <- paste(mex_dir_path, "GSM4829410_healthy_atac_filtered_matrix.mtx", sep = '/')
feature_path <- paste(mex_dir_path, "GSM4829410_healthy_atac_peaks.bed", sep = '/')
barcode_path <- paste(mex_dir_path, "GSM4829410_healthy_atac_filtered_barcodes.tsv", sep = '/')

features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)

mtx <- Matrix::readMM(mtx_path) %>%
  magrittr::set_rownames(features$feature) %>%
  magrittr::set_colnames(barcodes$barcode)

atac_counts <- as.matrix(mtx)
rs <- rowSums(atac_counts)
atac_counts <- atac_counts[which(rs > 3), ] 

fcounts <- as.table(as.matrix(atac_counts))
frag <- as.data.frame(fcounts)
head(frag)

frag$Var1 <- as.character(frag$Var1)
frag$Var2 <- as.character(frag$Var2)
# split rownames
head(frag)
frag$Freq <- floor(frag$Freq)
frag <- frag[which(frag$Freq > 0),]
peaks <- frag$Var1

chrs <- lapply(peaks, function(x) strsplit(x,"[_]")[[1]][1])
chrs <- unlist(chrs)
starts <- lapply(peaks, function(x) strsplit(x,"[_]")[[1]][2])
starts <- as.numeric(unlist(starts))
ends <- lapply(peaks, function(x) strsplit(x,"[_]")[[1]][3])
ends <- as.numeric(unlist(ends))

frag_new <- data.frame(chrs = chrs, starts = starts, ends = ends, cells = frag$Var2, nums = frag[,3])
frag_new[1:5,1:5]
frag_new <- frag_new[order(frag_new$chrs, frag_new$starts, frag_new$ends), ] #需要按照peak name排序

library(readr)
write_tsv(frag_new,"BMMC_scATAC_fragments.tsv", col_names = FALSE, quote = "none")
# zip and tbi
#bgzip BMMC_scATAC_fragments.tsv
# tabix -p bed BMMC_scATAC_fragments.tsv.gz
