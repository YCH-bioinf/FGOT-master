rm(list=ls())
# setwd("/home/nas2/biod/yangchenghui/my_ych_project_PBMC/output/")
setwd("Z:/yangchenghui/my_ych_project_PBMC/output/")

library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(tidyr)

TSS_full <- read.table(file = "Z:/yangchenghui/my_ych_project_PBMC/PBMC_10k/hg38.promoter.regions.txt") 
names(TSS_full) <- c("Chrom","Starts","Ends","genes")
genes <- lapply(TSS_full$genes, function(x) strsplit(x,"[|]")[[1]][1])
TSS_full$genes <- unlist(genes)
unik <- !duplicated(genes)# filter out different transcript
TSS <- TSS_full[unik,]

# import regulatory score from FGOT
Res <- read.csv2("PBMC_regulatory_scores.txt", sep = "\t", row.names = 1, header = TRUE)
Res <- t(Res)
re_g <- colnames(Res)


# for each cell state, build a dataframe
num_clust = nrow(Res)
regulatory_results <- list()
gene_flag <- c()
for (k in 1:num_clust) {
  results_k <- list()
  for (i in 1:length(re_g)) {
    s_i <- unlist(strsplit(re_g[[i]], "-", fixed = TRUE)) ### note _ or -
    gene_i <- s_i[[1]]
    id_i <- which(TSS$genes == gene_i)
    if (length(id_i) > 0){
      # build 
      gene_flag[i] = id_i
      peak_1 = paste0(TSS$Chrom[id_i],"_",TSS$Starts[id_i],"_",TSS$Ends[id_i])
      #peak_2 = gsub('-','_',s_i[2])
      peak_2 = paste0(s_i[2],"_",s_i[3],"_",s_i[4])
      conns_i <- data.frame(Peak1 = peak_1, Peak2 = peak_2, Score = as.double(Res[k,i]), Genes = gene_i)
      results_k[[i]] <- conns_i
    }
    regulatory_results[[k]] <- results_k
  }
}
  
gene_index <- unique(gene_flag[which(gene_flag > 0)])
### focus on these genes
focus_genes <- TSS$genes[gene_index]
TSS <- TSS[gene_index,]
Chr <- TSS$Chrom
Starts <- TSS$Starts


result_all <- list()
for (i in 1:length(regulatory_results)) {
  result_i <- do.call("rbind", regulatory_results[[i]])
  result_i$cluster <- rownames(Res)[i]
  result_all[[i]] <- result_i
}

result_df <- do.call("rbind", result_all)
# write out the dataframe
write.table(result_df, file = "PBMC_cell_cluster_specific_regulatory_score.txt",sep = "\t")



result_df <- read.csv("PBMC_cell_cluster_specific_regulatory_score.txt",sep = "\t")

## plot the regulatory links, echo=TRUE

marker = 'FCMR'

cell_type <- 'Naive.B'
ID <- match(marker, TSS$genes)
locus <- TSS[ID, ,drop = FALSE]
conns <- result_df[intersect(which(result_df$Genes == marker), which(result_df$cluster == cell_type)),,drop = FALSE]
coaccess <- conns$Score
# normalized to 0-1
Tscores <- (coaccess - min (coaccess)) / (max (coaccess) - min (coaccess)) 
conns$coaccess <- Tscores
library(cicero)
data(gene_annotation_sample)
g1 <- plot_connections(conns, locus$Chrom, locus$Starts-250000, locus$Starts+250000,
                       gene_model = gene_annotation_sample,
                       coaccess_cutoff = 0.2,
                       connection_width = 1,
                       connection_color = "#99000d",
                       peak_color = "#B4656F",
                       collapseTranscripts = "longest")

HiC_conns <- readRDS("Z:/yangchenghui/my_ych_project_PBMC/PCHiC_link_ground_truth/PCHiC_ground_truth_nB.rds")
fil_HiC_conns <- HiC_conns %>% filter(Genes == marker)
coaccess <- fil_HiC_conns$Score
Tscores <- (coaccess - min (coaccess)) / (max (coaccess) - min (coaccess)) 
fil_HiC_conns$coaccess <- Tscores

HiC_plot <- plot_connections(fil_HiC_conns, locus$Chrom, locus$Starts-250000, locus$Starts+250000,
                             gene_model = gene_annotation_sample, 
                             coaccess_cutoff = 0.5, 
                             connection_width = 1, 
                             collapseTranscripts = "longest")

fil_HiC_conns <- fil_HiC_conns %>% filter(coaccess > 0.5)
library(GenomicRanges)
conns$in_HiC <- compare_connections(conns, fil_HiC_conns[,c("Peak1","Peak2")], maxgap = 2000)
conns$connection_color <- ifelse(conns$in_HiC, "#92c5de", "lightgray")
g_plot <- plot_connections(conns, locus$Chrom, locus$Starts - 250000, locus$Starts + 250000,
                           gene_model = gene_annotation_sample, 
                           coaccess_cutoff = 0.2, 
                           connection_width = 1,
                           connection_color = "connection_color",  # 使用刚才添加的列
                           peak_color = "#B4656F",
                           collapseTranscripts = "longest")

