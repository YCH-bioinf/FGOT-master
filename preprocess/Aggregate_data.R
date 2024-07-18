#' Create aggregated data
#'
#' Function to generate aggregated inputs for DIRECT-NET. \code{generate_aggregated_data}
#' takes as input sparse data. This function will aggregate binary accessibility scores (or gene expression)
#' per cell cluster, if they do not overlap any existing group with more than 80% cells.
#'
#' @param seurat_object Seurat object.
#' @param k_neigh Number of cells to aggregate per group.
#' @param size_factor_normalize Logical, should accessibility values be
#'   normalized by size factor?
#' @param verbose Logical, should warning and info messages be printed?
#' @param max_overlap The maximum overlapping ratio of two groups.
#'
#'
#' @return Aggregated Seurat object.
#' @export
#'
Aggregate_data <- function (seurat_object, group = NULL, k_neigh = 50, atacbinary = TRUE,
                            snn_matrix = NULL, size_factor_normalize = TRUE, verbose = TRUE) 
{
  # if (is.null(cell_coord)) {
  #   if ("ATAC" %in% names(seurat_object@assays)) {
  #     cell_coord <- seurat_object@reductions$umap@cell.embeddings
  #   } else {
  #     cell_coord <- seurat_object@reductions$umap.peaks@cell.embeddings
  #   }
  # } 
  
  # if ("RNA" %in% names(seurat_object@assays)) {
  # rna_new <- matrix(0,nrow = nrow(seurat_object@assays[["RNA"]]@layers[["counts"]]), ncol =1)
  # #rna_new <- matrix(0,nrow = nrow(seurat_object@assays$RNA@counts), ncol =1)
  # }
  # if ("ATAC" %in% names(seurat_object@assays)) {
  # #atac_new <- matrix(0,nrow = nrow(seurat_object@assays$ATAC@counts), ncol =1)
  # atac_new <- matrix(0,nrow = nrow(seurat_object@assays[["ATAC"]]@layers[["counts"]]), ncol =1)
  # }
  # if ("peaks" %in% names(seurat_object@assays)) {
  # # atac_new <- matrix(0,nrow = nrow(seurat_object@assays$peaks@counts), ncol =1)
  # atac_new <- matrix(0,nrow = nrow(seurat_object@assays[["ATAC"]]@layers[["counts"]]), ncol =1)
  # }
  # 
  
  aggregated_data <- generate_aggregated_data(seurat_object,snn_matrix,k_neigh,atacbinary, verbose)
  cell_sample <- aggregated_data$cell_sample
  if ("RNA" %in% names(seurat_object@assays)) {
    rna_new <- aggregated_data$rna
  }
  if ("ATAC" %in% names(seurat_object@assays)) {
    atac_new <- aggregated_data$atac
  }
  if ("peaks" %in% names(seurat_object@assays)) {
    atac_new <- aggregated_data$atac
  }

  # cell_sample <- matrix(0,nrow = 1,ncol = k_neigh)
  # for (i in 1:length(uniqgroup)) {
  #   subobject <- subset(seurat_object, idents = uniqgroup[i])
  #   sub_index <- which(group %in% uniqgroup[i])
  #   cell_coord_i <- cell_coord[sub_index,]
  #   sub_aggregated_data <- generate_aggregated_data(subobject,cell_coord_i,snn_matrix,k_neigh,atacbinary, verbose, max_overlap)
  #   
  #   sub_cell_sample <- sub_aggregated_data$cell_sample
  #   if ("RNA" %in% names(seurat_object@assays)) {
  #   rna_new <- cbind(rna_new,sub_aggregated_data$rna)
  #   }
  #   if ("ATAC" %in% names(seurat_object@assays)) {
  #   atac_new <- cbind(atac_new,sub_aggregated_data$atac)
  #   }
  #   if ("peaks" %in% names(seurat_object@assays)) {
  #   atac_new <- cbind(atac_new,sub_aggregated_data$atac)
  #   }
  #   
  #   if (ncol(sub_cell_sample) < k_neigh) {
  #     sub_cell_sample_new <- as.matrix(sub_cell_sample)
  #     sub_cell_sample_new <- cbind(sub_cell_sample_new,matrix(0,nrow = 1,ncol = k_neigh - ncol(sub_cell_sample_new)))
  #   } else {
  #     sub_cell_sample_new <- apply(sub_cell_sample, 2, function(x) {  
  #       sub_index[x]#for each column return original index
  #     })
  #     sub_cell_sample_new <- as.data.frame(sub_cell_sample_new)
  #     sub_cell_sample_new <- as.matrix(sub_cell_sample_new)
  #   }
  #   cell_sample <- rbind(cell_sample,sub_cell_sample_new)
  # }
  # 
  
  
  ######### normalization
  
  if (size_factor_normalize) {
    if ("RNA" %in% names(seurat_object@assays)) {
    rna_new <- log(rna_new+1)/estimateSizeFactorsForMatrix(rna_new)
    }
    if ("ATAC" %in% names(seurat_object@assays)) {
    atac_new <- log(atac_new+1)/estimateSizeFactorsForMatrix(atac_new)
    }
    if ("peaks" %in% names(seurat_object@assays)) {
    atac_new <- log(atac_new+1)/estimateSizeFactorsForMatrix(atac_new)
    }
  }
  new_data <- list()
  if ("RNA" %in% names(seurat_object@assays)) {
  new_data$rna <- rna_new
  }
  if ("ATAC" %in% names(seurat_object@assays)) {
  new_data$atac <- atac_new
  }
  if ("peaks" %in% names(seurat_object@assays)) {
  new_data$atac <- atac_new
  }
  new_data$cell_sample <- cell_sample
  return (new_data)
}


