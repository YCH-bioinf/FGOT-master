#' Create aggregated data
#'
#' Function to generate aggregated inputs for DIRECT-NET. \code{generate_aggregated_data}
#' takes as input sparse data. This function will aggregate binary accessibility scores (or gene expression)
#' per cell cluster, if they do not overlap any existing group with more than 80% cells.
#'
#' @param object Seurat object.
#' @param cell_coord similarity matrix or dimiension reductions.
#' @param k_neigh Number of cells to aggregate per group.
#' @param size_factor_normalize Logical, should accessibility values be
#'   normalized by size factor?
#' @param verbose Logical, should warning and info messages be printed?
#' @param max_overlap The maximum overlapping ratio of two groups.
#'
#'
#' @return Aggregated data.
#' @export


# generate_aggregated_data <- function (object,cell_coord, k_neigh = 50,atacbinary = TRUE, 
#                            size_factor_normalize = TRUE, verbose = TRUE, max_overlap=0.8) 
# { 
#   if (nrow(cell_coord) > k_neigh) {
#     # Create a k-nearest neighbors map
#     nn_map <- as.data.frame(FNN::knn.index(cell_coord, 
#                                            k = (k_neigh - 1)))
#     row.names(nn_map) <- row.names(cell_coord)
#     nn_map$agg_cell <- 1:nrow(nn_map)
#     good_choices <- 1:nrow(nn_map)
#     
#     if (verbose)  
#       message("Sample cells randomly.")
#     
#     #Sample cells randomly
#     choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
#     chosen <- good_choices[choice]
#     good_choices <- good_choices[good_choices != good_choices[choice]]
#     
#     it <- 0
#     ##Slow (contain calculating of overlapping between cell groups)
#     while (length(good_choices) > 0 & it < nrow(cell_coord)/((1-max_overlap)*k_neigh)) {
#       it <- it + 1
#       choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
#       new_chosen <- c(chosen, good_choices[choice])
#       good_choices <- good_choices[good_choices != good_choices[choice]]
#       cell_sample <- nn_map[new_chosen, ]
#       
#       #calculate overlapping between cell groups
#       combs <- data.frame(1:(nrow(cell_sample) - 1), nrow(cell_sample))
#       shared <- apply(combs, 1, function(x) {    #Slow
#         (k_neigh * 2) - length(unique(as.vector(as.matrix(cell_sample[x, 
#                                                                       ]))))
#       })
#       if (max(shared) < max_overlap * k_neigh) {
#         chosen <- new_chosen
#       }
#     }
#     
#     
#     #if (verbose) {
#     #  cell_sample <- nn_map[chosen, ]
#     #  combs <- t(combn(nrow(cell_sample), 2))
#     #  shared <- apply(combs, 1, function(x) {
#     #    (k_neigh * 2) - length(unique(as.vector(as.matrix(cell_sample[x, 
#     #                                                                  ]))))
#     #  })
#     #  message(paste0("Overlap QC metrics:\nCells per bin: ", 
#     #                 k_neigh, "\nMaximum shared cells bin-bin: ", max(shared), 
#     #                 "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ", 
#     #                 median(shared)))
#       #if (mean(shared)/k_neigh > 0.1) 
#       #  warning("On average, more than 10% of cells are\n                                   shared between paired bins.")
#     #}
#     
#     #aggregating both scRNA-seq and scATAC-seq counts of cells within one group
#     if ("RNA" %in% names(object@assays)) {
#       # rna_old <- as.matrix(object@assays$RNA@counts)
#       rna_old <- as.matrix(object@assays[["RNA"]]@layers[["counts"]])
#       rna_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(rna_old)) %in% cell_sample[x,,drop=FALSE])
#       rna_mask <- Matrix::Matrix(rna_mask)
#       rna_new <- rna_old %*% rna_mask
#       rna_new <- as.matrix(rna_new)
#     }
#     
#     if ("ATAC" %in% names(object@assays)) {
#         #atac_old <- object@assays$ATAC@counts
#         atac_old <- object@assays[["ATAC"]]@layers[["counts"]]
#         # binarize
#         if (atacbinary) {
#             # atac_old[atac_old > 0] <- 1 # if data is large it is not true
#             atac_old <- atac_old > 0
#         }
#         
#         atac_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(atac_old)) %in% cell_sample[x,,drop=FALSE])
#         atac_mask <- Matrix::Matrix(atac_mask)
#         atac_new <- atac_old %*% atac_mask
#         atac_new <- as.matrix(atac_new)
#         
#     }
#     if ("peaks" %in% names(object@assays)) {
#         #atac_old <- object@assays$peaks@counts
#         atac_old <- object@assays[["peaks"]]@layers[["counts"]]
#         # binarize
#         if (atacbinary) {
#             # atac_old[atac_old > 0] <- 1 # if data is large it is not true
#             atac_old <- atac_old > 0
#         }
#         
#         atac_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(atac_old)) %in% cell_sample[x,,drop=FALSE])
#         atac_mask <- Matrix::Matrix(atac_mask)
#         atac_new <- atac_old %*% atac_mask
#         atac_new <- as.matrix(atac_new)
#         
#     }
#     
#     
#   } else {
#     if ("RNA" %in% names(object@assays)) {
#     #rna_old <- as.matrix(object@assays$RNA@counts)
#     rna_old <- as.matrix(object@assays[["RNA"]]@layers[["counts"]])
#     rna_new <- rowSums(rna_old)
#     rna_new <- as.matrix(rna_new) }
#     
#     if ("ATAC" %in% names(object@assays)) {
#         #atac_old <- object@assays$ATAC@counts
#         atac_old <- object@assays[["ATAC"]]@layers[["counts"]]
#         # binarize
#         if (atacbinary) {
#             atac_old[atac_old > 0] <- 1 # if data is large it is not true
#         }
#         atac_new <-rowSums(atac_old)
#         atac_new <- as.matrix(atac_new)
#         }
#     if ("peaks" %in% names(object@assays)) {
#         #atac_old <- object@assays$peaks@counts
#         atac_old <- object@assays[["peaks"]]@layers[["counts"]]
#         # binarize
#         if (atacbinary) {
#             atac_old[atac_old > 0] <- 1 # if data is large it is not true
#         }
#         atac_new <-rowSums(atac_old)
#         atac_new <- as.matrix(atac_new)
#         }
#     
#     
#     cell_sample <- as.data.frame(t(matrix(seq(from = 1, to = nrow(cell_coord)))))
#   }
#   
#   
#   new_data <- list()
#   if ("RNA" %in% names(object@assays)) {
#   new_data$rna <- rna_new }
#   if ("ATAC" %in% names(object@assays)) {
#   new_data$atac <- atac_new
#   }
#   if ("peaks" %in% names(object@assays)) {
#   new_data$atac <- atac_new
#   }
#   
#   new_data$cell_sample <- cell_sample
#   return (new_data)
# }


generate_aggregated_data <- function (object, snn_matrix = NULL, k_neigh = 50,atacbinary = TRUE, 
                                      size_factor_normalize = TRUE, verbose = TRUE) 
{ 
  # Create a k-nearest neighbors map
  n_cell <- nrow(snn_matrix)
  diag(snn_matrix) <- -Inf
  nn_map <- matrix(0, nrow = n_cell, ncol = k_neigh - 1)
  nn_map <- t(apply(snn_matrix, 1, function(row) {
    order(row, decreasing = TRUE)[1:(k_neigh - 1)]
  }))
  nn_map <- as.data.frame(nn_map)
  row.names(nn_map) <- row.names(snn_matrix)

  #aggregating both scRNA-seq and scATAC-seq counts of cells within one group
  if ("RNA" %in% names(object@assays)) {
    # rna_old <- as.matrix(object@assays$RNA@counts)
    rna_old <- as.matrix(object@assays[["RNA"]]@layers[["counts"]])
    rna_mask <- matrix(0, nrow = ncol(rna_old), ncol = ncol(rna_old))
    # 填充掩码矩阵
    for (i in seq_len(ncol(rna_old))) {
      neighbors <- unlist(nn_map[i,])
      rna_mask[neighbors, i] <- 1
      rna_mask[i, i] <- 1 # 包含自身
    }
    # 将掩码矩阵转换为稀疏矩阵
    rna_mask <- Matrix::Matrix(rna_mask, sparse = TRUE)
    
    rna_new <- rna_old %*% rna_mask
    rna_new <- as.matrix(rna_new)
  }

  if ("ATAC" %in% names(object@assays)) {
      #atac_old <- object@assays$ATAC@counts
      atac_old <- object@assays[["ATAC"]]@layers[["counts"]]
      # binarize
      if (atacbinary) {
          # atac_old[atac_old > 0] <- 1 # if data is large it is not true
          atac_old <- atac_old > 0
      }

      atac_mask <- matrix(0, nrow = ncol(atac_old), ncol = ncol(atac_old))
      # 填充掩码矩阵
      for (i in seq_len(ncol(atac_old))) {
        neighbors <- unlist(nn_map[i,])
        atac_mask[neighbors, i] <- 1
        atac_mask[i, i] <- 1 # 包含自身
      }
      # 将掩码矩阵转换为稀疏矩阵
      atac_mask <- Matrix::Matrix(atac_mask, sparse = TRUE)
      
      atac_new <- atac_old %*% atac_mask
      atac_new <- as.matrix(atac_new)

  }
  if ("peaks" %in% names(object@assays)) {
      #atac_old <- object@assays$peaks@counts
      atac_old <- object@assays[["peaks"]]@layers[["counts"]]
      # binarize
      if (atacbinary) {
          # atac_old[atac_old > 0] <- 1 # if data is large it is not true
          atac_old <- atac_old > 0
      }

      atac_mask <- matrix(0, nrow = ncol(atac_old), ncol = ncol(atac_old))
      # 填充掩码矩阵
      for (i in seq_len(ncol(atac_old))) {
        neighbors <- unlist(nn_map[i,])
        atac_mask[neighbors, i] <- 1
        atac_mask[i, i] <- 1 # 包含自身
      }
      # 将掩码矩阵转换为稀疏矩阵
      atac_mask <- Matrix::Matrix(atac_mask, sparse = TRUE)
        
      atac_new <- atac_old %*% atac_mask
      atac_new <- as.matrix(atac_new)

  }

  cell_sample <- nn_map

  new_data <- list()
  if ("RNA" %in% names(object@assays)) {
  new_data$rna <- rna_new }
  if ("ATAC" %in% names(object@assays)) {
  new_data$atac <- atac_new
  }
  if ("peaks" %in% names(object@assays)) {
  new_data$atac <- atac_new
  }

  new_data$cell_sample <- cell_sample
  return (new_data)
}