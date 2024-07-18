#' identify confident cells in each cell cluster
#'
#' @param object a Seurat object
#' @param quantile.cutoff quantile cutoff (default = 0.75)
#' @param min.cluster the minimum number of cells in the cluster for identifying confident cells
#' @param assay Assay to use
#' @importFrom stats quantile
#' @import Seurat
#'
#' @return a list of Seurat object with identified confident cells in each dataset
#' @export

identifyConfidentCells <- function(object, graph.name, quantile.cutoff = 0.75, min.cluster = 5, assay = NULL) {
  if (is.null(assay)) {assay <- DefaultAssay(object = object)}
  data.norm <- GetAssayData(object = object, assay = assay, slot = "data")
  SNN <- object[[graph.name]]
  #cluster <- Idents(object)
  cluster <- object$putative_clusters
  cluster.uni <- as.character(unique(cluster))
    
  flag <- c()
  for (i in 1:length(cluster.uni)) {
    Index <- which(cluster == cluster.uni[i])
    if (length(Index) > min.cluster) {
      Si <- SNN[Index, Index]
      colsi <- Matrix::colSums(Si)/nrow(Si)
      thresh <- stats::quantile(colsi,quantile.cutoff)
      index <- which(colsi > thresh)
      flag <- c(flag, Index[index])
    } 
  }
  confident_index <- flag
  confident.cells <- colnames(data.norm)[confident_index]
  Misc(object, slot = 'confident.cells') <- confident.cells
  return(object)
}
