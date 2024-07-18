#' Function to calculate the size factor for the single-cell RNA-seq data
#'  
#' @importFrom stats median
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts
#' @param locfunc The location function used to find the representive value 
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default), 
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total". 
#'
estimateSizeFactorsForMatrix <- function(counts, locfunc = median, round_exprs=TRUE,  method="mean-geometric-mean-total")
{
  library("slam")
  #if (isSparseMatrix(counts)){
  if (isSparseMatrix(counts)[1]){
    estimateSizeFactorsForSparseMatrix(counts, locfunc = locfunc, round_exprs=round_exprs, method=method)
  }else{
    estimateSizeFactorsForDenseMatrix(counts, locfunc = locfunc, round_exprs=round_exprs,  method=method)
  }
  
}

# Convert a slam matrix to a sparseMatrix
#' @import slam
#' @import Matrix
asSparseMatrix = function (simpleTripletMatrix) {
  retVal = sparseMatrix(i=simpleTripletMatrix[["i"]],
                        j=simpleTripletMatrix[["j"]],
                        x=simpleTripletMatrix[["v"]],
                        dims=c(simpleTripletMatrix[["nrow"]],
                               simpleTripletMatrix[["ncol"]]))
  if (!is.null(simpleTripletMatrix[["dimnames"]]))
    dimnames(retVal) = simpleTripletMatrix[["dimnames"]]
  return(retVal)
}

# Convert a sparseMatrix from Matrix package to a slam matrix
#' @import slam
asSlamMatrix = function (sp_mat) {
  sp <- Matrix::summary(sp_mat)
  simple_triplet_matrix(sp[,"i"], sp[,"j"], sp[,"x"], ncol=ncol(sp_mat), nrow=nrow(sp_mat), dimnames=dimnames(sp_mat))
}

# Convert a sparseMatrix from Matrix package to a slam matrix
#' @import Matrix
isSparseMatrix <- function(x){
  class(x) %in% c("dgCMatrix", "dgTMatrix")
}


# Estimate size factors for each column, given a sparseMatrix from the Matrix
# package
#' @import slam
#' @importFrom stats median
estimateSizeFactorsForSparseMatrix <- function(counts, 
                                               locfunc = median, 
                                               round_exprs=TRUE, 
                                               method="mean-geometric-mean-total"){
  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  CM <- asSlamMatrix(CM)
  
  if (method == "weighted-median"){
    
    log_medians <- rowapply_simple_triplet_matrix(CM, function(cell_expr) { 
      log(locfunc(cell_expr))
    })
    
    weights <- rowapply_simple_triplet_matrix(CM, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })
    
    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowapply_simple_triplet_matrix(CM, function(x) { mean(log(CM)) })
    
    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    stop("Error: method 'median' not yet supported for sparse matrices")
  }else if(method == 'mode'){
    stop("Error: method 'mode' not yet supported for sparse matrices")
  }else if(method == 'geometric-mean-total') {
    cell_total <- col_sums(CM)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- col_sums(CM)
    sfs <- cell_total / exp(mean(log(cell_total)))
  } 
  
  sfs[is.na(sfs)] <- 1 
  sfs   
}

#' @importFrom stats median
estimateSizeFactorsForDenseMatrix <- function(counts, locfunc = median, round_exprs=TRUE, method="mean-geometric-mean-total"){
  
  CM <- counts
  if (round_exprs)
    CM <- round(CM)
  if (method == "weighted-median"){
    log_medians <- apply(CM, 1, function(cell_expr) { 
      log(locfunc(cell_expr))
    })
    
    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })
    
    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <-  weights * (log(cnts) -  log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( mean(norm_cnts) )
    })
  }else if (method == "median-geometric-mean"){
    log_geo_means <- rowMeans(log(CM))
    
    sfs <- apply( CM, 2, function(cnts) {
      norm_cnts <- log(cnts) -  log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      #print (head(norm_cnts))
      exp( locfunc( norm_cnts ))
    })
  }else if(method == "median"){
    row_median <- apply(CM, 1, median)
    sfs <- apply(Matrix::t(Matrix::t(CM) - row_median), 2, median)
  }else if(method == 'mode'){
    sfs <- estimate_t(CM)
  }else if(method == 'geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  }else if(method == 'mean-geometric-mean-total') {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  } 
  
  sfs[is.na(sfs)] <- 1 
  sfs  
}