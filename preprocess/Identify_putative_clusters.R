#' Identify cell clusters
#'
#' @param object scMC object
#' @param graph.name wsnn or snn
#' @param setIdent whether set the current clustering results as the default cell identity; default is TRUE when performing integrated clustering
#' @param resRange the range of resolution values in Leiden algorithm; if it is NULL, the default range of resoultion will be from 0.1 to 0.5
#' @param resolution the resolution in Leiden algorithm; if it is NULL, the optimal resoultion will be inferred based on eigen spectrum
#' @param do.parallel whther do parallel when automatically inferring the resolution
#'
#'
#'
#' @return object
#' @export

identifyClusters <- function(object, graph.name, setIdent = TRUE, resRange = NULL, resolution = NULL, do.parallel = TRUE){
    object$putative_clusters <- identifyClusters_internal(object, graph.name = graph.name, resRange = resRange, resolution = resolution, do.parallel = do.parallel)
    #object <- addMeta(object, meta = as.matrix(object$putative_clusters), meta.name = "putative clusters")
    #object <- setIdent(object, ident.use = "putative clusters")
    Idents(object) <- "putative_clusters"
  return(object)
}


#' identify cell clusters
#'
#' @param seu Seurat object
#' @param resRange the range of resolutions in Leiden algorithm; if it is NULL, the optimal range of resoultion will be from 0.1 to 0.5
#' @param resolution the resolution in Leiden algorithm; if it is NULL, the optimal resoultion will be inferred based on eigen spectrum
#' @param do.parallel whther do parallel when automatically inferring the resolution
#'
#' @importFrom rsvd rsvd
#' @importFrom stats as.dist hclust cutree
#' @importFrom future plan nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom Matrix Matrix
#'
identifyClusters_internal <- function(seu,graph.name, resRange = NULL, resolution = NULL, do.parallel = TRUE) {
    if (!is.null(resolution)) {
        seu <- suppressWarnings(Seurat::FindClusters(seu, graph.name = graph.name,resolution = resolution, verbose = FALSE))
        #idents <- runLeiden(SNN, resolution)
        idents <- Idents(seu)
    } else {
        N <- nrow(seu@meta.data)
        if (is.null(resRange)) {
            resRange <- seq(0.1,0.5,by = 0.1)
        }
        
        if (do.parallel) {
            #future::plan("multiprocess", workers = 1)
            future::plan("multisession", workers = 1)
        }
        my.sapply <- ifelse(
        test = future::nbrOfWorkers() == 1,
        yes = pbapply::pbsapply,
        no = future.apply::future_sapply
        )
        results = my.sapply(
        X = 1:length(resRange),
        FUN = function(x) {
            seu <- suppressWarnings(Seurat::FindClusters(seu,graph.name = graph.name, resolution = resRange[x], verbose = FALSE))
            idents <- Idents(seu)
            clusIndex <- as.numeric(as.character(idents))
            adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, clusIndex, FUN = "==")), nrow = N, ncol = N)
            return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
        },
        simplify = FALSE
        )
        adjMat <- lapply(results, "[[", 1)
        CM <- Reduce('+', adjMat)/length(resRange)
        if (length(unique(as.vector(CM))) > 1) {
            # dertermine K
            nPC <- min(30,nrow(CM))
            if (length(unique(as.vector(CM))) > nPC) {
                out <- rsvd::rsvd(CM,nPC)
            }
            else {
                out <- svd(CM,nPC)
            }            
            s <- out$d
            # compute ratio
            cs = cumsum(s)/sum(s)
            tol <- 0.01
            K <- min(which((s/sum(s) < tol) == 1))-1
            d <- stats::as.dist(1 - CM)
            hc <- stats::hclust(d, "ave")
            idents<- stats::cutree(hc,k = K)
        }
        else {
            seu <- suppressWarnings(Seurat::FindClusters(seu, graph.name = graph.name, resolution = resRange[1], verbose = FALSE))
            idents <- Idents(seu)
        }
        
    }
    names(idents) <- rownames(seu@meta.data)
    return(idents)
}

#' Group single cells that make up their own cluster in with the cluster they are most connected to.
#
#' @param idents  clustering result
#' @param SNN     SNN graph used in clustering
#' @return        Returns identity with all singletons merged with most connected cluster
#' @export
#'
assignSingletons <- function(idents, SNN) {
  # identify singletons
  singletons <- c()
  for (cluster in unique(idents)) {
    if (length(which(idents %in% cluster)) == 1) {
      singletons <- append(x = singletons, values = cluster)
    }
  }
  #singletons = names(table(idents))[which(table(idents)==1)]
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- unique(x = idents)
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode="numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  for (i in singletons) {
    print(i)
    for (j in cluster_names) {
      subSNN = SNN[
        idents %in% i,
        idents %in% j
        ]
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    # which(x = idents %in% i)[which(x = idents %in% i)] <- closest_cluster
    idents[idents %in% i] <- closest_cluster
  }
  if (length(x = singletons) > 0) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(idents)),
      "final clusters."
    ))
  }
  return(idents)
}


#' identify shared nearest neighbors
#'
#' @param object scMC object
#' @param mode "separate" or "integrate"
#' @param nDims the number of dimensions to use for building SNN
#' @param do.fast whether do fast PCA
#'
#'
#' @return object with SNN
#' @export

identifyNeighbors <- function(object, mode = c("separate","integrated"), nDims = 40, do.fast = TRUE){
  mode <- match.arg(mode)
  if (mode == "separate") {
    for (i in 1:length(object@data)){
      pcaEmbeddings <- runPCA(object@scale.data[[i]], do.fast = do.fast, dimPC = nDims)
      object@snn[[i]] <- buildSNN(pcaEmbeddings)
    }
    names(object@snn) <- names(object@data)
  } else if (mode == "integrated") {
    integratedEmbeddings <- object@dr$integrated
    object@snn[["integrated"]] <- buildSNN(integratedEmbeddings[, 1:min(nDims,ncol(integratedEmbeddings))])
  }
  return(object)
}


#' Build SNN matrix
# #' Adapted from swne (https://github.com/yanwu2014/swne)
#' @param data.use Features x samples matrix to use to build the SNN
#' @param k Defines k for the k-nearest neighbor algorithm
#' @param prune.SNN Sets the cutoff for acceptable Jaccard distances when
#'                  computing the neighborhood overlap for the SNN construction.
#'
#' @return Returns similarity matrix in sparse matrix format
#'
#' @importFrom RANN nn2
#' @export
#'
buildSNN <- function(data.use, k = 30, prune.SNN = 1/15) {
  ## find the k-nearest neighbors for each single cell
  knn <- RANN::nn2(data.use, k = k)
  snn <- ComputeSNN(knn$nn.idx, prune = prune.SNN)
  dimnames(snn) <- list(rownames(data.use), rownames(data.use))
  return(snn)
}


##' Run Leiden clustering algorithm
##  This code is modified from Tom Kelly (https://github.com/TomKellyGenetics/leiden), where we added more parameters (seed.use and n.iter) to run the Python version. In addition, we also take care of the singleton issue after running leiden algorithm.
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##' @param SNN An adjacency matrix compatible with \code{\link[igraph]{igraph}} object.
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
##' @param initial.membership,weights,node.sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None).
##' @param resolution A parameter controlling the coarseness of the clusters
##' @param seed.use rand seed
##' @param n.iter number of iteration
##' @return A parition of clusters as a vector of integers
##' @examples
##' #check if python is availble
##' modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
##' #generate example weighted matrix
# ##' adjacency_matrix[adjacency_matrix == 1] <- weights
# ##' partition <- leiden(adjacency_matrix)
# ##' table(partition)
##' }
##'
##'
##' @keywords graph network igraph mvtnorm simulation
##' @importFrom reticulate import r_to_py
##' @export

runLeiden <- function(SNN = matrix(), resolution = 1, partition_type = c(
  'RBConfigurationVertexPartition',
  'ModularityVertexPartition',
  'RBERVertexPartition',
  'CPMVertexPartition',
  'MutableVertexPartition',
  'SignificanceVertexPartition',
  'SurpriseVertexPartition'
),
seed.use = 42L,
n.iter = 10L,
initial.membership = NULL, weights = NULL, node.sizes = NULL) {
  if (!reticulate::py_module_available(module = 'leidenalg')) {
    stop("Cannot find Leiden algorithm, please install through pip (e.g. pip install leidenalg).")
  }

  #import python modules with reticulate
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  resolution_parameter <- resolution
  initial_membership <- initial.membership
  node_sizes <- node.sizes
  #convert matrix input (corrects for sparse matrix input)
  adj_mat <- as.matrix(SNN)

  #compute weights if non-binary adjacency matrix given
  is_pure_adj <- all(as.logical(adj_mat) == adj_mat)
  if (is.null(weights) && !is_pure_adj) {
    #assign weights to edges (without dependancy on igraph)
    weights <- t(adj_mat)[t(adj_mat)!=0]
    #remove zeroes from rows of matrix and return vector of length edges
  }

  ##convert to python numpy.ndarray, then a list
  adj_mat_py <- r_to_py(adj_mat)
  adj_mat_py <- adj_mat_py$tolist()

  #convert graph structure to a Python compatible object
  GraphClass <- if (!is.null(weights) && !is_pure_adj){
    ig$Graph$Weighted_Adjacency
  } else {
    ig$Graph$Adjacency
  }
  snn_graph <- GraphClass(adj_mat_py)

  #compute partitions
  partition_type <- match.arg(partition_type)
  part <- switch(
    EXPR = partition_type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial.membership, weights = weights,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial.membership, weights = weights,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      resolution_parameter = resolution_parameter,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial.membership,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial.membership, node_sizes = node.sizes,
      resolution_parameter = resolution,
      n_iterations = n.iter,
      seed = seed.use
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial.membership, weights = weights, node_sizes = node.sizes,
      n_iterations = n.iter,
      seed = seed.use
    ),
    stop("please specify a partition type as a string out of those documented")
  )
  partition <- part$membership+1
  idents <- partition

  if (min(table(idents)) == 1) {
    idents <- assignSingletons(idents, SNN)
  }
  idents <- factor(idents)
  names(idents) <- row.names(SNN)
  return(idents)
}



#' Identify marker genes associated with cell clusters
#'
#' @param object scMC object
#' @param mode "separate" or "integrated"; "separate" means identifying marker genes of cell clusters from each dataset separately, and "integrated" means identifying marker genes of cell clusters from the integrated data
#' @param features features used to perform statistical test
#' @param test.use which test to use ("bimod", "wilcox")
#' @param thresh.pc Threshold of the percent of cells enriched in one cluster
#' @param thresh.fc Threshold of Log Fold Change
#' @param thresh.p Threshold of p-values
#'
#'
#' @export
#'
#' @examples
identifyMarkers <- function(object, mode = c("separate","integrated"), features = NULL, test.use = "bimod",
                            thresh.pc = 0.25, thresh.fc = 0.25, thresh.p = 0.05) {
  mode <- match.arg(mode)
  if (mode == "separate") {
    for (i in 1:length(object@data)){
      confident.cells <-  object@confident.cells[[i]]
      subcluster <- object@identity[[i]][confident.cells]
      data.use <- object@data[[i]][, confident.cells]
      group <- as.factor(subcluster)
      object@markers[[i]] <- identifyMarkers_internal(data.use, group, features, test.use,
                                                      thresh.pc, thresh.fc, thresh.p)
    }
    names(object@markers) <- names(object@data)
  } else if (mode == "integrated") {
    data.use <- do.call(cbind, object@data)
    group <- as.factor(object@identity[["integrated"]] )
    object@markers[["integrated"]] <- identifyMarkers_internal(data.use, group, features, test.use,
                                                               thresh.pc, thresh.fc, thresh.p)
  }
  return(object)
}


#' Identify markers in each cell cluster
#'
#' @param data.use normalized data
#' @param group a vector of cell cluster
#' @param features features used for identifying marker genes. default use all features
#' @param test.use which test to use ("wilcox", "bimod")
#' @param thresh.pc Threshold of the percent of cells enriched in one cluster
#' @param thresh.fc Threshold of Log Fold Change
#' @param thresh.p Threshold of p-values
#'
#'
#' @importFrom future nbrOfWorkers
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom stats sd wilcox.test
#' @importFrom stats p.adjust
#'
#'
#' @examples
identifyMarkers_internal <- function(data.use, group, features = NULL, test.use = "bimod",
                                     thresh.pc = 0.25, thresh.fc = 0.25, thresh.p = 0.05) {
  X <- data.use
  if (is.null(features)) {
    features.use <- row.names(X)
  } else {
    features.use <- intersect(features, row.names(X))
  }
  data.use <- X[features.use,]
  data.use <- as.matrix(data.use)

  labels <- group
  level.use <- levels(labels)
  numCluster <- length(level.use)

  my.sapply <- ifelse(
    test = future::nbrOfWorkers() == 1,
    yes = pbapply::pbsapply,
    no = future.apply::future_sapply
  )

  mean.fxn <- function(x) {
    return(log(x = mean(x = expm1(x = x)) + 1))
  }
  genes.de <- vector("list", length = numCluster)
  for (i in 1:numCluster) {
    features <- features.use
    cell.use1 <- which(labels %in% level.use[i])
    cell.use2 <- base::setdiff(1:length(labels), cell.use1)

    # feature selection (based on percentages)
    thresh.min <- 0
    pct.1 <- round(
      x = rowSums(x = data.use[features, cell.use1, drop = FALSE] > thresh.min) /
        length(x = cell.use1),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(x = data.use[features, cell.use2, drop = FALSE] > thresh.min) /
        length(x = cell.use2),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2)
    colnames(x = data.alpha) <- c("pct.1", "pct.2")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > thresh.pc))
    if (length(x = features) == 0) {
      #stop("No features pass thresh.pc threshold")
      next
    }

    # feature selection (based on average difference)
    data.1 <- apply(X = data.use[features, cell.use1, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    data.2 <- apply(X = data.use[features, cell.use2, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
    FC <- (data.1 - data.2)
    features.diff <- names(x = which(x = FC > thresh.fc))
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      #  stop("No features pass thresh.fc threshold")
      next
    }

    data1 <- data.use[features, cell.use1, drop = FALSE]
    data2 <- data.use[features, cell.use2, drop = FALSE]

    if (test.use == "wilcox") {
      pvalues <- unlist(
        x = my.sapply(
          X = 1:nrow(x = data1),
          FUN = function(x) {
            return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
          }
        )
      )

    } else if (test.use == 'bimod') {
      pvalues <- unlist(
        x = my.sapply(
          X = 1:nrow(x = data1),
          FUN = function(x) {
            return(DifferentialLRT(
              x = as.numeric(x = data1[x,]),
              y = as.numeric(x = data2[x,])
            ))
          }
        )
      )
    }

    pval.adj = stats::p.adjust(
      p = pvalues,
      method = "bonferroni",
      n = nrow(X)
    )
    genes.de[[i]] <- data.frame(clusters = level.use[i], features = as.character(rownames(data1)), pvalues = pvalues, pval_adj = pval.adj, logFC = FC[features], data.alpha[features,, drop = F])
  }

  markers.all <- data.frame()
  for (i in 1:numCluster) {
    gde <- genes.de[[i]]
    if (!is.null(gde)) {
      gde <- gde[order(gde$pvalues, -gde$logFC), ]
      gde <- subset(gde, subset = pvalues < thresh.p)
      if (nrow(gde) > 0) {
        markers.all <- rbind(markers.all, gde)
      }
    }

  }
  markers.all$features <- as.character(markers.all$features)
  return(markers.all)
}


# function to run mcdavid et al. DE test
#' likelood ratio test
#'
#' @param x a vector
#' @param y a vector
#' @param xmin threshold for the values in the vector
#'
#'
#' @importFrom stats pchisq
DifferentialLRT <- function(x, y, xmin = 0) {
  lrtX <- bimodLikData(x = x)
  lrtY <- bimodLikData(x = y)
  lrtZ <- bimodLikData(x = c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
}

#' likelood ratio test
#' @importFrom stats sd dnorm
#' @param x a vector
#' @param xmin threshold for the values in the vector

bimodLikData <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- length(x = x2) / length(x = x)
  xal[xal > 1 - 1e-5] <- 1 - 1e-5
  xal[xal < 1e-5] <- 1e-5
  likA <- length(x = x1) * log(x = 1 - xal)
  if (length(x = x2) < 2) {
    mysd <- 1
  } else {
    mysd <- sd(x = x2)
  }
  likB <- length(x = x2) *
    log(x = xal) +
    sum(dnorm(x = x2, mean = mean(x = x2), sd = mysd, log = TRUE))
  return(likA + likB)
}


#' perform PCA
#' @param data.use input scale data
#' @param do.fast whether do fast PCA
#' @param dimPC the number of components to keep in PCA
#' @param weight.by.var whether use weighted pc.scores
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @importFrom irlba irlba
#' @export
#'
runPCA <- function(data.use, do.fast = T, dimPC = 40, seed.use = 42, weight.by.var = T) {
  set.seed(seed = seed.use)
  if (do.fast) {
    dimPC <- min(dimPC, nrow(data.use) - 1)
    pca.res <- irlba::irlba(t(data.use), nv = dimPC)
    sdev <- pca.res$d/sqrt(max(1, ncol(data.use) - 1))
    if (weight.by.var){
      pc.scores <- pca.res$u %*% diag(pca.res$d)
    } else {
      pc.scores <- pca.res$u
    }
  } else {
    dimPC <- min(dimPC, nrow(data.use) - 1)
    pca.res <- stats::prcomp(x = t(data.use), rank. = dimPC)
    sdev <- pca.res$sdev
    if (weight.by.var) {
      pc.scores <- pca.res$x %*% diag(pca.res$sdev[1:dimPC]^2)
    } else {
      pc.scores <- pca.res$x
    }
  }
  rownames(pc.scores) <- colnames(data.use)
  colnames(pc.scores) <- paste0('PC', 1:ncol(pc.scores))
  cell_coords <- pc.scores
  return(cell_coords)
}


#' Run UMAP
#' @param object input data matrix
#' @param reducedDims Which reduced dimensional space ("integrated" or PCA) to use for the UMAP input. Default is "integratedEmbeddings"
#' @param dims Which dimensions to use as input features
#' @param n.neighbors This determines the number of neighboring points used in
#' local approximations of manifold structure. Larger values will result in more
#' global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50.
#' @param n.components The dimension of the space to embed into.
#' @param distance This determines the choice of metric used to measure distance in the input space.
#' @param n.epochs the number of training epochs to be used in optimizing the low dimensional embedding. Larger values result in more accurate embeddings. If NULL is specified, a value will be selected based on the size of the input dataset (200 for large datasets, 500 for small).
#' @param learning.rate The initial learning rate for the embedding optimization.
#' @param min.dist This controls how tightly the embedding is allowed compress points together.
#' Larger values ensure embedded points are moreevenly distributed, while smaller values allow the
#' algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5.
#' @param spread he effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
#' @param set.op.mix.ratio Interpolate between (fuzzy) union and intersection as the set operation used to combine local fuzzy simplicial sets to obtain a global fuzzy simplicial sets.
#' @param local.connectivity The local connectivity required - i.e. the number of nearest neighbors
#' that should be assumed to be connected at a local level. The higher this value the more connected
#' the manifold becomes locally. In practice this should be not more than the local intrinsic dimension of the manifold.
#' @param repulsion.strength Weighting applied to negative samples in low dimensional embedding
#' optimization. Values higher than one will result in greater weight being given to negative samples.
#' @param negative.sample.rate The number of negative samples to select per positive sample in the
#' optimization process. Increasing this value will result in greater repulsive force being applied, greater optimization cost, but slightly more accuracy.
#' @param a More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param b More specific parameters controlling the embedding. If NULL, these values are set automatically as determined by min. dist and spread.
#' @param seed.use Set a random seed. By default, sets the seed to 42.
#' @param metric.kwds,angular.rp.forest,verbose other parameters used in UMAP
#' @import reticulate
#' @export
#'
runUMAP <- function(
  object,
  reducedDims = "integrated",
  dims = 1:40,
  n.neighbors = 30L,
  n.components = 2L,
  distance = "correlation",
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42L,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  verbose = FALSE){
  if (!reticulate::py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn or reticulate::py_install(packages = 'umap-learn')).")
  }
  options(warn=-1)
  dims <- intersect(dims, 1:ncol(object@dr[[reducedDims]]))
  data.use <- object@dr[[reducedDims]][, dims]
  set.seed(seed.use)
  reticulate::py_set_seed(seed.use)
  umap_import <- reticulate::import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(n.neighbors),
    n_components = as.integer(n.components),
    metric = distance,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    verbose = verbose
  )
  Rumap <- umap$fit_transform
  umap_output <- Rumap(data.use)
  colnames(umap_output) <- paste0('UMAP', 1:ncol(umap_output))
  rownames(umap_output) <- rownames(data.use)
  object@dr$umap <- umap_output
  return(object)
}




