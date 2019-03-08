#' (Generalized) Random Dot Product Graph with Covariates
#'
#' Estimate (Generalized) Random Dot Product Graph (GRDPG) with effect of covariates.
#'
#' @import graphstats
#' @import mclust
#'
#' @param A An adjacency matrix.
#' @param covariates Observed covariates. Should be a `n` by `c` matrix or dataframe where `n` is the number of nodes and `c` is the number of covariates.
#' @param link Link function. Could be 'identity' (by default) or 'logit'.
#' @param clusterMethod Method to cluster the estimated latent position. Could be 'GMM' (by default, Gaussian Mixture Model) or 'kmeans'.
#' @param G `G` for \link[mclust]{Mclust} if \code{clusterMethod=='GMM'} or `centers` for \link[stats]{kmeans} if \code{clusterMethod=='kmeans'}. \code{G = 1:9} by default.
#' @param dmax Maximal embeded dimension. 10 by default.
#' @param dhat Embeded dimension. \code{NULL} by default. If \code{NULL}, will be chosen by \href{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf}{profile likelihood}.
#' @param maxit Maximum number of iterations for `\link[irlba]{irlba}`.
#' @param check Method to check probability matrix. Could be 'BF' (by default, see \link{BFcheck}) or 'Remove' (see \link{Removecheck}).
#' @param postAnalysis Whether to do some post analysis such as removing the effect of covariates. \code{TRUE} by default.
#' @param plot Whether to show scree plot and latent position. \code{TRUE} by default.
#' @param ... Additional parameters.
#'
#' @return A list containing the following:
#' \describe{
#' \item{`Xhat`}{An `n` by `d` matrix indicating the estimated latent positions, where `n` is the number nodes and `d` is the embeded dimension.}
#' \item{`betahat`}{A `c`-vector indicating the effect of covariates where `c` is the number of covariates.}
#' \item{`BXhat`}{A matrix indicating the estimated block probability matrix.}
#' \item{`clusters_cov`}{An `n`-vecotr indicating the block label of each nodes with the effect of covariates where `n` is the number nodes.}
#' \item{`Ipq`}{`Ipq` matrix for (G)RDPG, see \link{getIpq}.}
#' \item{`Xhatprime`}{If \code{postAnalysis==TRUE}, estimated latent positions after removing the effect of covariates.}
#' \item{`clusters`}{If \code{postAnalysis==TRUE}, block label of each nodes without the effect of covariates.}
#' \item{`pp1`}{Screeplot with covariates.}
#' \item{`pp2`}{Latent position in 2D.}
#' \item{`pp3`}{If \code{postAnalysis==TRUE}, screeplot without covariates.}
#' }
#'
#' @examples
#' ## Parameters
#' addCovariates <- TRUE
#' n <- 2000                  # Number of nodes
#' K <- 2                     # Number of blocks
#' d <- 1                     # Dimension of latent position
#' latent <- cbind(0.1, 0.7)  # Latent position
#'
#' # Balanced case
#' pi <- rep(1/K, K)
#' block_size <- round(pi * n)
#'
#' # Block label without effect of covariates
#' blocks <- c()
#' for (k in 1:length(block_size)) {
#'   blocks <- c(blocks, rep(k, block_size[k]))
#' }
#'
#' ## One binary covariate
#' beta <- 0.3
#' cov <- 2
#' pi_cov <- rep(1/(K*prod(cov)), K*prod(cov))
#' block_size_cov <- round(pi_cov * n)
#' covariates <- matrix(rep(1:cov, each = 500, times = K))
#'
#' # Block label with effect of covariates
#' blocks_cov <- c()
#' for (k in 1:length(block_size_cov)) {
#'   blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
#' }
#'
#' ## Generate network (adjacency matrix)
#' B <- generateB(latent, K, d, addCovariates, cov, beta)
#' P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
#' A <- generateA(n, P)
#'
#' ## Estimation
#' result <- GRDPGwithCovariates(A, covariates)
#'
#' ## Evaluation
#' print(adjustedRandIndex(blocks_cov, result$clusters_cov))    # ARI with covariates
#' print(adjustedRandIndex(blocks, result$clusters))            # ARI without covariates
#' print(beta)                     # True beta
#' print(result$betahat)           # Estimated beta
#' print(B)                        # True B matrix
#' print(result$BXhat)             # Estimated B matrix
#'
#' ## Visualization
#' pp2 <- plotLatentPosition(result$Xhat, blocks, withCovariates = TRUE, dhat = ncol(result$Xhat), covariates)
#' pp4 <- plotLatentPosition(result$Xhatprime, blocks, withCovariates = FALSE, latent, K, d)
#' multiplot(result$pp1, result$pp3, pp2[[1]], pp4, cols = 2)
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


GRDPGwithCovariates <- function(A, covariates, link = 'identity', clusterMethod = 'GMM', G = 1:9, dmax = 10, dhat = NULL, maxit = 1000, check = 'BF', postAnalysis = TRUE, plot = TRUE, ...) {
  if (nrow(A) != nrow(covariates) || ncol(A) != nrow(covariates)) {
    stop("The number of rows/columns in `A` should equal to the number of rows in `covariates`.")
  }
  if (!(link %in% c('identity', 'logit'))) {
    print("Unrecognized `link`, would use 'identity' by default.")
  }
  if (!(clusterMethod %in% c('GMM', 'kmeans'))) {
    print("Unrecognized `clusterMethod`, would use 'GMM' by default.")
  }
  if (!(check %in% c('BF', 'Remove'))) {
    print("Unrecognized `check`, would use 'BF' by default.")
  }


  result <- list()

  cat('\n\n', 'Embedding...')
  temp <- eigen(A)
  embedding <- embed(A, dmax, maxit = maxit)
  s <- embedding$D
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow+1, dhat)
  Ipq <- getIpq(temp$values, dhat)

  # if (link == 'logit') {
  #   Yhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  #   Qhat <- Yhat %*% Ipq %*% t(Yhat)
  #   if (check == 'BF') {
  #     What <- logit(BFcheck(Qhat))
  #   } else {
  #     What <- logit(Removecheck(Qhat))
  #   }
  #   embedding2 <- embed(What, dmax, maxit = maxit)
  #   s2 <- embedding2$D
  #   dhat2 <- ifelse(is.null(dhat), dimselect(s2)$elbow+1, dhat)
  #   Xhat <- embedding2$X[,1:dhat2] %*% sqrt(diag(s2[1:dhat2], nrow=dhat2, ncol=dhat2))
  #   Ipq <- getIpq(temp$values, dhat2)
  # } else {
  #   Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  # }

  Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  result$Xhat <- Xhat
  result$Ipq <- Ipq

  cat('\n\n', 'Estimating beta...')
  cov <- apply(covariates, 2, function(x) length(unique(x)))
  if (clusterMethod == 'GMM') {
    model <- Mclust(Xhat, G, verbose = FALSE)
    muhats <- model$parameters$mean
    BXhat <- t(muhats) %*% Ipq %*% muhats
    if (link == 'logit') {
      if (check == 'BF') {
        BXhat <- logit(BFcheck(BXhat))
      } else {
        BXhat <- logit(Removecheck(BXhat))
      }
    }
    clusters_cov <- getClusters(data.frame(model$z))
    covariates_block <- getBlockCovariates(covariates, clusters_cov)
  } else {
    centers <- max(G)
    model <- kmeans(Xhat, centers)
    muhats <- model$centers
    BXhat <- muhats %*% Ipq %*% t(muhats)
    if (link == 'logit') {
      if (check == 'BF') {
        BXhat <- logit(BFcheck(BXhat))
      } else {
        BXhat <- logit(Removecheck(BXhat))
      }
    }
    clusters_cov <- model$cluster
    covariates_block <- getBlockCovariates(covariates, clusters_cov)
  }
  betahats <- estimatebeta(BXhat, cov, covariates_block)
  betahat <- sapply(betahats, mean)

  result$BXhat <- BXhat
  result$betahat <- betahat
  result$clusters_cov <- clusters_cov

  if (postAnalysis) {
    cat('\n\n', 'Post Analysis...')
    if (link == 'logit') {
      Yhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
      Qhat <- Yhat %*% Ipq %*% t(Yhat)
      if (check == 'BF') {
        What <- logit(BFcheck(Qhat))
      } else {
        What <- logit(Removecheck(Qhat))
      }
      Aprime <- getAwithoutCovariates(What, betahat, covariates)
    } else {
      Aprime <- getAwithoutCovariates(A, betahat, covariates)
    }
    embedprime <- embed(Aprime, dmax, maxit = maxit)
    sprime <- embedprime$D
    dhatprime <- dimselect(sprime)$elbow
    Xhatprime <- embedprime$X[,1:dhatprime] %*% sqrt(diag(sprime[1:dhatprime], nrow=dhatprime, ncol=dhatprime))
    if (clusterMethod == 'GMM') {
      model2 <- Mclust(Xhatprime, G, verbose = FALSE)
      clusters <- getClusters(data.frame(model2$z))
    } else {
      model2 <- kmeans(Xhatprime, centers)
      clusters <- model2$cluster
    }
    result$Xhatprime <- Xhatprime
    result$clusters <- clusters
  }

  if (plot) {
    pp1 <- screeplot(s, 'Screeplot (with Covariates)')
    pp2 <- plotLatentPosition(Xhat, withCovariates = TRUE, dhat = dhat, covariates = covariates)
    if (postAnalysis) {
      pp3 <- screeplot(sprime, 'Screeplot (without Covariates)')
      result$pp3 <- pp3
    }
    multiplot(plotlist = pp2, cols = ceiling(length(pp2)/2))
    result$pp1 <- pp1
    result$pp2 <- pp2
  }

  cat('\n****************************************************************************\n')
  return(result)
}
