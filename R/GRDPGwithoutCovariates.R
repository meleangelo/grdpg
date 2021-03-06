#' (Generalized) Random Dot Product Graph without Covariates
#'
#' Estimate (Generalized) Random Dot Product Graph (GRDPG) without effect of covariates.
#'
#' @import graphstats
#' @import mclust
#'
#' @param A An adjacency matrix.
#' @param link Link function. Could be 'identity' (by default) or 'logit'.
#' @param clusterMethod Method to cluster the estimated latent position. Could be 'GMM' (by default, Gaussian Mixture Model) or 'kmeans'.
#' @param G `G` for \link[mclust]{Mclust} if \code{clusterMethod=='GMM'} or `centers` for \link[stats]{kmeans} if \code{clusterMethod=='kmeans'}. \code{G = 1:9} by default.
#' @param dmax Maximal embeded dimension. 10 by default.
#' @param dhat Embeded dimension. \code{NULL} by default. If \code{NULL}, will be chosen by \href{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf}{profile likelihood}.
#' @param maxit Maximum number of iterations for `\link[irlba]{irlba}`.
#' @param work Working subspace dimension for `\link[irlba]{irlba}`.
#' @param tol Stopping tolerance for `\link[irlba]{irlba}`.
#' @param check Method to check probability matrix. Could be 'BF' (by default, see \link{BFcheck}) or 'Remove' (see \link{Removecheck}).
#' @param plot Whether to show scree plot and latent position. \code{TRUE} by default.
#' @param ... Additional parameters.
#'
#' @return A list containing the following:
#' \describe{
#' \item{`Xhat`}{An `n` by `d` matrix indicating the estimated latent positions, where `n` is the number nodes and `d` is the embeded dimension.}
#' \item{`BXhat`}{A matrix indicating the estimated block probability matrix.}
#' \item{`clusters_cov`}{An `n`-vecotr indicating the block label of each nodes with the effect of covariates where `n` is the number nodes.}
#' \item{`Ipq`}{`Ipq` matrix for (G)RDPG, see \link{getIpq}.}
#' \item{`iter`}{The number of Lanczos iterations carried out. See \link[irlba]{irlba}.}
#' \item{`mprod`}{The total number of matrix vector products carried out. See \link[irlba]{irlba}.}
#' \item{`pp1`}{Screeplot with covariates.}
#' }
#'
#' @examples
#' ## Parameters
#' addCovariates <- FALSE
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
#' ## Generate network (adjacency matrix)
#' B <- generateB(latent, K, d, addCovariates)
#' P <- generateP(latent, d, block_size, addCovariates)
#' A <- generateA(n, P)
#'
#' ## Estimation
#' result <- GRDPGwithoutCovariates(A)
#'
#' ## Evaluation
#' print(adjustedRandIndex(blocks, result$clusters))            # ARI without covariates
#' print(B)                        # True B matrix
#' print(result$BXhat)             # Estimated B matrix
#'
#' ## Visualization
#' plotLatentPosition(matrix(result$Xhat[,1]), blocks, FALSE, latent, K, d)
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


GRDPGwithoutCovariates <- function(A, link = 'identity', clusterMethod = 'GMM', G = 1:9, dmax = 10, dhat = NULL, maxit = 1000, work = 12, tol = 1e-05, check = 'BF', plot = TRUE, ...) {
  if (!(link %in% c('identity', 'logit'))) {
    print("Unrecognized `link`, would use 'identity' by default.")
  }
  if (!(clusterMethod %in% c('GMM', 'kmeans'))) {
    print("Unrecognized `clusterMethod`, would use 'GMM' by default.")
  }
  if (!(check %in% c('BF', 'Remove'))) {
    print("Unrecognized `check`, would use 'BF' by default.")
  }
  if (clusterMethod == 'kmeans' && length(list(...)) < 1) {
    stop("Need to provide `centers` if use kmeans.")
  }

  result <- list()

  cat('\n\n', 'Embedding...')
  embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)
  s <- embedding$D
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1]+1, dhat)
  Ipq <- getIpq(A, dhat)

  result$iter <- embedding$iter
  result$mprod <- embedding$mprod

  # if (link == 'logit') {
  #   Yhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  #   Qhat <- Yhat %*% Ipq %*% t(Yhat)
  #   if (check == 'BF') {
  #     What <- logit(BFcheck(Qhat))
  #   } else {
  #     What <- logit(Removecheck(Qhat))
  #   }
  #   embedding2 <- SpectralEmbedding(What, dmax, maxit = maxit)
  #   s2 <- embedding2$D
  #   dhat2 <- ifelse(is.null(dhat), dimselect(s2)$elbow[1]+1, dhat)
  #   Xhat <- embedding2$X[,1:dhat2] %*% sqrt(diag(s2[1:dhat2], nrow=dhat2, ncol=dhat2))
  #   Ipq <- getIpq(A, dhat2)
  # } else {
  #   Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  # }

  Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  result$Xhat <- Xhat
  result$Ipq <- Ipq

  cat('\n\n', 'Estimating ...')
  if (clusterMethod == 'GMM') {
    model <- Mclust(Xhat, G, verbose = FALSE)
    muhats <- matrix(model$parameters$mean, nrow = dhat)
    BXhat <- t(muhats) %*% Ipq %*% muhats
    if (link == 'logit') {
      if (check == 'BF') {
        BXhat <- logit(BFcheck(BXhat))
      } else {
        BXhat <- logit(Removecheck(BXhat))
      }
    }
    clusters_cov <- getClusters(data.frame(model$z))
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
  }

  result$BXhat <- BXhat
  result$clusters_cov <- clusters_cov

  if (plot) {
    cols <- ncol(A)
    temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
    s1 <- temp1$values
    temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
    s2 <- temp2$values
    tempdat <- data.frame(raw = c(s1,s2)) %>%
      mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
      arrange(desc(s))
    pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
    result$pp1 <- pp1
  }

  cat('\n****************************************************************************\n')
  return(result)
}
