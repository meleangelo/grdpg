#' (Generalized) Random Dot Product Graph with Covariates
#'
#' Simulate and estimate (Generalized) Random Dot Product Graph (GRDPG) with effect of covariates.
#'
#' @import graphstats
#' @import mclust
#'
#' @param n Number of nodes.
#' @param K Number of blocks.
#' @param d Dimension of latent position.
#' @param latent True latent position. Suould be a `d` by `K` matrix where `d` is the dimension of latent position and `K` is the number of blocks.
#' @param block_size The size of each block without the effect of covariates.
#' @param beta The effect of covariates. Should be a `c`-vector where `c` is the number of covariates.
#' @param cov Possible values of each covariate could take (assuming discrete here). Should be a `c`-vector where `c` is the number of covariates. For example, if two binary covariates, then \code{cov <- c(2, 2)}.
#' @param block_size_cov The size of each block with the effect of covariates.
#' @param covariates Observed covariates. Should be a `n` by `c` matrix or dataframe where `n` is the number of nodes and `c` is the number of covariates.
#' @param link Link function. Could be 'identity' (by default) or 'logit'.
#' @param clusterMethod Method to cluster the estimated latent position. Could be 'GMM' (by default, Gaussian Mixture Model) or 'kmeans'.
#' @param G `G` for \link[mclust]{Mclust} if \code{clusterMethod=='GMM'} or `centers` for \link[stats]{kmeans} if \code{clusterMethod=='kmeans'}. \code{G = 1:9} by default.
#' @param dmax Maximal embeded dimension. 10 by default.
#' @param dhat Embeded dimension. \code{NULL} by default. If \code{NULL}, will be chosen by \href{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf}{profile likelihood}.
#' @param maxit Maximum number of iterations for `\link[irlba]{irlba}`.
#' @param work Working subspace dimension for `\link[irlba]{irlba}`.
#' @param tol Stopping tolerance for `\link[irlba]{irlba}`.
#' @param check Method to check probability matrix. Could be 'BF' (by default, see \link{BFcheck}) or 'Remove' (see \link{Removecheck}).
#' @param postAnalysis Whether to do some post analysis such as removing the effect of covariates. \code{TRUE} by default.
#' @param plot Whether to show scree plot and latent position. \code{TRUE} by default.
#' @param seed Random seed for reproductivity. 2019 by default.
#' @param ... Additional parameters.
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
#' ## One binary covariate
#' beta <- 0.3
#' cov <- 2
#' pi_cov <- rep(1/(K*prod(cov)), K*prod(cov))
#' block_size_cov <- round(pi_cov * n)
#' covariates <- matrix(rep(1:cov, each = 500, times = K))
#'
#' ## Simulation
#' simulation(n, K, d, latent, block_size, beta, cov, block_size_cov, covariates)
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


simulation <- function(n, K, d, latent, block_size, beta, cov, block_size_cov, covariates, link = 'identity', clusterMethod = 'GMM', G = 1:9, dmax = 10, dhat = NULL, maxit = 1000, work = 12, tol = 1e-05, check = 'BF', postAnalysis = TRUE, plot = TRUE, seed = 2019, ...) {
  cat('\n\n', 'Simulation: (G)RDPG with Covariates', '\n\n\n', 'Setting Up....')

  ## Set Up
  set.seed(seed)
  blocks <- c()
  for (k in 1:length(block_size)) {
    blocks <- c(blocks, rep(k, block_size[k]))
  }
  blocks_cov <- c()
  for (k in 1:length(block_size_cov)) {
    blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
  }
  addCovariates <- TRUE

  ## Generate network (adjacency matrix)
  cat('\n\n', 'Sampling...')
  B <- generateB(latent, K, d, addCovariates, cov, beta)
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
  if (link == 'logit') {
    A <- generateA(n, sigmoid(P))
  } else {
    A <- generateA(n, P)
  }

  ## Estimation
  result <- GRDPGwithCovariates(A, covariates, link, clusterMethod, G, dmax, dhat, maxit, work, tol, check, postAnalysis, plot)

  ## Visualization
  pp2 <- plotLatentPosition(result$Xhat, blocks, withCovariates = TRUE, dhat = ncol(result$Xhat), covariates)
  pp4 <- plotLatentPosition(result$Xhatprime, blocks, withCovariates = FALSE, latent, K, d)
  multiplot(result$pp1, result$pp3, pp2[[1]], pp4, cols = 2)

  ## Evaluation
  cat('\n\nSummary\n\n')
  cat('****************************************************************************\n')
  cat('Latent:\n')
  print(latent)
  cat('K:', K, '\nn:', n, '\nbeta:', beta, '\nbetahat:', result$betahat, '\n')
  cat('B:\n')
  print(B)
  cat('BXhat:\n')
  print(result$BXhat)
  cat('ARI with covariates:', adjustedRandIndex(blocks_cov, result$clusters_cov), '\n')
  cat('ARI without covariates:', adjustedRandIndex(blocks, result$clusters), '\n')
  cat('****************************************************************************\n')
}


