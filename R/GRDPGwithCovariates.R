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
#' @param work Working subspace dimension for `\link[irlba]{irlba}`.
#' @param tol Stopping tolerance for `\link[irlba]{irlba}`.
#' @param check Method to check probability matrix. Could be 'BF' (by default, see \link{BFcheck}) or 'Remove' (see \link{Removecheck}).
#' @param sd Whether to compute standard errors of the estimate of beta. \code{TRUE} by default.
#' @param rho Sparsity coefficient. Coule be `1` (by default) or `0`.
#' @param postAnalysis Whether to do some post analysis such as removing the effect of covariates. \code{TRUE} by default.
#' @param plot Whether to show scree plot and latent position. \code{TRUE} by default.
#' @param ... Additional parameters.
#'
#' @return A list containing the following:
#' \describe{
#' \item{`Xhat`}{An `n` by `d` matrix indicating the estimated latent positions, where `n` is the number nodes and `d` is the embeded dimension.}
#' \item{`betahat_simple`}{A `c`-vector indicating the effect of covariates where `c` is the number of covariates using simple procedure.}
#' \item{`betahat_simple_unbiased`}{A `c`-vector indicating the effect (unbiased estimate) of covariates where `c` is the number of covariates using simple procedure.}
#' \item{`sd2_simple`}{A `c`-vector indicating the variance of `betahat_simple` where `c` is the number of covariates using simple procedure.}
#' \item{`betahat_weighted`}{A `c`-vector indicating the effect of covariates where `c` is the number of covariates using weighted procedure.}
#' \item{`betahat_weighted_unbiased`}{A `c`-vector indicating the effect (unbiased estimate) of covariates where `c` is the number of covariates using weighted procedure.}
#' \item{`sd2_weighted`}{A `c`-vector indicating the variance of `betahat_weighted` where `c` is the number of covariates using weighted procedure.}
#' \item{`BXhat`}{A matrix indicating the estimated block probability matrix.}
#' \item{`clusters_cov`}{An `n`-vecotr indicating the block label of each nodes with the effect of covariates where `n` is the number nodes.}
#' \item{`Ipq`}{`Ipq` matrix for (G)RDPG, see \link{getIpq}.}
#' \item{`iter`}{The number of Lanczos iterations carried out. See \link[irlba]{irlba}.}
#' \item{`mprod`}{The total number of matrix vector products carried out. See \link[irlba]{irlba}.}
#' \item{`Xhatprime_simple`}{If \code{postAnalysis==TRUE}, estimated latent positions after removing the effect of covariates associated with `betahat_simple`.}
#' \item{`clusters_simple`}{If \code{postAnalysis==TRUE}, block label of each nodes without the effect of covariates associated with `betahat_simple`.}
#' \item{`iterprime_simple`}{If \code{postAnalysis==TRUE}, the number of Lanczos iterations carried out after removing the effect of covariates associated with `betahat_simple`. See \link[irlba]{irlba}.}
#' \item{`mprodprime_simple`}{If \code{postAnalysis==TRUE}, the total number of matrix vector products carried out after removing the effect of covariates associated with `betahat_simple`. See \link[irlba]{irlba}.}
#' \item{`Xhatprime_weighted`}{If \code{postAnalysis==TRUE}, estimated latent positions after removing the effect of covariates associated with `betahat_weighted`.}
#' \item{`clusters_weighted`}{If \code{postAnalysis==TRUE}, block label of each nodes without the effect of covariates associated with `betahat_weighted`.}
#' \item{`iterprime_weighted`}{If \code{postAnalysis==TRUE}, the number of Lanczos iterations carried out after removing the effect of covariates associated with `betahat_weighted`. See \link[irlba]{irlba}.}
#' \item{`mprodprime_weighted`}{If \code{postAnalysis==TRUE}, the total number of matrix vector products carried out after removing the effect of covariates associated with `betahat_weighted`. See \link[irlba]{irlba}.}
#' \item{`pp1`}{Screeplot with covariates.}
#' \item{`pp2`}{Latent position in 2D.}
#' \item{`pp3_simple`}{If \code{postAnalysis==TRUE}, screeplot without covariates associated with `betahat_simple`.}
#' \item{`pp3_weighted`}{If \code{postAnalysis==TRUE}, screeplot without covariates associated with `betahat_weighted`.}
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
#' print(adjustedRandIndex(blocks_cov, result$clusters_cov))          # ARI with covariates
#' print(adjustedRandIndex(blocks, result$clusters_simple))           # ARI without covariates associated with beta_simple
#' print(adjustedRandIndex(blocks, result$clusters_weighted))         # ARI without covariates associated with beta_weighted
#' print(beta)                                       # True beta
#' print(result$betahat_simple)                      # Estimated beta using simple procedure
#' print(result$betahat_simple_unbiased)             # Estimated beta (unbiased) using simple procedure
#' print(result$sd2_simple)                          # Variance of estimated beta using simple procedure
#' print(result$betahat_weighted)                    # Estimated beta using weighted procedure
#' print(result$betahat_weighted_unbiased)           # Estimated beta (unbiased) using weighted procedure
#' print(result$sd2_weighted)                        # Variance of estimated beta using weighted procedure
#' print(B)                                          # True B matrix
#' print(result$BXhat)                               # Estimated B matrix
#'
#' ## Visualization
#' pp2 <- plotLatentPosition(result$Xhat, blocks, withCovariates = TRUE, dhat = ncol(result$Xhat), covariates)
#' pp4_simple <- plotLatentPosition(result$Xhatprime_simple, blocks, withCovariates = FALSE, latent, K, d)
#' multiplot(result$pp1, result$pp3_simple, pp2[[1]], pp4_simple, cols = 2)
#'
#' pp4_weighted <- plotLatentPosition(as.matrix(result$Xhatprime_weighted[,nrow(latent)]), blocks, withCovariates = FALSE, latent, K, d)
#' multiplot(result$pp1, result$pp3_weighted, pp2[[1]], pp4_weighted, cols = 2)
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


GRDPGwithCovariates <- function(A, covariates, link = 'identity', clusterMethod = 'GMM', G = 1:9, dmax = 10, dhat = NULL, maxit = 1000, work = 12, tol = 1e-05, check = 'BF', sd = TRUE, rho = 1, postAnalysis = TRUE, plot = TRUE, ...) {
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

  cat('\n\n', 'Estimating beta...')
  cov <- apply(covariates, 2, function(x) length(unique(x)))
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
  result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link, check, sd, rho)
  betahat1 <- sapply(result1$betahats, mean)
  betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
  sd21 <- sapply(result1$sd2s, mean)

  result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link, check, sd, rho)
  betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
  betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
  sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

  result$BXhat <- BXhat
  result$betahat_simple <- betahat1
  result$betahat_simple_unbiased <- betahat1_unbiased
  result$sd2_simple <- sd21
  result$betahat_weighted <- betahat2
  result$betahat_weighted_unbiased <- betahat2_unbiased
  result$sd2_weighted <- sd22
  result$clusters_cov <- clusters_cov

  if (postAnalysis) {
    cat('\n\n', 'Post Analysis using beta_simple and beta_weighted...')
    if (link == 'logit') {
      Yhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
      Qhat <- Yhat %*% Ipq %*% t(Yhat)
      if (check == 'BF') {
        What <- logit(BFcheck(Qhat))
      } else {
        What <- logit(Removecheck(Qhat))
      }
      Aprime_simple <- getAwithoutCovariates(What, betahat1, covariates)
      Aprime_weighted <- getAwithoutCovariates(What, betahat2, covariates)
    } else {
      Aprime_simple <- getAwithoutCovariates(A, betahat1, covariates)
      Aprime_weighted <- getAwithoutCovariates(A, betahat2, covariates)
    }
    embedprime_simple <- SpectralEmbedding(Aprime_simple, dmax, maxit = maxit, work = work, tol = tol)
    embedprime_weighted <- SpectralEmbedding(Aprime_weighted, dmax, maxit = maxit, work = work, tol = tol)
    sprime_simple <- embedprime_simple$D
    sprime_weighted <- embedprime_weighted$D
    dhatprime_simple <- dimselect(sprime_simple)$elbow[1]
    dhatprime_weighted <- dimselect(sprime_weighted)$elbow[1]
    Xhatprime_simple <- embedprime_simple$X[,1:dhatprime_simple] %*% sqrt(diag(sprime_simple[1:dhatprime_simple], nrow=dhatprime_simple, ncol=dhatprime_simple))
    Xhatprime_weighted <- embedprime_weighted$X[,1:dhatprime_weighted] %*% sqrt(diag(sprime_weighted[1:dhatprime_weighted], nrow=dhatprime_weighted, ncol=dhatprime_weighted))
    if (clusterMethod == 'GMM') {
      model2_simple <- Mclust(Xhatprime_simple, G, verbose = FALSE)
      model2_weighted <- Mclust(Xhatprime_weighted, G, verbose = FALSE)
      clusters_simple <- getClusters(data.frame(model2_simple$z))
      clusters_weighted <- getClusters(data.frame(model2_weighted$z))
    } else {
      model2_simple <- kmeans(Xhatprime_simple, centers)
      model2_weighted <- kmeans(Xhatprime_weighted, centers)
      clusters_simple <- model2_simple$cluster
      clusters_weighted <- model2_weighted$cluster
    }
    result$iterprime_simple <- embedprime_simple$iter
    result$mprodprime_simple <- embedprime_simple$mprod
    result$Xhatprime_simple <- Xhatprime_simple
    result$clusters_simple <- clusters_simple
    result$iterprime_weighted <- embedprime_weighted$iter
    result$mprodprime_weighted <- embedprime_weighted$mprod
    result$Xhatprime_weighted <- Xhatprime_weighted
    result$clusters_weighted <- clusters_weighted
  }

  if (plot) {
    cols <- ncol(A)
    temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
    s1 <- temp1$values
    temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
    s2 <- temp2$values
    tempdat <- data.frame(raw = c(s1,s2)) %>%
      mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
      arrange(desc(s))
    pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (with Covariates)')
    pp2 <- plotLatentPosition(Xhat, withCovariates = TRUE, dhat = dhat, covariates = covariates)
    if (postAnalysis) {
      cols <- ncol(Aprime_simple)
      temp1 <- eigs_sym(matrix(as.numeric(Aprime_simple), ncol = cols), dmax, 'LA')
      s1 <- temp1$values
      temp2 <- eigs_sym(matrix(as.numeric(Aprime_simple), ncol = cols), dmax, 'SA')
      s2 <- temp2$values
      tempdat <- data.frame(raw = c(s1,s2)) %>%
        mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
        arrange(desc(s))
      pp3_simple <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
      result$pp3_simple <- pp3_simple
      cols <- ncol(Aprime_weighted)
      temp1 <- eigs_sym(matrix(as.numeric(Aprime_weighted), ncol = cols), dmax, 'LA')
      s1 <- temp1$values
      temp2 <- eigs_sym(matrix(as.numeric(Aprime_weighted), ncol = cols), dmax, 'SA')
      s2 <- temp2$values
      tempdat <- data.frame(raw = c(s1,s2)) %>%
        mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
        arrange(desc(s))
      pp3_weighted <- scree(sprime_weighted, 'Screeplot (without Covariates)')
      result$pp3_weighted <- pp3_weighted
    }
    multiplot(plotlist = pp2, cols = ceiling(length(pp2)/2))
    result$pp1 <- pp1
    result$pp2 <- pp2
  }


  cat('\n****************************************************************************\n')
  return(result)
}
