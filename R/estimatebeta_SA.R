#' Estimate beta
#'
#' Estimate beta (the effect of covariates) using simple procedure.
#'
#' @import mclust
#'
#' @param Xhat Estimated latent positions, should be an `n` by `d` matrix where `n` is the number nodes and `d` is the embeded dimension.
#' @param muhats The center of the latent position for each block. Should be a `d` by `K` matrix where `d` is the embed dimension and `K` is the number of blocks.
#' @param Ipq `Ipq` matrix for GRDPG. See \link{getIpq}.
#' @param cov A vector specifying possible value that each covariate could take. For example, if two binary covariates, then \code{cov <- c(2, 2)}.
#' @param covariates_block Estimated covariates for each block. Should be a `k` by `c` matrix or dataframe where `k` is the number of blocks and `c` is the number of covariates.
#' @param clusters_cov An `n`-vecotr indicating the block label of each nodes with the effect of covariates where `n` is the number nodes.
#' @param link Link function. Could be 'identity' (by default) or 'logit'.
#' @param check Method to check probability matrix. Could be 'BF' (by default, see \link{BFcheck}) or 'Remove' (see \link{Removecheck}).
#' @param sd Whether to compute standard errors of the estimate of beta. \code{TRUE} by default.
#' @param rho Sparsity coefficient. Coule be `1` (by default) or `0`.
#'
#' @return A list containing the following:
#' \describe{
#' \item{`betahats`}{A list containing all estimated beta. The length of the list equals to the number of covariates. Using \code{sapply(betahats, mean)} to get the estimate of each beta.}
#' \item{`bias`}{If \code{sd==TRUE}, A list containing all bias terms for `betahats` in the CLT. The length of the list equals to the number of covariates. Using \code{sapply(bias, mean)} to get the bias of each beta and using \code{sapply(betahats, mean) - sapply(bias, mean)} to get the unbiased estimate of each beta.}
#' \item{`sd2s`}{If \code{sd==TRUE}, A list containing all variances of `betahats`. The length of the list equals to the number of covariates. Using \code{sapply(sd2s, mean)/n^2} to get the variance of each beta where `n` is the number of nodes.}
#' \item{`...`}{If \code{sd==TRUE}, Lists containing all `psi`, `sigma2` and `covariances`. The length of the list equals to the number of covariates.}
#' }
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


estimatebeta_SA <- function(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'identity', check = 'BF', sd = TRUE, rho = 1) {

  result <- list()

  covariates_block <- data.frame(covariates_block)
  if (length(cov) != ncol(covariates_block)) {
    stop("The length of `cov` should equal to the number of columns in `covariates_block` (both equal to the number of covariates).")
  }
  if (!(link %in% c('identity', 'logit'))) {
    print("Unrecognized `link`, would use 'identity' by default.")
  }
  if (!(check %in% c('BF', 'Remove'))) {
    print("Unrecognized `check`, would use 'BF' by default.")
  }

  K <- length(unique(clusters_cov))
  BXhat <- t(muhats) %*% Ipq %*% muhats
  if (link == 'logit') {
    if (check == 'BF') {
      BXhat <- logit(BFcheck(BXhat))
    } else {
      BXhat <- logit(Removecheck(BXhat))
    }
  }
  betahats <- vector('list', ncol(covariates_block))
  if (sd) {
    theta <- t(muhats) %*% Ipq %*% muhats
    theta <- BFcheck(theta)
    eta <- getWeight(clusters_cov)$freq
    delta <- 0
    for (kk in 1:length(eta)) {
      delta <- delta + eta[kk] * muhats[,kk,drop=FALSE] %*% t(muhats[,kk,drop=FALSE])
    }
    deltainv <- solve(delta)
    zeta <- t(muhats) %*% deltainv %*% muhats
    E <- vector('list', length(unique(clusters_cov)))
    for (alpha in unique(clusters_cov)) {
      for (n in 1:nrow(Xhat)) {
        if (rho == 0) {
          E[[alpha]] <- c(E[[alpha]], theta[alpha,clusters_cov[n]]*(Xhat[n,,drop=FALSE]%*%t(Xhat[n,,drop=FALSE])))
        } else {
          E[[alpha]] <- c(E[[alpha]], theta[alpha,clusters_cov[n]]*(1-theta[alpha,clusters_cov[n]])*(Xhat[n,,drop=FALSE]%*%t(Xhat[n,,drop=FALSE])))
        }
      }
    }
    bias <- vector('list', ncol(covariates_block))
    sd2s <- vector('list', ncol(covariates_block))
    psiil1s <- vector('list', ncol(covariates_block))
    psiil2s <- vector('list', ncol(covariates_block))
    sigma2il1s <- vector('list', ncol(covariates_block))
    sigma2il2s <- vector('list', ncol(covariates_block))
    sigma2l1l1s <- vector('list', ncol(covariates_block))
    sigma2l2l2s <- vector('list', ncol(covariates_block))
    covil1il2s <- vector('list', ncol(covariates_block))
  }
  model2 <- Mclust(diag(BXhat), ncol(BXhat)/prod(cov), verbose = FALSE)
  c <- getClusters(data.frame(model2$z))
  for (i in 1:nrow(BXhat)) {
    for (k in 1:ncol(covariates_block)) {
      ind1 <- which(covariates_block[,k]==covariates_block[i,k])
      ind2 <- which(covariates_block[,k]!=covariates_block[i,k])
      for (l1 in ind1) {
        for (l2 in ind2) {
          if (c[l1] == c[l2]) {
            temp <- setdiff(1:ncol(covariates_block), k)
            ind <- c()
            if (length(temp) > 0) {
              for (l in temp) {
                ind <- c(ind, covariates_block[l1,l] == covariates_block[l2,l])
              }
            } else {
              ind <- TRUE
            }
            if (all(ind)) {
              betahats[[k]] <- c(betahats[[k]], BXhat[i,l1] - BXhat[i,l2])
              if (sd) {
                meanE <- sapply(E, mean)
                covs_temps <- compute_cov(i, l1, l2, K, meanE, muhats, deltainv, eta, theta, zeta, Ipq)
                covs <- covs_temps[[1]]
                temps <- covs_temps[[2]]
                psiil1 <- temps[1]
                psiil2 <- temps[2]
                sigma2il1 <- temps[3]
                sigma2il2 <- temps[4]
                sigma2l1l1 <- temps[5]
                sigma2l2l2 <- temps[6]
                covil1il2 <- covs[1] + covs[2] + covs[3] + covs[4] - covs[5] - covs[6] - covs[7] - covs[8] + covs[9]
                if (link == 'logit') {
                  psiil1 <- psiil1 * logitderivative(BFcheck(theta[i,l1]))
                  psiil2 <- psiil2 * logitderivative(BFcheck(theta[i,l2]))
                  sigma2il1 <- sigma2il1 * logitderivative(BFcheck(theta[i,l1]))^2
                  sigma2il2 <- sigma2il2 * logitderivative(BFcheck(theta[i,l2]))^2
                  sigma2l1l1 <- sigma2l1l1 * logitderivative(BFcheck(theta[l1,l1]))^2
                  sigma2l2l2 <- sigma2l2l2 * logitderivative(BFcheck(theta[l2,l2]))^2
                  covil1il2 <- covil1il2 * logitderivative(BFcheck(theta[i,l1])) * logitderivative(BFcheck(theta[i,l2]))
                }
                psiil1s[[k]] <- c(psiil1s[[k]], psiil1)
                psiil2s[[k]] <- c(psiil1s[[k]], psiil2)
                sigma2il1s[[k]] <- c(sigma2il1s[[k]], sigma2il1)
                sigma2il2s[[k]] <- c(sigma2il2s[[k]], sigma2il2)
                sigma2l1l1s[[k]] <- c(sigma2l1l1s[[k]], sigma2l1l1)
                sigma2l2l2s[[k]] <- c(sigma2l2l2s[[k]], sigma2l2l2)
                covil1il2s[[k]] <- c(covil1il2s[[k]], covil1il2)
                if (i == l1) {
                  tempsd <- sigma2l1l1 + sigma2il2 - 2 * covil1il2
                } else if (i == l2) {
                  tempsd <- sigma2il1 + sigma2l2l2 - 2 * covil1il2
                } else {
                  tempsd <- sigma2il1 + sigma2il2 - 2 * covil1il2
                }
                tempbias <- (psiil1 - psiil2) / length(clusters_cov)
                bias[[k]] <- c(bias[[k]], tempbias)
                sd2s[[k]] <- c(sd2s[[k]], tempsd)
              }
            }
          }
        }
      }
    }
  }
  result$betahats <- betahats
  if (sd) {
    result$bias <- bias
    result$sd2s <- sd2s
    result$psiil1s <- psiil1s
    result$psiil2s <- psiil2s
    result$sigma2il1s <- sigma2il1s
    result$sigma2il2s <- sigma2il2s
    result$sigma2l1l1s <- sigma2l1l1s
    result$sigma2l2l2s <- sigma2l2l2s
    result$covil1il2s <- covil1il2s
  }

  return(result)

}





