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
#' \item{`sd2s`}{If \code{sd==TRUE}, A list containing all variances of `betahats`. The length of the list equals to the number of covariates. Using \code{sapply(sd2s, mean)} to get the variance of each beta.}
#' }
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


estimatebeta <- function(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'identity', check = 'BF', sd = TRUE, rho = 1) {

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
  theta <- t(muhats) %*% Ipq %*% muhats
  eta <- getWeight(clusters_cov)$freq
  delta <- 0
  for (kk in 1:length(eta)) {
    delta <- delta + eta[kk] * muhats[,kk,drop=FALSE] %*% t(muhats[,kk,drop=FALSE])
  }
  zeta <- t(muhats) %*% solve(delta) %*% muhats
  betahats <- vector('list', ncol(covariates_block))
  if (sd) {
    bias <- vector('list', ncol(covariates_block))
    sd2s <- vector('list', ncol(covariates_block))
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
                cov1 <- mean(E[[i]])*t(muhats[,l1,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,l2,drop=FALSE] / eta[i]
                cov2 <- ifelse(i == l2, mean(E[[i]])*t(muhats[,l1,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE] / eta[i], 0)
                cov3 <- ifelse(i == l1, mean(E[[l1]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,l2,drop=FALSE] / eta[l1], 0)
                cov4 <- ifelse(l1 == l2, mean(E[[l1]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE] / eta[l1], 0)
                cov5 <- mean(E[[i]])*t(muhats[,l1,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]%*%t(muhats[,l2,drop=FALSE])%*%solve(delta)%*%muhats[,i,drop=FALSE]
                cov6 <- mean(E[[l1]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]%*%t(muhats[,l2,drop=FALSE])%*%solve(delta)%*%muhats[,l1,drop=FALSE]
                cov7 <- mean(E[[i]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%muhats[,l1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,l2,drop=FALSE]
                cov8 <- mean(E[[l2]])*t(muhats[,l2,drop=FALSE])%*%solve(delta)%*%muhats[,l1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]
                cov9 <- 0
                psiil11 <- psiil21 <- 0
                sigma2l1l12 <- sigma2l2l22 <- 0
                sigma2il13 <- sigma2il14 <- sigma2il15 <- 0
                sigma2il23 <- sigma2il24 <- sigma2il25 <- 0
                for (r in 1:K) {
                  psiil11 <- psiil11 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[l1,r]*(1-theta[l1,r]))*(t(muhats[,i,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%muhats[,l1,drop=FALSE])
                  psiil21 <- psiil21 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[l2,r]*(1-theta[l2,r]))*(t(muhats[,i,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%muhats[,l2,drop=FALSE])
                  sigma2l1l12 <- sigma2l1l12 + eta[r]*theta[l1,r]*(1-theta[l1,r])*zeta[l1,r]^2*(1/eta[l1]-2*zeta[l1,l1])
                  sigma2l2l22 <- sigma2l2l22 + eta[r]*theta[l2,r]*(1-theta[l2,r])*zeta[l2,r]^2*(1/eta[l2]-2*zeta[l2,l2])
                  sigma2il13 <- sigma2il13 + eta[r]*theta[i,r]*(1-theta[i,r])*zeta[l1,r]^2*(1/eta[i]-2*zeta[i,i])
                  sigma2il14 <- sigma2il14 + eta[r]*theta[l1,r]*(1-theta[l1,r])*zeta[i,r]^2*(1/eta[l1]-2*zeta[l1,l1])
                  sigma2il15 <- sigma2il15 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[l1,r]*(1-theta[l1,r]))*zeta[i,r]*zeta[r,l1]*zeta[i,l1]
                  sigma2il23 <- sigma2il23 + eta[r]*theta[i,r]*(1-theta[i,r])*zeta[l2,r]^2*(1/eta[i]-2*zeta[i,i])
                  sigma2il24 <- sigma2il24 + eta[r]*theta[l2,r]*(1-theta[l2,r])*zeta[i,r]^2*(1/eta[l2]-2*zeta[l2,l2])
                  sigma2il25 <- sigma2il25 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[l2,r]*(1-theta[l2,r]))*zeta[i,r]*zeta[r,l2]*zeta[i,l2]
                  cov9 <- cov9 + eta[r]*mean(E[[r]])*t(muhats[,r,drop=FALSE])%*%solve(delta)%*%muhats[,l1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]%*%t(muhats[,l2,drop=FALSE])%*%solve(delta)%*%muhats[,r,drop=FALSE]
                }
                psiil12 <- psiil22 <- 0
                sigma2l1l13 <- sigma2l2l23 <- 0
                sigma2il16 <- sigma2il26 <- 0
                for (r in 1:K) {
                  for (s in 1:K) {
                    psiil12 <- psiil12 + eta[r]*eta[s]*theta[s,r]*(1-theta[s,r])*(t(muhats[,s,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%(muhats[,l1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])+muhats[,i,drop=FALSE]%*%t(muhats[,l1,drop=FALSE]))%*%solve(delta)%*%muhats[,s,drop=FALSE])
                    psiil22 <- psiil22 + eta[r]*eta[s]*theta[s,r]*(1-theta[s,r])*(t(muhats[,s,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%(muhats[,l2,drop=FALSE]%*%t(muhats[,i,drop=FALSE])+muhats[,i,drop=FALSE]%*%t(muhats[,l2,drop=FALSE]))%*%solve(delta)%*%muhats[,s,drop=FALSE])
                    sigma2l1l13 <- sigma2l1l13 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*zeta[l1,r]^2*zeta[l1,s]^2
                    sigma2l2l23 <- sigma2l2l23 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*zeta[l2,r]^2*zeta[l2,s]^2
                    sigma2il16 <- sigma2il16 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*(zeta[i,r]*zeta[l1,s]+zeta[l1,r]*zeta[i,s])^2
                    sigma2il26 <- sigma2il26 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*(zeta[i,r]*zeta[l2,s]+zeta[l2,r]*zeta[i,s])^2
                  }
                }
                psiil1 <- psiil11 - psiil12
                psiil2 <- psiil21 - psiil22
                sigma2il1 <- (theta[i,i]*(1-theta[i,i])+theta[l1,l1]*(1-theta[l1,l1]))*zeta[i,l1]^2+2*theta[i,l1]*(1-theta[i,l1])*zeta[i,i]*zeta[l1,l1]+sigma2il13+sigma2il14-2*sigma2il15+sigma2il16/2
                sigma2il2 <- (theta[i,i]*(1-theta[i,i])+theta[l2,l2]*(1-theta[l2,l2]))*zeta[i,l2]^2+2*theta[i,l2]*(1-theta[i,l2])*zeta[i,i]*zeta[l2,l2]+sigma2il23+sigma2il24-2*sigma2il25+sigma2il26/2
                sigma2l1l1 <- 4*theta[l1,l1]*(1-theta[l1,l1])*zeta[l1,l1]^2+4*sigma2l1l12+2*sigma2l1l13
                sigma2l2l2 <- 4*theta[l2,l2]*(1-theta[l2,l2])*zeta[l2,l2]^2+4*sigma2l2l22+2*sigma2l2l23
                covil1il2 <- cov1 + cov2 + cov3 + cov4 + cov5 + cov6 + cov7 + cov8 + cov9
                if (link == 'logit') {
                  psiil1 <- psiil1 * logitderivative(BFcheck(theta[i,l1]))
                  psiil2 <- psiil2 * logitderivative(BFcheck(theta[i,l2]))
                  sigma2il1 <- sigma2il1 * logitderivative(BFcheck(theta[i,l1]))^2
                  sigma2il2 <- sigma2il2 * logitderivative(BFcheck(theta[i,l2]))^2
                  sigma2l1l1 <- sigma2l1l1 * logitderivative(BFcheck(theta[l1,l1]))^2
                  sigma2l2l2 <- sigma2l2l2 * logitderivative(BFcheck(theta[l2,l2]))^2
                  covil1il2 <- covil1il2 * logitderivative(BFcheck(theta[i,l1])) * logitderivative(BFcheck(theta[i,l2]))
                }
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
  }

  return(result)

}







