#' Estimate beta
#'
#' Estimate beta (the effect of covariates) using weighted procedure.
#'
#' @import mclust
#'
#' @param Xhat Estimated latent positions, should be an `n` by `d` matrix where `n` is the number nodes and `d` is the embeded dimension.
#' @param muhats The center of the latent position for each block. Should be a `d` by `K` matrix where `d` is the embed dimension and `K` is the number of blocks.
#' @param Ipq `Ipq` matrix for GRDPG. See \link{getIpq}.
#' @param cov A vector specifying possible value that each covariate could take. For example, if two binary covariates, then \code{cov <- c(2, 2)}.
#' @param covariates Observed covariates. Should be a `n` by `c` matrix or dataframe where `n` is the number of nodes and `c` is the number of covariates.
#' @param clusters_cov An `n`-vecotr indicating the block label of each nodes with the effect of covariates where `n` is the number nodes.
#' @param link Link function. Could be 'identity' (by default) or 'logit'.
#' @param check Method to check probability matrix. Could be 'BF' (by default, see \link{BFcheck}) or 'Remove' (see \link{Removecheck}).
#' @param sd Whether to compute standard errors of the estimate of beta. \code{TRUE} by default.
#' @param rho Sparsity coefficient. Coule be `1` (by default) or `0`.
#'
#' @return A list containing the following:
#' \describe{
#' \item{`betahats`}{A list containing all estimated beta. The length of the list equals to the number of covariates. Using \code{sapply(Map('*',betahats,pis), sum)} to get the estimate of each beta.}
#' \item{`pis`}{A list containing all weight of each pair. The length of the list equals to the number of covariates.}
#' \item{`bias`}{If \code{sd==TRUE}, A list containing all bias terms for `betahats` in the CLT. The length of the list equals to the number of covariates. Using \code{sapply(Map('*',bias,pis), sum)} to get the bias of each beta and using \code{sapply(Map('*',betahats,pis), sum) - sapply(Map('*',bias,pis), sum)} to get the unbiased estimate of each beta.}
#' \item{`sd2s`}{If \code{sd==TRUE}, A list containing all variances of `betahats`. The length of the list equals to the number of covariates. Using \code{sapply(Map('*',sd2s,pis), sum)} to get the variance of each beta.}
#' \item{`...`}{If \code{sd==TRUE}, Lists containing all `psi`, `sigma2` and `covariances`. The length of the list equals to the number of covariates.}
#' }
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


estimatebeta2 <- function(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'identity', check = 'BF', sd = TRUE, rho = 1) {

  result <- list()

  K <- length(unique(clusters_cov))
  BXhat <- t(muhats) %*% Ipq %*% muhats
  if (link == 'logit') {
    if (check == 'BF') {
      BXhat <- logit(BFcheck(BXhat))
    } else {
      BXhat <- logit(Removecheck(BXhat))
    }
  }
  betahats <- vector('list', ncol(covariates))
  pis <- vector('list', length(cov))
  if (sd) {
    theta <- t(muhats) %*% Ipq %*% muhats
    theta <- BFcheck(theta)
    eta <- getWeight(clusters_cov)$freq
    delta <- 0
    for (kk in 1:length(eta)) {
      delta <- delta + eta[kk] * muhats[,kk,drop=FALSE] %*% t(muhats[,kk,drop=FALSE])
    }
    zeta <- t(muhats) %*% solve(delta) %*% muhats
    bias <- vector('list', ncol(covariates))
    sd2s <- vector('list', ncol(covariates))
    psiij1s <- vector('list', ncol(covariates))
    psiij2s <- vector('list', ncol(covariates))
    sigma2ij1s <- vector('list', ncol(covariates))
    sigma2ij2s <- vector('list', ncol(covariates))
    sigma2j1j1s <- vector('list', ncol(covariates))
    sigma2j2j2s <- vector('list', ncol(covariates))
    covij1ij2s <- vector('list', ncol(covariates))
  }
  model2 <- Mclust(diag(BXhat), ncol(BXhat)/prod(cov), verbose = FALSE)
  c <- getClusters(data.frame(model2$z))
  weight_cov <- getWeight(clusters_cov)
  weight <- data.frame()
  for (clusters in unique(c)) {
    ind <- which(c==clusters)
    count <- sum(weight_cov$count[weight_cov$clusters %in% ind])
    freq <- sum(weight_cov$freq[weight_cov$clusters %in% ind])
    weight <- rbind(weight, data.frame(clusters,count,freq))
  }
  weights_cov <- vector('list', ncol(covariates))
  for (k in 1:ncol(covariates)) {
    temp <- data.frame(clusters_cov, covariates[,k])
    names(temp)[2] <- 'covariates'
    weights_cov[[k]] <- temp %>%
      group_by(clusters_cov, covariates) %>%
      summarise(count = n()) %>%
      left_join(weight_cov, by = c('clusters_cov'='clusters')) %>%
      mutate(freq = count.x/count.y)
  }
  for (k in 1:ncol(covariates)) {
    for (i in 1:nrow(BXhat)) {
      for (j1 in 1:(ncol(BXhat)-1)) {
        for (j2 in (j1+1):ncol(BXhat)) {
          if (c[j1] == c[j2]) {
            t1 <- filter(weights_cov[[k]], clusters_cov == j1)
            t2 <- filter(weights_cov[[k]], clusters_cov == j2)
            for (l1 in 1:nrow(t1)) {
              ind1 <- setdiff(unique(t2$covariates), t1$covariates[l1])
              if (length(ind1) > 0) {
                for (l2 in ind1) {
                  # pi_1 <- weight$freq[weight$clusters==c[i]]
                  # pi_2 <- weight$freq[weight$clusters==c[j1]]
                  pi_cov_1 <- weight_cov$freq[weight_cov$clusters==j1]
                  pi_cov_2 <- weight_cov$freq[weight_cov$clusters==j2]
                  pi_z_1 <- t1$freq[l1]
                  pi_z_2 <- t2$freq[t2$covariates==l2]
                  temp <- setdiff(1:ncol(covariates), k)
                  pi_0 <- 1
                  if (length(temp) > 0) {
                    for (l in temp) {
                      temp_pi_0 <- 0
                      t3 <- filter(weights_cov[[l]], clusters_cov == j1)
                      t4 <- filter(weights_cov[[l]], clusters_cov == j2)
                      for (t_cov in unique(t3$covariates)) {
                        if (t_cov %in% t4$covariates) {
                          temp_pi_0 <- temp_pi_0 + sum(t3$freq[t3$covariates==t_cov]) * sum(t4$freq[t4$covariates==t_cov])
                        } else {
                          temp_pi_0 <- temp_pi_0 + 0
                        }
                      }
                      pi_0 <- c(pi_0, temp_pi_0)
                    }
                  }
                  betahats[[k]] <- c(betahats[[k]], BXhat[i,j1] - BXhat[i,j2])
                  pis[[k]] <- c(pis[[k]], prod(pi_0) * pi_cov_1 * pi_cov_2 * pi_z_1 * pi_z_2)
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
                    cov1 <- mean(E[[i]])*t(muhats[,j1,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,j2,drop=FALSE] / eta[i]
                    cov2 <- ifelse(i == j2, mean(E[[i]])*t(muhats[,j1,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE] / eta[i], 0)
                    cov3 <- ifelse(i == j1, mean(E[[j1]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,j2,drop=FALSE] / eta[j1], 0)
                    cov4 <- ifelse(j1 == j2, mean(E[[j1]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE] / eta[j1], 0)
                    cov5 <- mean(E[[i]])*t(muhats[,j1,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]%*%t(muhats[,j2,drop=FALSE])%*%solve(delta)%*%muhats[,i,drop=FALSE]
                    cov6 <- mean(E[[j1]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]%*%t(muhats[,j2,drop=FALSE])%*%solve(delta)%*%muhats[,j1,drop=FALSE]
                    cov7 <- mean(E[[i]])*t(muhats[,i,drop=FALSE])%*%solve(delta)%*%muhats[,j1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,j2,drop=FALSE]
                    cov8 <- mean(E[[j2]])*t(muhats[,j2,drop=FALSE])%*%solve(delta)%*%muhats[,j1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]
                    cov9 <- 0
                    psiij11 <- psiij21 <- 0
                    sigma2j1j12 <- sigma2j2j22 <- 0
                    sigma2ij13 <- sigma2ij14 <- sigma2ij15 <- 0
                    sigma2ij23 <- sigma2ij24 <- sigma2ij25 <- 0
                    for (r in 1:K) {
                      psiij11 <- psiij11 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[j1,r]*(1-theta[j1,r]))*(t(muhats[,i,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%muhats[,j1,drop=FALSE])
                      psiij21 <- psiij21 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[j2,r]*(1-theta[j2,r]))*(t(muhats[,i,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%muhats[,j2,drop=FALSE])
                      sigma2j1j12 <- sigma2j1j12 + eta[r]*theta[j1,r]*(1-theta[j1,r])*zeta[j1,r]^2*(1/eta[j1]-2*zeta[j1,j1])
                      sigma2j2j22 <- sigma2j2j22 + eta[r]*theta[j2,r]*(1-theta[j2,r])*zeta[j2,r]^2*(1/eta[j2]-2*zeta[j2,j2])
                      sigma2ij13 <- sigma2ij13 + eta[r]*theta[i,r]*(1-theta[i,r])*zeta[j1,r]^2*(1/eta[i]-2*zeta[i,i])
                      sigma2ij14 <- sigma2ij14 + eta[r]*theta[j1,r]*(1-theta[j1,r])*zeta[i,r]^2*(1/eta[j1]-2*zeta[j1,j1])
                      sigma2ij15 <- sigma2ij15 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[j1,r]*(1-theta[j1,r]))*zeta[i,r]*zeta[r,j1]*zeta[i,j1]
                      sigma2ij23 <- sigma2ij23 + eta[r]*theta[i,r]*(1-theta[i,r])*zeta[j2,r]^2*(1/eta[i]-2*zeta[i,i])
                      sigma2ij24 <- sigma2ij24 + eta[r]*theta[j2,r]*(1-theta[j2,r])*zeta[i,r]^2*(1/eta[j2]-2*zeta[j2,j2])
                      sigma2ij25 <- sigma2ij25 + eta[r]*(theta[i,r]*(1-theta[i,r])+theta[j2,r]*(1-theta[j2,r]))*zeta[i,r]*zeta[r,j2]*zeta[i,j2]
                      cov9 <- cov9 + eta[r]*mean(E[[r]])*t(muhats[,r,drop=FALSE])%*%solve(delta)%*%muhats[,j1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])%*%solve(delta)%*%solve(delta)%*%muhats[,i,drop=FALSE]%*%t(muhats[,j2,drop=FALSE])%*%solve(delta)%*%muhats[,r,drop=FALSE]
                    }
                    psiij12 <- psiij22 <- 0
                    sigma2j1j13 <- sigma2j2j23 <- 0
                    sigma2ij16 <- sigma2ij26 <- 0
                    for (r in 1:K) {
                      for (s in 1:K) {
                        psiij12 <- psiij12 + eta[r]*eta[s]*theta[s,r]*(1-theta[s,r])*(t(muhats[,s,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%(muhats[,j1,drop=FALSE]%*%t(muhats[,i,drop=FALSE])+muhats[,i,drop=FALSE]%*%t(muhats[,j1,drop=FALSE]))%*%solve(delta)%*%muhats[,s,drop=FALSE])
                        psiij22 <- psiij22 + eta[r]*eta[s]*theta[s,r]*(1-theta[s,r])*(t(muhats[,s,drop=FALSE])%*%solve(delta)%*%Ipq%*%solve(delta)%*%(muhats[,j2,drop=FALSE]%*%t(muhats[,i,drop=FALSE])+muhats[,i,drop=FALSE]%*%t(muhats[,j2,drop=FALSE]))%*%solve(delta)%*%muhats[,s,drop=FALSE])
                        sigma2j1j13 <- sigma2j1j13 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*zeta[j1,r]^2*zeta[j1,s]^2
                        sigma2j2j23 <- sigma2j2j23 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*zeta[j2,r]^2*zeta[j2,s]^2
                        sigma2ij16 <- sigma2ij16 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*(zeta[i,r]*zeta[j1,s]+zeta[j1,r]*zeta[i,s])^2
                        sigma2ij26 <- sigma2ij26 + eta[r]*eta[s]*theta[r,s]*(1-theta[r,s])*(zeta[i,r]*zeta[j2,s]+zeta[j2,r]*zeta[i,s])^2
                      }
                    }
                    psiij1 <- psiij11 - psiij12
                    psiij2 <- psiij21 - psiij22
                    sigma2ij1 <- (theta[i,i]*(1-theta[i,i])+theta[j1,j1]*(1-theta[j1,j1]))*zeta[i,j1]^2+2*theta[i,j1]*(1-theta[i,j1])*zeta[i,i]*zeta[j1,j1]+sigma2ij13+sigma2ij14-2*sigma2ij15+sigma2ij16/2
                    sigma2ij2 <- (theta[i,i]*(1-theta[i,i])+theta[j2,j2]*(1-theta[j2,j2]))*zeta[i,j2]^2+2*theta[i,j2]*(1-theta[i,j2])*zeta[i,i]*zeta[j2,j2]+sigma2ij23+sigma2ij24-2*sigma2ij25+sigma2ij26/2
                    sigma2j1j1 <- 4*theta[j1,j1]*(1-theta[j1,j1])*zeta[j1,j1]^2+4*sigma2j1j12+2*sigma2j1j13
                    sigma2j2j2 <- 4*theta[j2,j2]*(1-theta[j2,j2])*zeta[j2,j2]^2+4*sigma2j2j22+2*sigma2j2j23
                    covij1ij2 <- cov1 + cov2 + cov3 + cov4 + cov5 + cov6 + cov7 + cov8 + cov9
                    if (link == 'logit') {
                      psiij1 <- psiij1 * logitderivative(BFcheck(theta[i,j1]))
                      psiij2 <- psiij2 * logitderivative(BFcheck(theta[i,j2]))
                      sigma2ij1 <- sigma2ij1 * logitderivative(BFcheck(theta[i,j1]))^2
                      sigma2ij2 <- sigma2ij2 * logitderivative(BFcheck(theta[i,j2]))^2
                      sigma2j1j1 <- sigma2j1j1 * logitderivative(BFcheck(theta[j1,j1]))^2
                      sigma2j2j2 <- sigma2j2j2 * logitderivative(BFcheck(theta[j2,j2]))^2
                      covij1ij2 <- covij1ij2 * logitderivative(BFcheck(theta[i,j1])) * logitderivative(BFcheck(theta[i,j2]))
                    }
                    psiij1s[[k]] <- c(psiij1s[[k]], psiij1)
                    psiij2s[[k]] <- c(psiij1s[[k]], psiij2)
                    sigma2ij1s[[k]] <- c(sigma2ij1s[[k]], sigma2ij1)
                    sigma2ij2s[[k]] <- c(sigma2ij2s[[k]], sigma2ij2)
                    sigma2j1j1s[[k]] <- c(sigma2j1j1s[[k]], sigma2j1j1)
                    sigma2j2j2s[[k]] <- c(sigma2j2j2s[[k]], sigma2j2j2)
                    covij1ij2s[[k]] <- c(covij1ij2s[[k]], covij1ij2)
                    if (i == j1) {
                      tempsd <- sigma2j1j1 + sigma2ij2 - 2 * covij1ij2
                    } else if (i == j2) {
                      tempsd <- sigma2ij1 + sigma2j2j2 - 2 * covij1ij2
                    } else {
                      tempsd <- sigma2ij1 + sigma2ij2 - 2 * covij1ij2
                    }
                    tempbias <- (psiij1 - psiij2) / length(clusters_cov)
                    bias[[k]] <- c(bias[[k]], tempbias)
                    sd2s[[k]] <- c(sd2s[[k]], tempsd)
                  }
                }
              }
            }
          }
        }
      }
    }
    symbol <- ifelse(BXhat[1,1]-BXhat[1,which(c==c[1])[2]]>0, 1, -1)
    pis[[k]] <- pis[[k]] / sum(pis[[k]])
    betahats[[k]] <- abs(betahats[[k]]) * symbol
  }
  result$betahats <- betahats
  result$pis <- pis
  if (sd) {
    result$bias <- bias
    result$sd2s <- sd2s
    result$psiij1s <- psiij1s
    result$psiij2s <- psiij2s
    result$sigma2ij1s <- sigma2ij1s
    result$sigma2ij2s <- sigma2ij2s
    result$sigma2j1j1s <- sigma2j1j1s
    result$sigma2j2j2s <- sigma2j2j2s
    result$covij1ij2s <- covij1ij2s
  }

  return(result)

}














