#' Estimate beta
#'
#' Estimate beta (the effect of covariates) using weighted procedure.
#'
#' @import mclust
#'
#' @param BXhat Estimated block probability matrix. Should be a `k` by `k` matrix where `k` is the number of blocks.
#' @param cov A vector specifying possible value that each covariate could take. For example, if two binary covariates, then \code{cov <- c(2, 2)}.
#' @param covariates Observed covariates. Should be a `n` by `c` matrix or dataframe where `n` is the number of nodes and `c` is the number of covariates.
#' @param clusters_cov An `n`-vecotr indicating the block label of each nodes with the effect of covariates where `n` is the number nodes.
#'
#' @return A list containing all estimated beta. The length of the list equals to the number of covariates. Using \code{sapply(betahats, sum)} to get the estimate of each beta.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


estimatebeta2 <- function(BXhat, cov, covariates, clusters_cov) {
  betahats <- vector('list', length(cov))
  pis <- vector('list', length(cov))
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
                  pi_0 <- 0
                  if (length(temp) > 1) {
                    for (l in temp) {
                      t3 <- filter(weights_cov[[l]], clusters_cov == j1)
                      t4 <- filter(weights_cov[[l]], clusters_cov == j2)
                      for (t_cov in unique(t3$covariates)) {
                        if (t_cov %in% t4$covariates) {
                          pi_0 <- pi_0 + sum(t3$freq[t3$covariates==t_cov]) * sum(t4$freq[t4$covariates==t_cov])
                        } else {
                          pi_0 <- pi_0 + 0
                        }
                      }
                    }
                  } else {
                    pi_0 <- 1
                  }
                  betahats[[k]] <- c(betahats[[k]], BXhat[i,j1] - BXhat[i,j2])
                  pis[[k]] <- c(pis[[k]], pi_0 * pi_cov_1 * pi_cov_2 * pi_z_1 * pi_z_2)
                }
              }
            }
          }
        }
      }
    }
    symbol <- ifelse(BXhat[1,1]-BXhat[1,which(c==c[1])[2]]>0, 1, -1)
    pis[[k]] <- pis[[k]] / sum(pis[[k]])
    betahats[[k]] <- abs(betahats[[k]]) * symbol * pis[[k]]
  }
  return(betahats)
}




