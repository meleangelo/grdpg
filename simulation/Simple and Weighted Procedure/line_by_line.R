require(dplyr)
require(mclust)
require(graphstats)
require(grdpg)

## Some hyperparameter
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 12               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 

## Embedding (giving adjacency matrix A)
embedding <- embed(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- ifelse(is.null(dhat), dimselect(s)$elbow+1, dhat)

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
if (dhat == 1) {
  Ipq <- matrix(1)
} else {
  temp <- eigen(A)
  Ipq <- getIpq(temp$values, dhat)
}

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats

## Calculate possible value the covariate could take (e.x. cov = 2 for one binary covariate, giving covariate)
cov <- apply(covariates, 2, function(x) length(unique(x)))

## Create a list to save all betahats
betahats <- vector('list', length(cov))

## Cluster the diagonal of the estimated probability matrix
model2 <- Mclust(diag(BXhat), ncol(BXhat)/prod(cov), verbose = FALSE)

## Get the block label if we remove the effect of covariate (ideally from 2K block to K block)
c <- getClusters(data.frame(model2$z))

## Estimate the size of each block using frequency
weight_cov <- getWeight(clusters_cov)
weight <- data.frame()
for (clusters in unique(c)) {
  ind <- which(c==clusters)
  count <- sum(weight_cov$count[weight_cov$clusters %in% ind])
  freq <- sum(weight_cov$freq[weight_cov$clusters %in% ind])
  weight <- rbind(weight, data.frame(clusters,count,freq))
}

## Estimate the porportion of each covariate in different block using frequency
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

## Loop over all the possible pair in estimated probability matrix
for (k in 1:ncol(covariates)) {
  for (i in 1:nrow(BXhat)) {
    for (j1 in 1:ncol(BXhat)) {
      for (j2 in 1:ncol(BXhat)) {
        # Only consider pair that should belong to the same block without covariate
        if (c[j1] == c[j2]) {
          t1 <- filter(weights_cov[[k]], clusters_cov == j1)
          t2 <- filter(weights_cov[[k]], clusters_cov == j2)
          for (l1 in 1:nrow(t1)) {
            ind1 <- setdiff(unique(t2$covariates), t1$covariates[l1])
            if (length(ind1) > 0) {
              for (l2 in ind1) {
                # Estimated size of each block for this pair
                pi_1 <- weight$freq[weight$clusters==c[i]]
                pi_2 <- weight$freq[weight$clusters==c[j1]]
                # Estimated proportion of each covariate for this pair
                pi_z_1 <- t1$freq[t1$covariates==l1]
                pi_z_2 <- t2$freq[t2$covariates==l2]
                # Weighted by the estimated size of the block and proportion of each covariate
                betahats[[k]] <- c(betahats[[k]], pi_1 * pi_2 * pi_z_1 * pi_z_2 * (BXhat[i,j1] - BXhat[i,j2]) / cov[k])
              }
            }
          }
        }
      }
    }
  }
  symbol <- ifelse(BXhat[1,1]-BXhat[1,which(c==c[1])[2]]>0, 1, -1)
  betahats[[k]] <- abs(betahats[[k]]) * symbol
}

## Estimated betahat
betahat <- sapply(betahats, sum)






