require(dplyr)
require(mclust)
require(graphstats)
require(grdpg)



## Estimate weight using frequency
getWeight <- function(clusters) {
  n <- length(clusters)
  weight <- data.frame(clusters) %>%
    group_by(clusters) %>%
    summarise(count = n()) %>%
    mutate(freq = count/n)
  return(weight)
}


## Estimate beta using weighted procedure
estimatebeta2 <- function(BXhat, cov, covariates, clusters_cov) {
  betahats <- vector('list', length(cov))
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
      for (j1 in 1:ncol(BXhat)) {
        for (j2 in 1:ncol(BXhat)) {
          if (c[j1] == c[j2]) {
            t1 <- filter(weights_cov[[k]], clusters_cov == j1)
            t2 <- filter(weights_cov[[k]], clusters_cov == j2)
            for (l1 in 1:nrow(t1)) {
              ind1 <- setdiff(unique(t2$covariates), t1$covariates[l1])
              if (length(ind1) > 0) {
                for (l2 in ind1) {
                  pi_1 <- weight$freq[weight$clusters==c[i]]
                  pi_2 <- weight$freq[weight$clusters==c[j1]]
                  pi_z_1 <- t1$freq[t1$covariates==l1]
                  pi_z_2 <- t2$freq[t2$covariates==l2]
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
  return(betahats)
}



## Hyperparameter
G <- 1:9 
dmax <- 10
dhat <- NULL
maxit <- 1000
work <- 12
tol <- 1e-05


addCovariates <- TRUE
seed <- 2019

n <- 2000
K <- 2
d <- 1
p <- 0.1
q <- 0.7
latent <- cbind(p, q)


## Change the size of each block
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Change the proportion of each covariate
beta <- 0.3
cov <- 2
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}


## Generate network
B <- generateB(latent, K, d, addCovariates, cov, beta)
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)

# allresults <- data.frame()

## Simulation
for (i in 101:105) {
  cat(i, '\n')
  
  A <- generateA(n, P, seed = i)
  
  embedding <- embed(A, dmax, maxit = maxit, work = work, tol = tol)
  s <- embedding$D
  dhat <- ifelse(is.null(dhat), dimselect(s)$elbow+1, dhat)
  
  temp <- eigen(A)
  Ipq <- getIpq(temp$values, dhat)
  
  Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
  
  model <- Mclust(Xhat, G, verbose = FALSE)
  muhats <- matrix(model$parameters$mean, nrow = dhat)
  BXhat <- t(muhats) %*% Ipq %*% muhats
  
  clusters_cov <- getClusters(data.frame(model$z))
  covariates_block <- getBlockCovariates(covariates, clusters_cov)
  
  betahats1 <- estimatebeta(BXhat, cov, covariates_block)
  betahat1 <- sapply(betahats1, mean)
  
  betahats2 <- estimatebeta2(BXhat, cov, covariates, clusters_cov)
  betahat2 <- sapply(betahats2, sum)
  
  allresults <- rbind(allresults, data.frame(n, pi_1, pi_2, pi_z_11, pi_z_12, pi_z_21, pi_z_22, beta, betahat1, betahat2))
}





