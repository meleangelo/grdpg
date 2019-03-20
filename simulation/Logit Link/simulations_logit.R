require(dplyr)
require(mclust)
require(graphstats)
require(grdpg)


## Hyperparameter
addCovariates <- TRUE

K <- 2
d <- 1
p <- -1.5

cov <- 2


## Change n
n <- 2000


## Change q
q <- 1.5
latent <- cbind(p, q)


## Change size of each block
pi <- c(0.5, 0.5)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Change beta
beta <- 1    


## Change proportion of one binary covariate
pi_z <- c(0.5, 0.5)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[1], pi[2]*pi_z[2])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}


## Generate B and P matrix
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   


## Simulation
allresults <- data.frame()

imax <- 10
for (i in 1:imax) {
  
  tryCatch(
    {
      cat(i, '\n')
      
      A <- generateA(n, sigmoid(P), seed = i)
      
      result <- GRDPGwithCovariates(A, covariates, link = 'logit', dhat = 5, plot = FALSE)
      
      ARI <- adjustedRandIndex(blocks_cov, result$clusters_cov)
      model2 <- Mclust(rotate(result$Xhatprime[,1,drop=FALSE], latent, K), K, verbose = FALSE)
      
      phat <- model2$parameters$mean[1]
      qhat <- model2$parameters$mean[2]
      betahat <- result$betahat
      
      allresults <- rbind(allresults, data.frame(n, pi_z[1], pi_z[2], p, q, phat, qhat, beta, betahat, ARI))
    }, 
    error = function(e) { cat("\n", "ERROR: ", conditionMessage(e), "\n\n") }
  )
  
}




