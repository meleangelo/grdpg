require(dplyr)
require(mclust)
require(graphstats)
require(grdpg)
require(Rmpi)


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
pi_1 <- 0.5
pi_2 <- 1 - pi_2
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Change beta
beta <- 1    


## Change proportion of one binary covariate
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


## Generate B and P matrix
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   


## Simulation
allresults <- data.frame()

imax <- 100
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
      
      allresults <- rbind(allresults, data.frame(n, pi_1, pi_2, pi_z_11, pi_z_21, p, q, phat, qhat, beta, betahat, ARI))
    }, 
    error = function(e) { cat("\n", "ERROR: ", conditionMessage(e), "\n\n") }
  )
  
}


## Using MPI
simulation_logit <- function(seed) {
  
  results <- data.frame()
  
  for (i in seed) {
    tryCatch(
      {
        A <- generateA(n, sigmoid(P), seed = i)
        
        result <- GRDPGwithCovariates(A, covariates, link = 'logit', dhat = 5, plot = FALSE)
        
        ARI <- adjustedRandIndex(blocks_cov, result$clusters_cov)
        model2 <- Mclust(rotate(result$Xhatprime[,1,drop=FALSE], latent, K), K, verbose = FALSE)
        
        phat <- model2$parameters$mean[1]
        qhat <- model2$parameters$mean[2]
        betahat <- result$betahat
        
        results <- rbind(results, data.frame(n, pi_1, pi_2, pi_z_11, pi_z_12, pi_z_21, pi_z_22, p, q, phat, qhat, beta, betahat, ARI))
      }, 
      error = function(e) { cat("\n", "ERROR: ", conditionMessage(e), "\n\n") }
    )
  }
  
  return(results)
  
}


mpi.spawn.Rslaves(nslaves = 3)

mpi.remote.exec(library(dplyr))
mpi.remote.exec(library(mclust))
mpi.remote.exec(library(grdpg))

mpi.bcast.cmd(n <- mpi.bcast.Robj())
mpi.bcast.Robj(n)
mpi.bcast.cmd(pi_1 <- mpi.bcast.Robj())
mpi.bcast.Robj(pi_1)
mpi.bcast.cmd(pi_2 <- mpi.bcast.Robj())
mpi.bcast.Robj(pi_2)
mpi.bcast.cmd(pi_z_11 <- mpi.bcast.Robj())
mpi.bcast.Robj(pi_z_11)
mpi.bcast.cmd(pi_z_12 <- mpi.bcast.Robj())
mpi.bcast.Robj(pi_z_12)
mpi.bcast.cmd(pi_z_21 <- mpi.bcast.Robj())
mpi.bcast.Robj(pi_z_21)
mpi.bcast.cmd(pi_z_22 <- mpi.bcast.Robj())
mpi.bcast.Robj(pi_z_22)
mpi.bcast.cmd(p <- mpi.bcast.Robj())
mpi.bcast.Robj(p)
mpi.bcast.cmd(q <- mpi.bcast.Robj())
mpi.bcast.Robj(q)
mpi.bcast.cmd(latent <- mpi.bcast.Robj())
mpi.bcast.Robj(latent)
mpi.bcast.cmd(K <- mpi.bcast.Robj())
mpi.bcast.Robj(K)
mpi.bcast.cmd(beta <- mpi.bcast.Robj())
mpi.bcast.Robj(beta)
mpi.bcast.cmd(blocks_cov <- mpi.bcast.Robj())
mpi.bcast.Robj(blocks_cov)
mpi.bcast.cmd(covariates <- mpi.bcast.Robj())
mpi.bcast.Robj(covariates)
mpi.bcast.cmd(P <- mpi.bcast.Robj())
mpi.bcast.Robj(P)

mpi.bcast.Robj2slave(simulation_logit)

seed_list <- list(10001:10010, 10011:10020, 10021:10030)
mpi.scatter.Robj2slave(seed_list)

temp <- mpi.remote.exec(simulation_logit(seed_list))

tempresult <- bind_rows(temp[[1]], temp[[2]], temp[[3]])
allresults <- bind_rows(allresults, tempresult)















