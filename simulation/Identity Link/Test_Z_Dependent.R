require(dplyr)
require(graphstats)
require(mclust)
require(ggplot2)




## B matrix
generateBwithCovariate <- function(latent, beta, K, d) {
  X <- matrix(rep(latent[,1]), nrow = cov, ncol = d, byrow = TRUE)
  for (k in 2:ncol(latent)) {
    X <- rbind(X, matrix(rep(latent[,k]), nrow = cov, ncol = d, byrow = TRUE))
  }
  
  Z <- rep(c(1,-1), times = K)
  ones=rep(1,times=K*2)
  Beta <- beta*(Z %*% t(Z) + ones %*% t(ones))/2
  B <- X %*% t(X) + Beta
  return(B)
}


## P matrix
generatePwithCovariate <- function(latent, beta, K, d, n, block_size, cov, covariates) {
  X <- matrix(rep(latent[,1]), nrow = block_size[1], ncol = d, byrow = TRUE)
  for (k in 2:ncol(latent)) {
    X <- rbind(X, matrix(rep(latent[,k]), nrow = block_size[k], ncol = d, byrow = TRUE))
  }
  # Z <- rep(1:cov, each = n/(K*cov), times = K)
  # Beta <- Z %*% t(Z)
  Beta <- covariates %*% t(covariates)
  Beta[Beta %in% (1:cov)^2] <- 1
  Beta[!(Beta %in% (1:cov)^2)] <- 0
  Beta <- Beta * beta
  P <- X %*% t(X) + Beta
  return(P)
}


## Generate adjacency matrix
generateA <- function(n, P, seed) {
  set.seed(seed)
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      A[i,j] <- rbinom(1, 1, P[i,j])
    }
  }
  A <- A + t(A)
  return(A)
}


## Get Ipq for GRDPG
getIpq <- function(s, d) {
  p <- which(s<abs(s[length(s)]))[1]
  if (p > d) {
    Ipq <- diag(rep(1,d))
  } else {
    Ipq <- diag(c(rep(1,p-1), rep(-1, d-p+1)))
  }
  return(Ipq)
}


## Estimate beta
estimatebeta <- function(BXhat, cov) {
  betahats <- c()
  model2 <- Mclust(diag(BXhat), ncol(BXhat)/cov, verbose = FALSE)
  c <- getClusters(model2)
  for (i in 1:(ncol(BXhat)/cov)) {
    ind <- which(c==i)
    for (j in 1:(length(ind)-1)) {
      betahats <- c(betahats, abs(BXhat[ind[j],]-BXhat[ind[j+1],]))
    }
  }
  symbol <- ifelse(BXhat[1,1]-BXhat[1,which(c==c[1])[2]]>0, 1, -1)
  betahat <- mean(betahats) * symbol
  return(betahat)
}





## Get the label of clusters
getClusters <- function(model) {
  block_assignment_probs <- data.frame(model$z)
  block_assignment_probs <- mutate(block_assignment_probs, cluster = apply(block_assignment_probs, 1, which.max))
  assignments <- block_assignment_probs$cluster
  return(assignments)
}




## Simulation
simulation_GRDPGwithCovariates <- function(latent, beta, K, d, n, block_size, block_size_cov, cov, dmax, sd=FALSE, seed=2018) {
  
  
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
  covariates <- c()
  ind <- 1
  for (k in 1:length(block_size_cov)) {
    if(k %% cov == 0) {
      ind <- 2
    }
    covariates <- c(covariates, rep(ind, block_size_cov[k]))
    ind <- 1
  }
  B <- generateBwithCovariate(latent, beta, K, d)
  P <- generatePwithCovariate(latent, beta, K, d, n, block_size, cov, covariates)
  
  cat('\n\n', 'Sampling...')
  
  
  ## Sample Network
  # g <- igraph::sample_sbm(n, B, block_size_cov)
  # A <- as_adj(g)
  
  A <- generateA(n, P, seed)
  
  cat('\n\n', 'Embedding...')
  
  ptm <- proc.time()
  
  ## ASE
  s <- gs.embed.ase(A, dmax)$D
  
  
  ## Embed
  dhat=4
  embed <- gs.embed.ase(A, dhat)
  Xhat <- embed$X %*% sqrt(diag(embed$D, nrow=dhat, ncol=dhat))
  
  
  
 
  
  ## Estiamte beta
  model <- Mclust(Xhat, verbose = FALSE) 
  muhats <- model$parameters$mean
  
  temp <- eigen(A)
  Ipq <- getIpq(temp$values, dhat)
  BXhat <- t(muhats) %*% Ipq %*% muhats
  
  betahat <- estimatebeta(BXhat, cov) 
  
  
  

  
  
  
  
  
  
  return(betahat)
}



## Example
seed <- 2018
latent <- cbind(0.425, 0.525)
#latent <- cbind(0.2, 0.4, 0.7, 0.8)                # d = 1
#latent <- cbind(c(0.63, -0.14), c(0.69, 0.13))     # d = 2
beta <- 0.15
K <- 2
d <- 1
n <- 2000
pi <- rep(1/K, K)  # Balanced 
#pi <- c(0.3, 0.7)  # Unbalanced
#pi_cov <- rep(1/(K*2), K*2)                      # Balanced 
#pi_cov <- c(0.3*0.4, 0.3*0.6, 0.7*0.4, 0.7*0.6)  # Unbalanced
block_size <- round(pi * n)
#block_size_cov <- round(pi_cov * n)
cov <- 2     # Possible value of covariate, e.g. binary = 2
dmax <- 5


## Simulation
pi_cov0=c(0.25,0.25,0.05,0.5-0.05)
bandwidth=0.025
block_size_cov <- round(pi_cov0 * n)
betahat=simulation_GRDPGwithCovariates(latent, beta, K, d, n, block_size, block_size_cov, cov, dmax)
error=abs(betahat-beta)
for (pi_cov_ in seq(0.05+bandwidth,0.45,bandwidth)){
  pi_cov=c(0.25,0.25,pi_cov_,0.5-pi_cov_)
  block_size_cov <- round(pi_cov * n)
  betahat=simulation_GRDPGwithCovariates(latent, beta, K, d, n, block_size, block_size_cov, cov, dmax)
  error=c(error,abs(betahat-beta))
}
plot(seq(0.05,0.45,0.025)*2,error,xlab='Proportion of first covariate in second block',ylab = 'beta error',main = 'n=2000, p=0.425, beta=0.15')
