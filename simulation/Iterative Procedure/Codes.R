require(dplyr)
require(graphstats)
require(mclust)
require(ggplot2)
require(grdpg)

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

get_beta <- function(Phat,covariates,Xhat){
  Z=covariates
  n=dim(A)[1]
  fz=0
  for (i in sequence(n)){
    for (j in sequence(n)){
      fz=fz+Z[i]*Z[j]*(Phat[i,j]-t(Xhat[i,])%*%Xhat[j,])
    }
  }
  fm=norm(Z%*%t(Z),"f")^2
  beta=fz/fm
  return(beta)
}



latent <- cbind(0.35, 0.65)
#latent <- cbind(0.2, 0.4, 0.7, 0.8)                # d = 1
#latent <- cbind(c(0.63, -0.14), c(0.69, 0.13))     # d = 2
beta <- 0.15
K <- 2
d <- 1
n <- 2000
pi <- rep(1/K, K)  # Balanced 
#pi <- c(0.3, 0.7)  # Unbalanced
#pi_cov <- rep(1/(K*2), K*2)                      # Balanced 
pi_cov <- c(0.5*0.2, 0.5*0.8, 0.5*0.4, 0.5*0.6)  # Unbalanced
block_size <- round(pi * n)
block_size_cov <- round(pi_cov * n)
cov <- 2     # Possible value of covariate, e.g. binary = 2
dmax <- 5
seed <- 2017

B=generateB(latent,K,d,addCovariates=TRUE,cov,beta)

covariates <- c()
ind <- 1
for (k in 1:length(block_size_cov)) {
  if(k %% cov == 0) {
    ind <- 2
  }
  covariates <- c(covariates, rep(ind, block_size_cov[k]))
  ind <- 1
}

covariates <- as.matrix(covariates)


P=generateP(latent, d, block_size, addCovariates=TRUE,covariates,beta)

A <- generateA(n, P, type = 'bernoulli', seed)
dhat=4
embed <- gs.embed.ase(A, dhat)
yhat <- embed$X %*% sqrt(diag(embed$D, nrow=dhat, ncol=dhat))

# temp <- eigen(A)
# Ipq <- getIpq(temp$values, dhat)
Phat=yhat %*% t(yhat)

set.seed(2019)
betahat=runif(1)
beta_c=c(betahat)
iter=25
for (k in 1:iter){
  embed <- gs.embed.ase(Phat-betahat*covariates%*%t(covariates), dhat)
  Xhat <- embed$X %*% sqrt(diag(embed$D, nrow=dhat, ncol=dhat))
  betahat <- get_beta(Phat,covariates,Xhat)
  betahat = as.numeric(betahat)
  beta_c = c(beta_c,betahat)
}
print(beta_c)

dat=data.frame(beta_c)
pp1 <- ggplot(dat, aes(x=0:iter, y=beta_c)) + geom_line() + geom_point()
pp1 <- pp1 + labs(title = 'True beta=0.15,pi_cov0=c(0.1,0.4,0.2,0.3)', x = 'iterations', y = 'estimated beta')
plot(pp1)
