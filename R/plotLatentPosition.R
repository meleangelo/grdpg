#' Plot Latent Position
#'
#' Plot latent position (with the effect of covariates).
#'
#' @import ggplot2
#'
#' @param Xhat Estimated latent position. Should be a `n` by `d` matrix where `n` is the number of nodes and `d` is the embeded dimension.
#' @param blocks Block label of each nodes. Should be a `n`-vector where `n` is the number of nodes. \code{NULL} by default.
#' @param withCovariates Whether with effect of covariates. `FALSE` suggests without effect of covariates.
#' @param ... Additional parameters based on the value of `withCovariates`
#' \itemize{
#' \item{With covariates}
#' \describe{
#' \item{`dhat`}{Embeded dimension}
#' \item{`covariates`}{Observed covariates}
#' }
#' \item{Without covariates}
#' \describe{
#' \item{`latent`}{True latent position}
#' \item{`K`}{Number of blocks}
#' \item{`d`}{Dimension of latent position}
#' }
#' }
#'
#' @return  A `ggplot2` object.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


plotLatentPosition <- function(Xhat, blocks = NULL, withCovariates, ...) {
  if (withCovariates) {
    dhat <- list(...)[[1]]
    covariates <- list(...)[[2]]
    pp2 <- vector('list', dhat-1)
    if (is.null(blocks)) {
      for (k in 1:ncol(covariates)) {
        for (l in 2:dhat) {
          dat <- data.frame(Xhat[,c(1,l)], Covariate = as.factor(covariates[,k]))
          pp2[[l-1]] <- ggplot(dat) + geom_point(aes(x=X1, y=X2, shape=Covariate), alpha = 0.5)
          pp2[[l-1]] <- pp2[[l-1]] + labs(title = 'Latent Position (with Covariates)', x = 'PC1', y = paste0('PC',l), shape = paste0('Covariate', k))
          print(pp2[[l-1]])
        }
      }
    } else {
      for (k in 1:ncol(covariates)) {
        for (l in 2:dhat) {
          dat <- data.frame(Xhat[,c(1,l)], Blocks = as.factor(blocks), Covariate = as.factor(covariates[,k]))
          pp2[[l-1]] <- ggplot(dat) + geom_point(aes(x=X1, y=X2, color=Blocks, shape=Covariate), alpha = 0.5)
          pp2[[l-1]] <- pp2[[l-1]] + labs(title = 'Latent Position (with Covariates)', x = 'PC1', y = paste0('PC',l), shape = paste0('Covariate', k))
          print(pp2[[l-1]])
        }
      }
    }
    multiplot(plotlist = pp2, cols = ceiling(length(pp2)/2))
    return(pp2)
  } else {
    latent <- list(...)[[1]]
    K <- list(...)[[2]]
    d <- list(...)[[3]]
    dat <- data.frame(rotate(Xhat, latent, K))
    if(d == 1) {
      names(dat) <- 'Xhatprime'
      if (is.null(blocks)) {
        pp4 <- ggplot(dat, aes(Xhatprime)) + geom_histogram(bins = 100)
        for (k in 1:length(latent)) {
          pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
        }
        pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
        print(pp4)
      } else {
        Blocks <- as.factor(blocks)
        pp4 <- ggplot(dat, aes(Xhatprime)) + geom_histogram(aes(fill=Blocks), bins = 100)
        for (k in 1:length(latent)) {
          pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
        }
        pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
        print(pp4)
      }
    } else {
      if (is.null(blocks)) {
        latent_vecs <- data.frame(t(latent))
        pp4 <- ggplot(dat) + geom_point(aes(X1, X2), size = 0.3, alpha = 0.5)
        pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
        pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
        print(pp4)
      } else {
        Blocks <- as.factor(blocks)
        latent_vecs <- data.frame(t(latent))
        pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
        pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
        pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
        print(pp4)
      }
    }
    return(pp4)
  }
}

