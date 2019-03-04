#' Generate Block Probability Matrix
#'
#' Generate block probability matrix from latent position (with the effect of covariates).
#'
#' @param latent True latent position. Suould be a `d` by `K` matrix where `d` is the dimension of latent position and `K` is the number of blocks.
#' @param K Number of blocks.
#' @param d Dimension of latent position.
#' @param addCovariates Whether add the effect of covariates. \code{FALSE} suggests without the effect of covariates.
#' @param ... If `addCovariates` is \code{TRUE}, then pass two additional parameters `cov` (a vector specifying possible value that each covariate could take) and `beta` (effect of covariates).
#'
#' @return A \code{K*prod(cov)} by \code{K*prod(cov)} matrix (assuming discrete covariates) where `K` is the number of blocks and `cov` is a vector specifying possible value that each covariate could take (1 if without effect of covariates).
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' ## Without covariates
#' addCovariates <- FALSE
#' K <- 2                                           # Number of blocks
#' d <- 2                                           # Dimension of latent position
#' latent <- cbind(c(0.63, -0.14), c(0.69, 0.13))   # Latent position
#'
#' B <- generateB(latent, K, d, addCovariates)
#'
#'
#' ## With covariates
#' # Hyperparameter
#' addCovariates <- TRUE
#' K <- 2                        # Number of blocks
#' d <- 1                        # Dimension of latent position
#' latent <- cbind(0.1, 0.7)     # Latent position
#'
#' # One binary covariate
#' beta <- 0.3           # If two covariates: c(0,1, 0.3), etc
#' cov <- 2              # If two binary: c(2, 2), etc
#'
#' B <- generateB(latent, K, d, addCovariates, cov, beta)
#'
#' @export


generateB <- function(latent, K, d, addCovariates, ...) {
  if (nrow(latent) != d || ncol(latent) != K) {
    stop("`latent` should be a d by K matrix.")
  }
  if (addCovariates && length(list(...)) != 2) {
    stop("There should be two more parameters (`cov` and `beta`) if add covariates.")
  }
  if (addCovariates) {
    cov <- list(...)[[1]]
    beta <- list(...)[[2]]
    if (length(cov) != length(beta)) {
      stop("The length of `cov` should equal to the length of `beta`.")
    }
  } else {
    cov <- 1
  }
  X <- matrix(rep(latent[,1]), nrow = prod(cov), ncol = d, byrow = TRUE)
  for (k in 2:ncol(latent)) {
    X <- rbind(X, matrix(rep(latent[,k]), nrow = prod(cov), ncol = d, byrow = TRUE))
  }
  B <- X %*% t(X)
  if (addCovariates) {
    for (k in 1:length(beta)) {
      Beta <- matrix(0, nrow = nrow(B), ncol = ncol(B))
      Z <- rep(1:cov[k], each = prod(cov[k:length(cov)])/cov[k], times = K*prod(cov[1:k])/cov[k])
      for (i in 1:nrow(Beta)) {
        for (j in 1:ncol(Beta)) {
          Beta[i,j] <- ifelse(Z[i]==Z[j], 1, 0)
        }
      }
      Beta <- Beta * beta[k]
      B <- B + Beta
    }
  }
  return(B)
}
