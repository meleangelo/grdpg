#' Generate Probability Matrix
#'
#' Generate probability matrix from latent position (with the effect of covariates).
#'
#' @param latent True latent position. Suould be a `d` by `K` matrix where `d` is the dimension of latent position and `K` is the number of blocks.
#' @param d Dimension of latent position.
#' @param block_size Size of each block.
#' @param addCovariates Whether add the effect of covariates. `FALSE` suggests without the effect of covariates.
#' @param ... If `addCovariates` is `TRUE`, then pass two additional parameters `covariates` (observed covariates, `n` by `c` matrix or dataframe where `n` is the number of nodes and `c` is the number of covariates) and `beta` (effect of covariates).
#'
#' @return A `n` by `n` matrix where `n` is the number of nodes.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' ## Without covariates
#' # Hyperparameter
#' addCovariates <- FALSE
#' n <- 2000                                       # Number of nodes
#' K <- 2                                          # Number of blocks
#' d <- 2                                          # Dimension of latent position
#' latent <- cbind(c(0.63, -0.14), c(0.69, 0.13))  # Latent position
#'
#' # Size of each block is balanced
#' pi <- rep(1/K, K)
#' block_size <- round(pi * n)
#'
#' P <- generateP(latent, d, block_size, addCovariates)
#'
#'
#' ## With covariates
#' # Hyperparameter
#' addCovariates <- TRUE
#' n <- 2000                     # Number of nodes
#' K <- 2                        # Number of blocks
#' d <- 1                        # Dimension of latent position
#' latent <- cbind(0.1, 0.7)     # Latent position
#'
#' # Size of each block is balanced
#' pi <- rep(1/K, K)
#' block_size <- round(pi * n)
#'
#' # One binary covariate
#' beta <- 0.3           # If two covariates: c(0,1, 0.3), etc
#' cov <- 2              # If two binary: c(2, 2), etc
#' covariates <- matrix(rep(1:cov, each = 500, times = K))
#'
#' P <- generateP(latent, d, block_size, addCovariates, covariates, beta)
#'
#' @export


generateP <- function(latent, d, block_size, addCovariates, ...) {
  if (nrow(latent) != d) {
    stop("The number of rows in `latent` should equal to `d`.")
  }
  if (addCovariates && length(list(...)) != 2) {
    stop("There should be two more parameters (`covariates` and `beta`) if add covariates.")
  }
  X <- matrix(rep(latent[,1]), nrow = block_size[1], ncol = d, byrow = TRUE)
  for (k in 2:ncol(latent)) {
    X <- rbind(X, matrix(rep(latent[,k]), nrow = block_size[k], ncol = d, byrow = TRUE))
  }
  P <- X %*% t(X)
  if (addCovariates) {
    covariates <- list(...)[[1]]
    beta <- list(...)[[2]]
    for (k in 1:ncol(covariates)) {
      Beta <- matrix(0, nrow = nrow(P), ncol = ncol(P))
      for (i in 1:(nrow(Beta)-1)) {
        for (j in (i+1):ncol(Beta)) {
          Beta[i,j] <- ifelse(covariates[i,k]==covariates[j,k], 1, 0)
        }
      }
      Beta <- Beta + t(Beta)
      diag(Beta) <- 1
      Beta <- Beta * beta[k]
      P <- P + Beta
    }
  }
  return(P)
}
