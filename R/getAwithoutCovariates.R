#' Get the Adjacency Matrix without Covariates
#'
#' Get the adjacency matrix without the effect of covariates.
#'
#' @param A Adjacency matrix with the effect of covariates. Should be an `n` by `n` matrix where `n` is the number of nodes.
#' @param betahat Estimated beta. Should be a `c`-vector where `c` is the number of covariates.
#' @param covariates Observed covariates. Should be a `n` by `c` matrix or dataframe where `n` is the number of nodes and `c` is the number of covariates.
#'
#' @return An adjacency matrix removing the effect of covariates.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


getAwithoutCovariates <- function(A, betahat, covariates) {
  Aprime <- A
  covariates <- data.frame(covariates)
  if (nrow(A) != nrow(covariates) || ncol(A) != nrow(covariates)) {
    stop("The number of rows/columns in `A` should equal to the number of rows in `covariates`.")
  }
  if (length(betahat) != ncol(covariates)) {
    stop("The length of `betahat` should equal to the number of columns in `covariates`.")
  }
  for (k in 1:ncol(covariates)) {
    Beta <- matrix(0, nrow = nrow(A), ncol = ncol(A))
    for (i in 1:(nrow(Beta)-1)) {
      for (j in (i+1):ncol(Beta)) {
        Beta[i,j] <- ifelse(covariates[i,k]==covariates[j,k], 1, 0)
      }
    }
    Beta <- Beta + t(Beta)
    diag(Beta) <- 1
    Beta <- Beta * betahat[k]
    Aprime <- Aprime - Beta
  }
  return(Aprime)
}


