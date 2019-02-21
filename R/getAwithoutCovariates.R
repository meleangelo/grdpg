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
  for (k in 1:ncol(covariates)) {
    Beta <- matrix(0, nrow = nrow(A), ncol = ncol(A))
    for (i in 1:nrow(Beta)) {
      for (j in 1:ncol(Beta)) {
        Beta[i,j] <- ifelse(covariates[i,k]==covariates[j,k], 1, 0)
      }
    }
    Beta <- Beta * betahat[k]
    Aprime <- Aprime - Beta
  }
  return(Aprime)
}


