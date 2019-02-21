#' Get the Covariates for Each Block
#'
#' Get the covariates for each block based on the frequency of each covariate within each block.
#'
#' @param covariates Observed covariates. Should be an `n` by `c` matrix or dataframe where `n` is the number of nodes and `c` is the number of covariates.
#' @param clusters_cov Block label for each nodes with the effect of covariates. Should be a `n`-vector with `k` unique values where `n` is the number of nodes and `k` is the number of blocks with the effect of covariates.
#'
#' @return A `k` by `c` matrix where `k` is the number of blocks with the effect of covariates, `c` is the number of covariates and entry `[i,j]` is the covariate `j` for block `i`.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


getBlockCovariates <- function(covariates, clusters_cov) {
  covariates_block <- matrix(0, nrow = length(unique(clusters_cov)), ncol = ncol(covariates))
  for (k in 1:ncol(covariates_block)) {
    for (l in 1:nrow(covariates_block)) {
      temp <- which(clusters_cov==l)
      covariates_block[l,k] <- as.integer(names(which.max(table(covariates[temp,k]))))
    }
  }
  return(covariates_block)
}


