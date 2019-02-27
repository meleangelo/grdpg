#' Estimate beta
#'
#' Estimate beta (the effect of covariates).
#'
#' @import mclust
#'
#' @param BXhat Estimated block probability matrix. Should be a `k` by `k` matrix where `k` is the number of blocks.
#' @param cov A vector specifying possible value that each covariate could take. For example, if two binary covariates, then \code{cov <- c(2, 2)}.
#' @param covariates_block Estimated covariates for each block. Should be a `k` by `c` matrix or dataframe where `k` is the number of blocks and `c` is the number of covariates.
#'
#' @return A list containing all estimated beta. Each element of the list, i.e. the length of the list equals to the number of covariates.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


estimatebeta <- function(BXhat, cov, covariates_block) {
  covariates_block <- data.frame(covariates_block)
  if (length(cov) != ncol(covariates_block)) {
    stop("The length of `cov` should equal to the number of columns in `covariates_block` (both equal to the number of covariates).")
  }
  betahats <- vector('list', ncol(covariates_block))
  model2 <- Mclust(diag(BXhat), ncol(BXhat)/prod(cov), verbose = FALSE)
  c <- getClusters(data.frame(model2$z))
  for (i in 1:(nrow(BXhat)-1)) {
    for (j in (i+1):ncol(BXhat)) {
      for (k in 1:ncol(covariates_block)) {
        if (c[i] == c[j] & covariates_block[i,k] != covariates_block[j,k]) {
          temp <- setdiff(1:ncol(covariates_block),k)
          ind <- c()
          if (length(temp) > 0) {
            for (l in temp) {
              ind <- c(ind, covariates_block[i,l] == covariates_block[j,l])
            }
          } else {
            ind <- TRUE
          }
          if (all(ind)) {
            betahats[[k]] <- c(betahats[[k]], BXhat[i,i] - BXhat[i,j])
          }
        }
      }
    }
  }
  return(betahats)
}


