#' Rotation of Estimated Latent Position
#'
#' Rotate estimated latent position to compare with the truth.
#'
#' @importFrom mclust Mclust
#'
#' @param Xhat Estimated latent position. Should be an `n` by `d` matrix where `n` is the number of nodes and `d` is the embeded dimension.
#' @param latent True latent position. Suould be a `d` by `K` matrix where `d` is the dimension of latent position and `K` is the number of blocks.
#' @param K Number of blocks.
#'
#' @return Rotated estimated latent position that could be compared with the truth.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


rotate <- function(Xhat, latent, K) {
  if (ncol(Xhat) != nrow(latent)) {
    stop("The number of columns in `Xhat` should equal to the number of rows in `latent`.")
  }
  if (ncol(latent) != K) {
    stop("The number of columns in `latent` should equal to K.")
  }
  model <- Mclust(Xhat, K, verbose = FALSE)
  means <- model$parameters$mean
  M <- svd(means %*% t(latent))
  R <- M$u %*% t(M$v)
  X_R <- Xhat %*% R
  return(X_R)
}

