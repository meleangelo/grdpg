#' Remove Check of Probability Matrix
#'
#' Check whether all entries of the estimated probabiltiy matrix are between 0 and 1.
#' Remove all rows and columns that contain values that are not greater than 0 and not less than 1.
#'
#' @param P Estimated probability matrix.
#'
#' @return A matrix of which all entries are between 0 and 1.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' P <- matrix(runif(100,-1,2), nrow = 10)
#' Removecheck(P)
#'
#' @export


Removecheck <- function(P) {
  if (!is.numeric(P)) {
    stop("Input need to be numeric.")
  }
  ind <- which(P<=0|P>=1, arr.ind = T)
  removed <- unique(c(ind[,1],ind[,2]))
  return(P[-removed,-removed])
}

