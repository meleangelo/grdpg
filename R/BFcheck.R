#' Brute-Force Check of Probability Matrix
#'
#' Check whether all entries of the estimated probabiltiy matrix are between 0 and 1.
#' Set the ones that not greater than 0 to be a lower bound and the ones that not less than 1 to be a upper bound.
#'
#' @param P Estimated probability matrix.
#' @param l Lower bound of the probability matrix, i.e. set all entries that not greater than 0 to `l`. 0.0001 by default.
#' @param u Upper bound of the probability matrix, i.e. set all entries that not less than 1 to be `u`. 0.9999 by default.
#'
#' @return A matrix of which all entries are between 0 and 1.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' P <- matrix(runif(25,-1,2), nrow = 5)
#' BFcheck(P)
#'
#' @export


BFcheck <- function(P, l = 0.0001, u = 0.9999) {
  if (!is.numeric(P)) {
    stop("Input need to be numeric.")
  }
  P[P<=0] <- l
  P[P>=1] <- u
  return(P)
}

