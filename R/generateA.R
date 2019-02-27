#' Generate Adjacency Matrix
#'
#' Generate adjacency matrix from probability matrix.
#'
#' @importFrom stats rbinom
#'
#' @param n Number of nodes.
#' @param P Probability matrix. Shoud be an `n` by `n` matrix.
#' @param directed Whether it is directed. `FALSE` by default (undirected).
#' @param type Type of network links. 'bernoulli' by default.
#' @param seed Random seed for reproductivity. 2019 by default.
#'
#' @return An `n` by `n` matrix where `n` is the number of nodes.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' n <- 4
#' P <- matrix(runif(n*n), nrow = n)
#' P[lower.tri(P)] <- t(P)[lower.tri(t(P))]
#' generateA(n, P)
#'
#' @export


generateA <- function(n, P, directed = FALSE, type = 'bernoulli', seed = 2019) {
  if (n <= 0 || !is.numeric(n)) {
    stop("`n` needs to be a positive number.")
  }
  if (!is.numeric(P)) {
    stop("`P` needs to be numeric.")
  }
  if (n != ncol(P) || n != nrow(P)) {
    stop("`n` should equal to the number of rows/columns in `P`.")
  }
  if (!(type %in% c('bernoulli', 'poisson'))) {
    print("Unrecognized `type`, would use 'bernoulli' by default.")
  }
  set.seed(seed)
  if (type == 'poisson') {
    A <- matrix(rpois(n*n,P), nrow = n)
  } else {
    A <- matrix(rbinom(n*n,1,P), nrow = n)
  }
  if (!directed) {
    A[lower.tri(A)] <- t(A)[lower.tri(t(A))]
  }
  return(A)
}

