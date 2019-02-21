#' Sigmoid Function
#'
#' Sigmoid function.
#'
#' \href{https://en.wikipedia.org/wiki/Sigmoid_function}{Sigmoid function}
#'
#' @param x Input of the sigmoid function. Could be a number, vector or matrix.
#'
#' @return Sigmoid of the input `x`. If `x` is a vector or matrix, it will return the element-wise sigmoid of `x`.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' X <- matrix(runif(4), nrow = 2)
#' sigmoid(X)
#'
#' @export


sigmoid <- function(x) {
  return(1/(1+exp(-x)))
}

