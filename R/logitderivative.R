#' Derivative of logit Function
#'
#' Derivative of logit Function.
#'
#' @param x Input of the derivative of logit function. Could be a number, vector or matrix.
#'
#' @return Derivative of logit of the input `x`. If `x` is a vector or matrix, it will return the element-wise derivative of logit of `x`.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' X <- matrix(runif(4), nrow = 2)
#' logitderivative(X)
#'
#' @export


logitderivative <- function(x) {
  return(1/x + 1/(1-x))
}

