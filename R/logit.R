#' Logit Function
#'
#' Logit function.
#'
#' \href{https://en.wikipedia.org/wiki/Logit}{Logit function}
#'
#' @param p Input of the logit function. Could be a number, vector or matrix.
#'
#' @return Logit of the input `p`. If `p` is a vector or matrix, it will return the element-wise logit of `p`.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' P <- matrix(runif(4), nrow = 2)
#' logit(P)
#'
#' @export


logit <- function(p) {
  if (any(p <= 0) | any(p >= 1)) {
    stop('All input need to lie between 0 and 1')
  }
  return(log(p/(1-p)))
}

