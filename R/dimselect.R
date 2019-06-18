#' Dimension Selection
#'
#' Select embeded dimension of latent position by finding the \sQuote{elbow} of the scree plot using \href{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf}{profile likelihood}.
#'
#' @param s A vector of ordered singular values.
#' @param n Number of elbows to be returned. 3 by default.
#'
#' @return A list containing the following:
#' \describe{
#' \item{`value`}{The singular value associated with the elbow.}
#' \item{`elbow`}{The embeded dimension (index of the elbow).}
#' }
#'
#' @references `\link[graphstats]{gs.dim.select}`
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


dimselect <- function(s, n = 3) {
  if (!is.numeric(s)) {
    stop("Input need to be numeric (a vector of ordered singular values).")
  }
  d <- s
  p <- length(s)
  lq <- rep(0.0, p)
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) / (p - 1 - (q < p))
    lq[q] <- sum(dnorm(d[1:q ], mu1, sqrt(sigma2), log=TRUE)) + sum(dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE))
  }
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + dimselect(d[(q+1):p], n=n-1)$elbow)
  }
  out <- list(value = d[q], elbow = q)
  return(out)
}

