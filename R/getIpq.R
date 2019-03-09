#' Get `Ipq` Matrix for GRDPG
#'
#' Construct `Ipq` matrix for Generalized Random Dot Product Graph.
#'
#' @param s Eigenvalues of a matrix.
#' @param d Embeded dimension.
#'
#' @return A `d` by `d` diagonal matrix where `d` is the embeded dimension with 1 and -1 on the diagonal.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @references Rubin-Delanchy, P., Priebe, C. E., Tang, M., & Cape, J. (2017). A statistical interpretation of spectral embedding: the generalised random dot product graph. \emph{arXiv preprint \href{https://arxiv.org/abs/1709.05506}{arXiv:1709.05506}}.
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


getIpq <- function(s, d) {
  if (length(s) < d) {
    stop("The length of `s` should be greater than `d`.")
  }
  p <- which(s<abs(s[length(s)]))[1]
  if (is.na(p)) {
    Ipq <- diag(rep(1,d))
  } else if (p > d) {
    Ipq <- diag(rep(1,d))
  } else {
    Ipq <- diag(c(rep(1,p-1), rep(-1, d-p+1)))
  }
  return(Ipq)
}


