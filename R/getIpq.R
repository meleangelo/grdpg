#' Get `Ipq` Matrix for GRDPG
#'
#' Construct `Ipq` matrix for Generalized Random Dot Product Graph.
#'
#' @import RSpectra
#'
#' @param A A square matrix.
#' @param d Embeded dimension.
#' @param method Method used to compute eigenvalues. Coule be `RSpectra` (default) or `eigen`.
#'
#' @return A `d` by `d` diagonal matrix where `d` is the embeded dimension with 1 and -1 on the diagonal.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @references Rubin-Delanchy, P., Priebe, C. E., Tang, M., & Cape, J. (2017). A statistical interpretation of spectral embedding: the generalised random dot product graph. \emph{arXiv preprint \href{https://arxiv.org/abs/1709.05506}{arXiv:1709.05506}}.
#'
#' @seealso \code{\link{grdpg}}, \code{\link{eigs}}, \code{\link{eigen}}
#'
#' @export


getIpq <- function(A, d, method = 'RSpectra') {
  if (!(method %in% c('RSpectra', 'eigen'))) {
    print("Unrecognized `method`, would use 'RSpectra' by default.")
  }

  if (d == 1) {
    Ipq <- matrix(1)
  } else if (method == 'eigen') {
    temp <- eigen(A)
    s <- temp$values
    if (length(s) < d) {
      stop("The length of `s` should be greater than `d`.")
    }
    tempdat <- data.frame(raw = s) %>%
      mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
      arrange(desc(s))
    p <- sum(tempdat$raw[1:d]>0)
    if (p == d) {
      Ipq <- diag(rep(1,d))
    } else {
      Ipq <- diag(c(rep(1,p), rep(-1,d-p)))
    }
  } else {
    cols <- ncol(A)
    temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), d, 'LA')
    s1 <- temp1$values
    temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), d, 'SA')
    s2 <- temp2$values
    tempdat <- data.frame(raw = c(s1,s2)) %>%
      mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
      arrange(desc(s))
    p <- sum(tempdat$raw[1:d]>0)
    if (p == d) {
      Ipq <- diag(rep(1,d))
    } else {
      Ipq <- diag(c(rep(1,p), rep(-1,d-p)))
    }
  }

  return(Ipq)
}









