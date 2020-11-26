#' Intermediate Function
#'
#' Intermediate function.
#'
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#'
#' @export

compute_cov <- function(i, l1, l2, K, meanE, muhats, deltainv, eta, theta, zeta, Ipq) {
  .Call("_grdpg_compute_cov", PACKAGE = "grdpg", i, l1, l2, K, meanE, muhats, deltainv, eta, theta, zeta, Ipq)
}
