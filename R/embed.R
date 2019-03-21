#' Spectral Embedding
#'
#' Spectral decomposition of the adjacency/Laplacian matrices of graphs.
#'
#' @import graphstats
#' @importFrom igraph as_adj
#' @importFrom irlba irlba
#'
#' @param g An `\link[igraph]{igraph}` object or an `n` by `n` adjacency matrix with `n` vertices.
#' @param k Embeded dimension. \code{NULL} by default. `k` should be less than the number of vertices in `g`. If \code{k==NULL}, defaults to \code{floor(log2(nvertices))}.
#' @param edge.attr If `g` is a `\link[igraph]{igraph}`, the name of the attribute to use for weights. Defaults to \code{NULL}, which assumes the graph is binary.
#' \itemize{
#' \item{\code{is.null(edge.attr)} Assumes unweighted.}
#' \item{\code{is.character(edge.attr)} Assumes weighted with weights given by `edge.attr`.}
#' }
#' @param maxit Maximum number of iterations for `\link[irlba]{irlba}`.
#' @param work Working subspace dimension for `\link[irlba]{irlba}`.
#' @param tol Stopping tolerance for `\link[irlba]{irlba}`.
#' @param ... Other parameters for `\link[irlba]{irlba}`.
#'
#' @return A list containing the following:
#' \describe{
#' \item{`X`}{An `n` by `k` matrix indicating the estimated latent positions, where `n` is the number of vertices of `g`.}
#' \item{`Y`}{\code{NULL} if `g` is undirected. If `g` is directed, `Y` is a `n` by `k` matrix indicating the second half of the latent positions.}
#' \item{`D`}{The eigenvalues (for undirected graphs) or the singular values (for directed graphs) associated with the latent positions.}
#' \item{`iter`}{The number of Lanczos iterations carried out. See \link[irlba]{irlba}.}
#' \item{`mprod`}{The total number of matrix vector products carried out. See \link[irlba]{irlba}.}
#' }
#'
#' @references `\link[graphstats]{gs.embed}`, `\link[graphstats]{gs.embed.ase}`, `\link[graphstats]{gs.embed.lse}`, `\link[irlba]{irlba}`
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


embed <- function(g, k = NULL, edge.attr = NULL, maxit = 1000, work = 12, tol = 1e-05, ...) {
  if (class(g) == 'igraph') {
    A <- as.matrix(as_adj(g, type = 'both', attr = edge.attr))
  } else {
    tryCatch({
      A <- as.matrix(g)
    }, error=function(e) stop("You have passed neither an 'igraph' object nor an object that can be cast to a symmetric 'matrix'."))
    dimA <- dim(A)
    if (dimA[1] != dimA[2]) {
      stop("You have not passed a symmetric matrix.")
    }
  }

  if (is.null(k)) {
    k <- floor(log2(dim(A)[1]))
  }
  if (length(k) > 1) { stop("The number of embedding dimensions 'k' has length > 1.") }
  if (class(k) != "numeric" && !is.integer(k)) { stop("The number of embedding dimensions 'k' is not a number.") }
  if (k%%1 != 0) { stop("The number of embedding dimensions 'k' must be an integer.") }
  if (k < 1) { stop("The number of embedding dimensions 'k' < 1.") }
  if (k > dim(A)[1]) { stop("The number of embedding dimensions 'k' is greater than number of vertices.") }

  nu <- k
  if (all(t(A) == A)) {
    nv <- 0
  } else {
    nv <- k
  }
  result <- irlba(A, nu = nu, nv = nv, maxit = maxit, work = work, tol = tol)
  if (nv != 0) {
    result$Y <- result$v
  }
  result$X <- result$u
  result$D <- result$d
  result$u <- NULL
  result$v <- NULL
  result$d <- NULL
  return(result)
}


