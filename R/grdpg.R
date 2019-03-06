#' grdpg Package
#'
#' Generate and estimate (Generalized) Random Dot Product Graph with observed covariates to model the networks.
#'
#' Functions in this package are summarised as follows.
#' \itemize{
#' \item{Utility functions}
#' \describe{
#' \item{\code{\link{logit}}}{Logit function.}
#' \item{\code{\link{sigmoid}}}{Sigmoid function.}
#' \item{\code{\link{BFcheck}}}{Check whether all entries of the estimated probabiltiy matrix are between 0 and 1. Set the ones that not greater than 0 to be a lower bound and the ones that not less than 1 to be a upper bound.}
#' \item{\code{\link{Removecheck}}}{Check whether all entries of the estimated probabiltiy matrix are between 0 and 1. Remove all rows and columns that contain values that are not greater than 0 and not less than 1.}
#' }
#' \item{Simulation functions}
#' \describe{
#' \item{\code{\link{generateB}}}{Generate block probability matrix from latent position (with the effect of covariates).}
#' \item{\code{\link{generateP}}}{Generate probability matrix from latent position (with the effect of covariates).}
#' \item{\code{\link{generateA}}}{Generate adjacency matrix from probability matrix.}
#' \item{\code{\link{rotate}}}{Rotate estimated latent position to compare with the truth.}
#' \item{\code{\link{simulation}}}{Simulate and estimate (Generalized) Random Dot Product Graph (GRDPG) with effect of covariates.}
#' }
#' \item{Estimation functions}
#' \describe{
#' \item{\code{\link{dimselect}}}{Select embeded dimension of latent position by finding the \sQuote{elbow} of the scree plot using \href{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.90.3768&rep=rep1&type=pdf}{profile likelihood}.}
#' \item{\code{\link{embed}}}{Spectral decomposition of the adjacency/Laplacian matrices of graphs.}
#' \item{\code{\link{getIpq}}}{Construct `Ipq` matrix for Generalized Random Dot Product Graph.}
#' \item{\code{\link{getBlockCovariates}}}{Get the covariates for each block based on the frequency of each covariate within each block.}
#' \item{\code{\link{estimatebeta}}}{Estimate beta (the effect of covariates).}
#' \item{\code{\link{getAwithoutCovariates}}}{Get the adjacency matrix without the effect of covariates.}
#' \item{\code{\link{getClusters}}}{Get the block label of each nodes according to the block assignment probability.}
#' \item{\code{\link{GRDPGwithoutCovariates}}}{Estimate (Generalized) Random Dot Product Graph (GRDPG) without effect of covariates.}
#' \item{\code{\link{GRDPGwithCovariates}}}{Estimate (Generalized) Random Dot Product Graph (GRDPG) with effect of covariates.}
#' }
#' \item{Visualization functions}
#' \describe{
#' \item{\code{\link{multiplot}}}{Show multiple plots on one page for `ggplot2`.}
#' \item{\code{\link{screeplot}}}{Scree plot of the singular values of a matrix.}
#' \item{\code{\link{plotLatentPosition}}}{Plot latent position (with the effect of covariates).}
#' }
#' }
#'
#'
#' @docType package
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @name grdpg
NULL
