#' Get the Block Label of Each Nodes
#'
#' Get the block label of each nodes according to the block assignment probability.
#'
#' @importFrom dplyr mutate
#'
#' @param block_assignment_probs Estimated block assignment probability. Should be an `n` by `k` matrix or dataframe where `n` is the number of nodes and `k` is the number of blocks.
#'
#' @return A vector with block labels of each nodes.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


getClusters <- function(block_assignment_probs) {
  block_assignment_probs <- data.frame(block_assignment_probs)
  block_assignment_probs <- mutate(block_assignment_probs, cluster = apply(block_assignment_probs, 1, which.max))
  assignments <- block_assignment_probs$cluster
  return(assignments)
}


