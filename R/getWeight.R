#' Get the Weight of Each Group
#'
#' Use frequency as the weight of each group.
#'
#' @import dplyr
#'
#' @param clusters An `n`-vector indicating the label of group assignments.
#'
#' @return A data frame summarizing the frequecny of each group.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @export


getWeight <- function(clusters) {
  n <- length(clusters)
  weight <- data.frame(clusters) %>%
    group_by(clusters) %>%
    summarise(count = n()) %>%
    mutate(freq = count/n)
  return(weight)
}

