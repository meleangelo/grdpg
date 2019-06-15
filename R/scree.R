#' Scree Plot of Singular Values
#'
#' Scree plot of the singular values of a matrix.
#'
#' @param s Singular values of a matrix.
#' @param title Title of the plot. 'Screeplot' by default.
#' @param xlab Label of x-axis. 'Rank' by default.
#' @param ylab Label of y-axis. 'Singular Value' by default.
#'
#' @return A `ggplo2` object.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' X <- matrix(runif(100), nrow = 10)
#' dicomp <- svd(X)
#' scree(dicomp$d)
#'
#' @export


scree <- function(s, title = 'Screeplot', xlab = 'Rank', ylab = 'Singular Value') {
  dat <- data.frame(s)
  pp1 <- ggplot(dat, aes(x=1:nrow(dat), y=s)) + geom_line() + geom_point()
  pp1 <- pp1 + labs(title = title, x = xlab, y = ylab) + scale_x_continuous(breaks = 1:nrow(dat))
  print(pp1)
  return(pp1)
}

