#' Scree Plot of Eigenvalues (Singular Values)
#'
#' Scree plot of the eigenvalues (singular values) of a matrix.
#'
#' @import RSpectra
#'
#' @param s Eigenvalues (singular values) of a matrix.
#' @param title Title of the plot. 'Screeplot' by default.
#' @param xlab Label of x-axis. 'Rank in Magnitude' by default.
#' @param ylab Label of y-axis. 'Eigenvalue in Magnitude' by default.
#'
#' @return A `ggplo2` object.
#'
#' @author Cong Mu \email{placid8889@gmail.com}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' X <- matrix(runif(100), nrow = 10)
#' X <- X + t(X)
#' dicomp <- svd(X)
#' dicomp2 <- eigen(X)
#' scree(dicomp$d, xlab = 'Rank', ylab = 'Singular Value')
#' scree(dicomp2$values)
#'
#' @export


scree <- function(s, title = 'Screeplot', xlab = 'Rank in Magnitude', ylab = 'Eigenvalue in Magnitude') {
  # dat <- data.frame(s)
  # pp1 <- ggplot(dat, aes(x=1:nrow(dat), y=s)) + geom_line() + geom_point()
  # pp1 <- pp1 + labs(title = title, x = xlab, y = ylab) + scale_x_continuous(breaks = 1:nrow(dat))
  # print(pp1)

  dat <- data.frame(raw = s) %>%
    mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
    arrange(desc(s))
  pp1 <- ggplot(dat) + geom_line(aes(x=1:nrow(dat), y=s), linetype='dotted')
  pp1 <- pp1 + geom_point(aes(x=1:nrow(dat), y=s, color=sign, shape=sign), size=2)
  pp1 <- pp1 + labs(title = title, x = xlab, y = ylab) + scale_x_continuous(breaks = 1:nrow(dat))
  print(pp1)

  return(pp1)
}

