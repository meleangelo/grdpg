#' Multiple Plots on One Page for `ggplot2`
#'
#' Show multiple plots on one page for `ggplot2`.
#'
#' @import grid
#' @import ggplot2
#'
#' @param ... A list of `ggplot2` objects.
#' @param plotlist A list of `ggplot2` objects. \code{NULL} by default.
#' @param cols Number of columns in layout. 1 by default.
#' @param layout A matrix specifying the layout. \code{NULL} by default. If present, `cols` is ignored. For example, if the layout is \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)}, then plot 1 will go in the upper left, 2 will go in the upper right, and 3 will go all the way across the bottom.
#'
#' @references \href{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}{Cookbook for R}
#'
#' @seealso \code{\link{grdpg}}
#'
#' @examples
#' require(ggplot2)
#' p1 <- ggplot(economics, aes(date, unemploy)) + geom_line()
#' p2 <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' multiplot(p1, p2)
#'
#' @export


multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {

  # Collect all plots
  plots <- c(list(...), plotlist)

  # Number of the plots
  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns in layout
    # nrow: Number of rows needed, calculated from number of columns
    layout <- matrix(seq(1, cols*ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the corresponding location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}
