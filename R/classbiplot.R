#' Grouped biplot for PLS scores and loadings
#'
#' This helper plots the individuals and predictors from a fitted
#' \code{\link{plsR}} or \code{\link{plsRglm}} model while coloring the
#' individuals according to a grouping variable.
#'
#' @param object an object containing score and loading matrices in
#' \code{object$tt} and \code{object$pp}, typically returned by
#' \code{\link{plsR}} or \code{\link{plsRglm}}.
#' @param group optional grouping vector for the individuals. When supplied,
#' observations are colored according to the levels of \code{group}.
#' @param comps integer vector of length 2 giving the components to display.
#' @param col colors for the individuals. If \code{group} is provided,
#' \code{col} can have length 1, the number of groups, or the number of
#' individuals.
#' @param colvar color used for variable labels, arrows and axes.
#' @param pch plotting character for the individuals.
#' @param cex character expansion. Length 1 is recycled to length 2: the first
#' value is used for the individuals and the second one for the variables.
#' @param xlabs optional labels for the individuals.
#' @param ylabs optional labels for the variables.
#' @param point.labels shall the individuals be displayed using text labels
#' instead of points? Defaults to \code{FALSE}.
#' @param show.legend shall a legend be added when \code{group} is provided?
#' Defaults to \code{TRUE}.
#' @param legendpos position of the legend as in
#' \code{\link[graphics:legend]{legend}}, defaults to \code{"topright"}.
#' @param var.axes shall arrows be drawn for the variables? Defaults to
#' \code{TRUE}.
#' @param expand expansion factor for the variables layer, as in
#' \code{\link[stats:biplot]{biplot}}. Defaults to \code{1}.
#' @param xlim,ylim limits for the scores panel. When both are missing they are
#' chosen symmetrically as in \code{biplot}.
#' @param arrow.len length of the arrows for the variables.
#' @param main,sub,xlab,ylab usual graphical parameters passed to
#' \code{\link[graphics:plot]{plot}}.
#' @param \dots further graphical parameters passed to
#' \code{\link[graphics:plot]{plot}} for the scores panel.
#'
#' @return Invisibly returns a list with the scores, loadings, colors, grouping
#' factor and scaling ratio used in the plot.
#'
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#'
#' @seealso \code{\link{plsR}}, \code{\link{plsRglm}},
#' \code{\link[stats:biplot]{biplot}}
#'
#' @keywords regression models
#'
#' @examples
#' data(Cornell)
#' modpls <- plsR(Y ~ ., data = Cornell, nt = 2)
#' grp <- factor(Cornell$Y > median(Cornell$Y), labels = c("Low", "High"))
#' classbiplot(modpls, group = grp, col = c("firebrick3", "steelblue3"))
#'
#' @export classbiplot
classbiplot <- function(object, group = NULL, comps = 1:2, col,
                        colvar = "gray30", pch = 19,
                        cex = rep(par("cex"), 2), xlabs = NULL, ylabs = NULL,
                        point.labels = FALSE, show.legend = TRUE,
                        legendpos = "topright", var.axes = TRUE, expand = 1,
                        xlim = NULL, ylim = NULL, arrow.len = 0.1,
                        main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
                        ...) {
  if (is.null(object$tt) || is.null(object$pp)) {
    stop("'object' must contain 'tt' and 'pp' components.")
  }

  if (length(comps) != 2L) {
    stop("'comps' should have length 2.")
  }

  comps <- as.integer(comps)
  if (any(is.na(comps)) || any(comps < 1L)) {
    stop("'comps' should contain two positive integers.")
  }

  tt <- as.matrix(object$tt)
  pp <- as.matrix(object$pp)

  if (ncol(tt) < max(comps) || ncol(pp) < max(comps)) {
    stop("Requested components are not available in 'object'. Fit at least max(comps) components first.")
  }

  tt <- tt[, comps, drop = FALSE]
  pp <- pp[, comps, drop = FALSE]

  n <- nrow(tt)
  p <- nrow(pp)

  if (is.null(group)) {
    group_factor <- NULL
  } else {
    if (length(group) != n) {
      stop("'group' should have the same length as the number of individuals.")
    }
    if (anyNA(group)) {
      stop("'group' should not contain missing values.")
    }
    group_factor <- factor(group)
  }

  if (length(cex) == 1L) {
    cex <- c(cex, cex)
  }

  if (is.null(xlabs)) {
    xlabs <- rownames(tt)
    if (is.null(xlabs)) {
      xlabs <- seq_len(n)
    }
  }
  xlabs <- as.character(xlabs)

  if (is.null(ylabs)) {
    ylabs <- rownames(pp)
    if (is.null(ylabs)) {
      ylabs <- paste("Var", seq_len(p))
    }
  }
  ylabs <- as.character(ylabs)

  if (missing(col)) {
    if (is.null(group_factor)) {
      point_col <- rep(par("col"), n)
      legend_col <- NULL
    } else {
      legend_col <- grDevices::palette()
      if (nlevels(group_factor) > length(legend_col)) {
        legend_col <- grDevices::rainbow(nlevels(group_factor))
      } else {
        legend_col <- legend_col[seq_len(nlevels(group_factor))]
      }
      point_col <- legend_col[as.integer(group_factor)]
    }
  } else if (length(col) == 1L) {
    point_col <- rep(col, n)
    legend_col <- if (is.null(group_factor)) NULL else rep(col, nlevels(group_factor))
  } else if (!is.null(group_factor) && length(col) == nlevels(group_factor)) {
    legend_col <- col
    point_col <- legend_col[as.integer(group_factor)]
  } else if (length(col) == n) {
    point_col <- col
    legend_col <- if (is.null(group_factor)) NULL else vapply(split(point_col, group_factor), `[`, character(1), 1L)
  } else {
    stop("'col' should have length 1, the number of groups, or the number of individuals.")
  }

  if (is.null(xlab)) {
    xlab <- colnames(tt)[1L]
    if (is.null(xlab) || is.na(xlab) || !nzchar(xlab)) {
      xlab <- paste("Comp", comps[1L])
    }
  }
  if (is.null(ylab)) {
    ylab <- colnames(tt)[2L]
    if (is.null(ylab) || is.na(ylab) || !nzchar(ylab)) {
      ylab <- paste("Comp", comps[2L])
    }
  }

  unsigned.range <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) {
      return(c(-1, 1))
    }
    c(-abs(min(x)), abs(max(x)))
  }

  rangx1 <- unsigned.range(tt[, 1L])
  rangx2 <- unsigned.range(tt[, 2L])
  rangy1 <- unsigned.range(pp[, 1L])
  rangy2 <- unsigned.range(pp[, 2L])

  if (is.null(xlim) && is.null(ylim)) {
    xlim <- ylim <- range(rangx1, rangx2)
  } else if (is.null(xlim)) {
    xlim <- rangx1
  } else if (is.null(ylim)) {
    ylim <- rangx2
  }

  ratio_candidates <- c(rangy1 / rangx1, rangy2 / rangx2)
  ratio_candidates <- ratio_candidates[is.finite(ratio_candidates) & ratio_candidates > 0]
  ratio <- if (length(ratio_candidates)) max(ratio_candidates) / expand else 1
  if (!is.finite(ratio) || ratio <= 0) {
    ratio <- 1
  }

  op <- par(pty = "s")
  on.exit(par(op))

  plot(tt, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       sub = sub, main = main, ...)
  if (point.labels) {
    text(tt, labels = xlabs, cex = cex[1L], col = point_col)
  } else {
    points(tt[, 1L], tt[, 2L], pch = pch, cex = cex[1L], col = point_col)
  }

  par(new = TRUE)
  plot(pp, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim * ratio,
       xlab = "", ylab = "", col = colvar)
  axis(3, col = colvar)
  axis(4, col = colvar)
  graphics::box(col = par("fg"))
  text(pp, labels = ylabs, cex = cex[2L], col = colvar)
  if (var.axes) {
    arrows(0, 0, pp[, 1L] * 0.8, pp[, 2L] * 0.8, col = colvar,
           length = arrow.len)
  }

  if (show.legend && !is.null(group_factor) && nlevels(group_factor) > 1L) {
    legend(legendpos, legend = levels(group_factor), col = legend_col,
           pch = pch, bty = "n")
  }

  invisible(list(scores = tt, loadings = pp, colours = point_col,
                 group = group_factor, ratio = ratio))
}
