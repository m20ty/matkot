#' Head/tail left/right
#'
#' Functions to view the top/bottom left/right corners of large matrices or data frames.
#' @param x A matrix or data frame.
#' @param m The number of rows to view. Default: 6.
#' @param n The number of columns to view. Default: 6.
#' @return A matrix of size \code{m} x \code{n}.
#' @export
#' @examples
#' mat <- matrix(1:900, 30, 30)
#' headl(mat)
#' tailr(mat)
#' headl(mat, 10, 5)
#' headl(as.data.frame(mat))
headl <- function(x, m = 6L, n = 6L) {x[1:min(m, nrow(x)), 1:min(n, ncol(x))]}

#' @rdname headl
#' @export
headr <- function(x, m = 6L, n = 6L) {x[1:min(m, nrow(x)), (ncol(x) - min(n, ncol(x)) + 1):ncol(x)]}

#' @rdname headl
#' @export
taill <- function(x, m = 6L, n = 6L) {x[(nrow(x) - min(m, nrow(x)) + 1):nrow(x), 1:min(n, ncol(x))]}

#' @rdname headl
#' @export
tailr <- function(x, m = 6L, n = 6L) {x[(nrow(x) - min(m, nrow(x)) + 1):nrow(x), (ncol(x) - min(n, ncol(x)) + 1):ncol(x)]}





#' lapply that uses names
#'
#' A wrapper for \code{sapply} with arguments \code{simplify = FALSE} and \code{USE.NAMES = TRUE}.
#' @param x A vector (atomic or list).
#' @param FUN The function to be applied to each element of \code{x}.
#' @param ... Optional arguments to \code{FUN}.
#' @return A list with names equal to the names of \code{x}.
#' @export
slapply <- function(X, FUN, ...) {sapply(X, FUN, ..., simplify = FALSE, USE.NAMES = TRUE)}





#' ggplot2 default colours
#'
#' Generates colours like the default colours used in ggplot2, by taking points evenly distributed around a colour wheel.
#' @param n The number of colours to be returned.
#' @return A character vector of hexadecimal colour codes.
#' @export
ggplot_colours <- function(n) {hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:n]}





#' Make a blank plot with ggplot2
#'
#' Generates a blank plot, useful for making white space in combined plots using e.g. cowplot.
#' @return A \code{ggplot} object.
#' @export
blank_plot <- function() {ggplot() + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), 'pt'))}
