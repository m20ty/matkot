#' No-frills heatmap
#'
#' Makes a heatmap with optional row/column ordering and aligned dendrograms.
#' @param mat A numeric matrix.
#' @param ord An optional ordering object with which to reorder both the rows and the columns of \code{mat}. Can be an \code{hclust} object or an integer permutation vector.
#' @param ord_x An optional ordering object with which to reorder the columns of \code{mat}. Can be an \code{hclust} object or an integer permutation vector. Ignored if \code{ord} is supplied.
#' @param ord_y An optional ordering object with which to reorder the rows of \code{mat}. Can be an \code{hclust} object or an integer permutation vector. Ignored if \code{ord} is supplied.
#' @param draw_dendro Logical or character indicating whether to draw dendrograms or which dendrograms to draw. Possible values are \code{TRUE}, \code{FALSE}, 'x', 'y', 'xy'. If \code{TRUE} or 'xy', all possible dendrograms will be drawn. Default: \code{FALSE}.
#' @param return_data Logical indicating whether to return a data frame (in long format) containing the data used to construct the heatmap.
#' @param show_legend Logical indicating whether to show the heatmap legend.
#' @param plot_title The title of the plot. Set to \code{NULL} if you don't want a title. Default: \code{NULL}.
#' @param axis_text_x Character vector of x axis labels. These will be re-ordered if \code{ord} or \code{ord_x} are supplied. Set to \code{NULL} if you want no x axis labels. Default: \code{colnames(mat)}.
#' @param axis_text_y Character vector of y axis labels. These will be re-ordered if \code{ord} or \code{ord_y} are supplied. Set to \code{NULL} if you want no y axis labels. Default: \code{rownames(mat)}.
#' @param axis_title_x The x axis title. Set to \code{NULL} if you want no x axis title.
#' @param axis_title_y The y axis title. Set to \code{NULL} if you want no y axis title.
#' @param colours The colours to be used in the heatmap. These will be converted to a colour scale by ggplot2's \code{scale_fill_gradientn}. Default: c('red', 'yellow').
#' @param limits Two-element numeric vector specifying the limits to use for the colour scale.
#' @param oob A function to handle values lying outside the range specified by \code{limits}. Default: \code{scales::squish}.
#' @param legend_breaks A numeric vector of legend breaks.
#' @param legend_labels A vector of legend labels.
#' @param legend_title The legend title. Set to \code{NULL} if you want no legend title.
#' @param axis_ticks As \code{axis.ticks} in \code{theme}.
#' @param axis_ticks_length As \code{axis.ticks.length} in \code{theme}.
#' @param axis_text_size The size of the axis text. Default: 7.
#' @param plot_margin Numeric vector of plot margins in pts, in the order top, right, bottom, left. Default: c(5.5, 5.5, 5.5, 5.5).
#' @param rel_widths Named vector of relative widths to be passed to cowplot's \code{plot_grid}. Names should include 'heatmap' and, if applicable, 'dend' and 'legend'.
#' @param rel_heights Named vector of relative heights to be passed to cowplot's \code{plot_grid}. Names should include 'heatmap' and, if applicable, 'dend'.
#' @param ... Additional arguments to be passed to ggplot2's \code{theme} function in the construction of the heatmap.
#' @return If \code{return_data} is \code{FALSE}, a \code{ggplot} object. If \code{return_data} is \code{TRUE}, a list with two elements: \code{plot}, a \code{ggplot} object; and \code{data}, a data frame.
#' @export

heat_map <- function(

    mat,
    ord = NULL,
    ord_x = NULL,
    ord_y = NULL,
    draw_dendro = FALSE,
    return_data = FALSE,
    show_legend = TRUE,
    plot_title = NULL,
    axis_text_x = colnames(mat),
    axis_text_y = rownames(mat),
    axis_title_x = waiver(),
    axis_title_y = waiver(),
    colours = c('red', 'yellow'),
    limits = NULL,
    oob = scales::squish,
    legend_breaks = waiver(),
    legend_labels = waiver(),
    legend_title = waiver(),
    axis_ticks = element_blank(),
    axis_ticks_length = 0,
    axis_text_size = 7,
    plot_margin = c(5.5, 5.5, 5.5, 5.5),
    rel_widths = c(heatmap = 1, dend = 0.2, legend = 0.2),
    rel_heights = c(heatmap = 1, dend = 0.2),
    ...

) {

    # I would like to change this so that we don't have to specify the relative space for the legend.  This should be possible by converting the
    # heatmap to a grob and extracting the legend width (seems to be in cm), and then the rest of the plot should be unit(1, 'null').  But I cannot
    # figure out how to make this work.

    # We could also remove some arguments, leaving them to the '...' argument, even when a value is set here.  E.g. I included axis_ticks because I
    # set axis.ticks in theme() and would like the user to be able to override that, but I could just check if 'axis.ticks' appears in '...', and only
    # if it doesn't set the value myself.

    # <ord>, <ord_x> and <ord_y> can be integer permutation vectors or 'hclust' objects.  In the latter case, the labels and order
    # are taken from the hclust objects (ignoring the row and column names of mat), and the parameter <draw_dendro> enables construction of
    # dendrograms from the hclust objects, which will be aligned and plotted with the heatmap.  If <draw_dendro> is set to TRUE, as many dendrograms
    # as possible will be drawn.  The user can explicitly choose which dendrograms to draw by setting <draw_dendro> to 'x', 'y' or 'xy' ('xy' has the
    # same effect as TRUE).

    # The '...' parameter is for extra arguments to the call theme() in the construction of the heatmap.

    if(!(draw_dendro %in% c(FALSE, TRUE, 'x', 'y', 'xy'))) {stop("<draw_dendro> must be one of FALSE, TRUE, 'x', 'y' or 'xy'.")}
    if(!(show_legend %in% c(TRUE, FALSE, 'top right'))) {stop("<show_legend> must be one of FALSE, TRUE or 'top right'.")}
    # draw_dendro <- match.arg(draw_dendro, choices = c(FALSE, TRUE, 'x', 'y', 'xy'))
    # show_legend <- match.arg(show_legend, choices = c(TRUE, FALSE, 'top right'))
    if(show_legend == 'top right' & !(draw_dendro %in% c(TRUE, 'xy'))) {
        warning("<show_legend> set to 'top right' but fewer than 2 dendrograms requested.  Plotting legend to the right of the heatmap.")
        show_legend <- TRUE
    }

    if(is.null(rownames(mat))) {rownames(mat) <- as.character(1:nrow(mat))}
    if(is.null(colnames(mat))) {colnames(mat) <- as.character(1:ncol(mat))}

    plot_data <- data.table::melt(as.data.table(mat, keep.rownames = 'idy'), id.vars = 'idy', variable.name = 'idx', value.name = 'value')
    setcolorder(plot_data, c('idx', 'idy', 'value'))

    if(!is.null(ord)) { # If value is supplied to <ord>, use the same ord for both axes
        if(dim(mat)[1] != dim(mat)[1]) { # Check that matrix is square
            stop('A non-NULL value has been supplied for <ord>, but <mat> is not square.  Use <ord_x> and/or <ord_y> instead.')
        } else {
            if('hclust' %in% class(ord)) {
                if(!all(ord$labels == rownames(mat)) || !all(ord$labels == colnames(mat))) {
                    stop('Labels from <ord> do not match row/column names of <mat>.')
                }
                plot_data[, c('idx', 'idy') := .(factor(idx, levels = with(ord, labels[order])), factor(idy, levels = with(ord, labels[order])))]
            } else {
                plot_data[, c('idx', 'idy') := .(factor(idx, levels = colnames(mat)[ord]), factor(idy, levels = rownames(mat)[ord]))]
            }
        }
    } else {
        if(!is.null(ord_x)) {
            if('hclust' %in% class(ord_x)) {
                if(!all(ord_x$labels == colnames(mat))) {stop('Labels from <ord_x> do not match column names of <mat>.')}
                plot_data[, idx := factor(idx, levels = with(ord_x, labels[order]))]
            } else {
                plot_data[, idx := factor(idx, levels = colnames(mat)[ord_x])]
            }
        }
        if(!is.null(ord_y)) {
            if('hclust' %in% class(ord_y)) {
                if(!all(ord_y$labels == rownames(mat))) {stop('Labels from <ord_y> do not match row names of <mat>.')}
                plot_data[, idy := factor(idy, levels = with(ord_y, labels[order]))]
            } else {
                plot_data[, idy := factor(idy, levels = rownames(mat)[ord_y])]
            }
        }
    }

    # First construct dendrograms, so we know what to align and how to set the heatmap plot margins.

    dend_list <- NULL

    if(draw_dendro %in% c(TRUE, 'xy')) {
        if(!is.null(ord)) {
            if(!('hclust' %in% class(ord))) {
                warning("Cannot draw dendrograms when <ord> does not have class 'hclust'.")
            } else {
                dend_list <- list(
                    x = dendro(ord, edge = 'bottom', plot_margin = c(plot_margin[1:2], 0, plot_margin[4])),
                    y = dendro(ord, edge = 'left', plot_margin = c(plot_margin[1:3], 0))
                )
            }
        } else {
            if((is.null(ord_x) || !('hclust' %in% class(ord_x))) & (is.null(ord_y) || !('hclust' %in% class(ord_y)))) {
                warning("Cannot draw dendrograms when <ord_x> and <ord_y> are NULL or do not have class 'hclust'.")
            } else {
                if(!is.null(ord_x)) {
                    if(!('hclust' %in% class(ord_x))) {
                        warning("Cannot draw x dendrogram when <ord_x> does not have class 'hclust'.")
                    } else {
                        dend_list <- c(dend_list, list(x = dendro(ord_x, edge = 'bottom', plot_margin = c(plot_margin[1:2], 0, plot_margin[4]))))
                    }
                }
                if(!is.null(ord_y)) {
                    if(!('hclust' %in% class(ord_y))) {
                        warning("Cannot draw y dendrogram when <ord_y> does not have class 'hclust'.")
                    } else {
                        dend_list <- c(dend_list, list(y = dendro(ord_y, edge = 'left', plot_margin = c(plot_margin[1:3], 0))))
                    }
                }
            }
        }
    } else if(draw_dendro == 'x') {
        if((is.null(ord) || !('hclust' %in% class(ord))) & (is.null(ord_x) || !('hclust' %in% class(ord_x)))) {
            warning("Cannot draw x dendrogram when <ord> and <ord_x> are NULL or do not have class 'hclust'.")
        } else if(!is.null(ord)) {
            if(!('hclust' %in% class(ord))) {
                warning("Cannot draw x dendrogram when <ord> does not have class 'hclust'.")
            } else {
                dend_list <- list(x = dendro(ord, edge = 'bottom', plot_margin = c(plot_margin[1:2], 0, plot_margin[4])))
            }
        } else if(!is.null(ord_x)) {
            if(!('hclust' %in% class(ord_x))) {
                warning("Cannot draw x dendrogram when <ord_x> does not have class 'hclust'.")
            } else {
                dend_list <- list(x = dendro(ord_x, edge = 'bottom', plot_margin = c(plot_margin[1:2], 0, plot_margin[4])))
            }
        }
    } else if(draw_dendro == 'y') {
        if((is.null(ord) || !('hclust' %in% class(ord))) & (is.null(ord_y) || !('hclust' %in% class(ord_y)))) {
            warning("Cannot draw y dendrogram when <ord> and <ord_y> are NULL or do not have class 'hclust'.")
        } else if(!is.null(ord)) {
            if(!('hclust' %in% class(ord))) {
                warning("Cannot draw y dendrogram when <ord> does not have class 'hclust'.")
            } else {
                dend_list <- list(y = dendro(ord, edge = 'left', plot_margin = c(plot_margin[1:3], 0)))
            }
        } else if(!is.null(ord_y)) {
            if(!('hclust' %in% class(ord_y))) {
                warning("Cannot draw y dendrogram when <ord_y> does not have class 'hclust'.")
            } else {
                dend_list <- list(y = dendro(ord_y, edge = 'left', plot_margin = c(plot_margin[1:3], 0)))
            }
        }
    }

    htmp <- ggplot(plot_data, aes(x = idx, y = idy, fill = value)) +
        geom_raster() +
        scale_x_discrete(labels = axis_text_x, expand = c(0, 0)) + # Change "breaks" to "labels"?
        scale_y_discrete(labels = axis_text_y, expand = c(0, 0)) +
        scale_fill_gradientn(colours = colours, limits = limits, oob = oob, breaks = legend_breaks, labels = legend_labels) +
        theme(
            axis.ticks = axis_ticks,
            axis.ticks.length = unit(axis_ticks_length, 'pt'),
            axis.text = element_text(size = axis_text_size),
            axis.text.y = element_text(vjust = 0.5), # vjust makes labels in centre of rows
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # likewise columns
            plot.margin = unit(
                c(
                    switch(('x' %in% names(dend_list)) + 1, plot_margin[1], 0),
                    switch(('y' %in% names(dend_list)) + 1, plot_margin[2], 0),
                    plot_margin[3:4]
                ),
                'pt'
            ),
            panel.background = element_blank(), # Eliminates the grey from the background (sometimes the grey is still visible without this)
            ...
        ) +
        labs(fill = legend_title, x = axis_title_x, y = axis_title_y, title = plot_title)

    # if(paste(axis_text_y, collapse = '') == '') {htmp <- htmp + theme(axis.text.y = element_blank())}
    # if(paste(axis_text_x, collapse = '') == '') {htmp <- htmp + theme(axis.text.x = element_blank())}

    if(!is.null(dend_list)) {
        if(length(dend_list) == 2) {
            aligned_plots_1 <- align_plots(
                dend_list$x + theme(plot.margin = unit(c(plot_margin[1], 0, 0, plot_margin[4]), 'pt')),
                htmp + theme(legend.position = 'none'),
                align = 'v'
            )
            aligned_plots_2 <- align_plots(
                htmp + theme(legend.position = 'none'),
                dend_list$y + theme(plot.margin = unit(c(0, plot_margin[2:3], 0), 'pt')),
                align = 'h'
            )
            if(show_legend == TRUE) {
                out <- list(
                    combined_plot = plot_grid(
                        plot_grid(
                            aligned_plots_1[[1]],
                            blank_plot(),
                            aligned_plots_2[[1]],
                            aligned_plots_2[[2]],
                            nrow = 2,
                            ncol = 2,
                            rel_widths = rel_widths[c('heatmap', 'dend')],
                            rel_heights = rel_heights[c('dend', 'heatmap')]
                        ),
                        get_legend(htmp),
                        nrow = 1,
                        ncol = 2,
                        rel_widths = c(sum(rel_widths[c('heatmap', 'dend')]), rel_widths['legend'])
                    )
                )
            } else if(show_legend == 'top right') {
                out <- list(
                    combined_plot = plot_grid(
                        aligned_plots_1[[1]],
                        get_legend(htmp),
                        aligned_plots_2[[1]],
                        aligned_plots_2[[2]],
                        nrow = 2,
                        ncol = 2,
                        rel_widths = rel_widths[c('heatmap', 'dend')],
                        rel_heights = rel_heights[c('dend', 'heatmap')]
                    )
                )
            } else {
                out <- list(
                    combined_plot = plot_grid(
                        aligned_plots_1[[1]],
                        blank_plot(),
                        aligned_plots_2[[1]],
                        aligned_plots_2[[2]],
                        nrow = 2,
                        ncol = 2,
                        rel_widths = rel_widths[c('heatmap', 'dend')],
                        rel_heights = rel_heights[c('dend', 'heatmap')]
                    )
                )
            }
            if(return_data) {
                return(c(out, list(plots = list(heatmap = htmp, dendx = dend_list$x, dendy = dend_list$y), data = plot_data)))
            } else {
                return(out$combined_plot)
            }
        } else if(names(dend_list) == 'x') {
            aligned_plots <- align_plots(dend_list$x, htmp + theme(legend.position = 'none'), align = 'v')
            if(show_legend == TRUE | show_legend == 'top right') {
                if(show_legend == 'top right') {warning("Setting <show_legend> to 'top right' is only valid when two dendrograms are drawn.")}
                out <- list(
                    combined_plot = plot_grid(
                        plot_grid(
                            aligned_plots[[1]],
                            aligned_plots[[2]],
                            nrow = 2,
                            ncol = 1,
                            rel_heights = rel_heights[c('dend', 'heatmap')]
                        ),
                        get_legend(htmp),
                        nrow = 1,
                        ncol = 2,
                        rel_widths = c(sum(rel_widths[c('heatmap', 'dend')]), rel_widths['legend'])
                    )
                )
            } else {
                out <- list(
                    combined_plot = plot_grid(
                        aligned_plots[[1]],
                        aligned_plots[[2]],
                        nrow = 2,
                        ncol = 1,
                        rel_heights = rel_heights[c('dend', 'heatmap')]
                    )
                )
            }
            if(return_data) {
                return(c(out, list(plots = list(heatmap = htmp, dend = dend_list$x), data = plot_data)))
            } else {
                return(out$combined_plot)
            }
        } else if(names(dend_list) == 'y') {
            aligned_plots <- align_plots(htmp + theme(legend.position = 'none'), dend_list$y, align = 'h')
            if(show_legend == TRUE | show_legend == 'top right') {
                if(show_legend == 'top right') {warning("Setting <show_legend> to 'top right' is only valid when two dendrograms are drawn.")}
                out <- list(
                    combined_plot = plot_grid(
                        aligned_plots[[1]],
                        aligned_plots[[2]],
                        get_legend(htmp),
                        nrow = 1,
                        ncol = 3,
                        rel_widths = rel_widths[c('heatmap', 'dend', 'legend')]
                    )
                )
            } else {
                out <- list(
                    combined_plot = plot_grid(
                        aligned_plots[[1]],
                        aligned_plots[[2]],
                        nrow = 1,
                        ncol = 2,
                        rel_widths = rel_widths[c('heatmap', 'dend')]
                    )
                )
            }
            if(return_data) {
                return(c(out, list(plots = list(heatmap = htmp, dend = dend_list$y), data = plot_data)))
            } else {
                return(out$combined_plot)
            }
        }
    } else {
        if(return_data) {
            return(list(plot = htmp, data = plot_data))
        } else {
            return(htmp)
        }
    }

}
