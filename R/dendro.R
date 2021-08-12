#' Make a dendrogram from an hclust object
#'
#' Makes a dendrogram aligned to a chosen edge of the plot area. This function is useful for making a dendrogram to align next to a heatmap. It is used in the \code{heat_map} function.
#' @param clust An \code{hclust} object.
#' @param edge The edge of the plot window to which to align the dendrogram, possible values being \code{'bottom'}, \code{'left'}, \code{'top'} and \code{'right'}. For example, if the dendrogram is to be plotted above a heatmap, use \code{edge = 'bottom'}. Default: \code{'bottom'}.
#' @param plot_margin A numeric vector of margin sizes in pts, in the order top, right, bottom, left. Default: 0 for \code{edge}, 5.5 for all other margins.
#' @param size The weight of the lines.
#' @return A ggplot object.
#' @export

dendro <- function(
    clust,
    edge = c('bottom', 'left', 'top', 'right'),
    plot_margin = sapply(c('top', 'right', 'bottom', 'left'), function(x) ifelse(x == edge, 0, 5.5)),
    size = 0.5
) {

    edge <- match.arg(edge)

    dend_data <- ggdendro::dendro_data(as.dendrogram(clust))

    if(edge == 'bottom') {

        ggplot(dend_data$segments) +
            geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = size) +
            scale_x_continuous(expand = c(0, 0.5)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )

    } else if(edge == 'left') {

        ggplot(dend_data$segments) +
            geom_segment(aes(x = y, y = x, xend = yend, yend = xend), size = size) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0.5)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )

    } else if(edge == 'top') {

        ggplot(dend_data$segments) +
            geom_segment(aes(x = x, y = -y, xend = xend, yend = -yend), size = size) +
            scale_x_continuous(expand = c(0, 0.5)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )

    } else if(edge == 'right') {

        ggplot(dend_data$segments) +
            geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend), size = size) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0.5)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )

    }

}
