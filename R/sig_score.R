#' Compute gene signature scores
#'
#' Compute gene signature scores for a matrix, returning one score per column. The scoring procedure is that used in Tirosh et al. 2016 (https://www.nature.com/articles/nature20123). This function is a simplified version of the \code{sigScores} function from the scalop package (https://github.com/jlaffy/scalop).
#' @param mat A matrix with genes for rows and samples for columns.
#' @param sig A character vector of signature genes.
#' @param nbin The number of expression bins into which to divide the genes.
#' @param n The number of control genes to sample for each gene. Should be less than or equal to \code{nrow(mat)/nbin}
#' @param replace Logical indicating whether to sample with replacement. Default: \code{FALSE}.
#' @param return_control_sets Logical indicating whether to return the control genes. Default: \code{FALSE}.
#' @return If \code{return_control_sets} is \code{FALSE}, a vector of scores with length equal to \code{ncol(mat)}. If \code{return_control_sets} is \code{TRUE}, a list with three elements: \code{scores}, a vector of scores with length equal to \code{ncol(mat)}; \code{controls}, a list of the control gene sets used for each gene in \code{sig}; and \code{comparable_gene_sets}, a list of length \code{n} and a rearrangement of \code{controls}, each element being a vector of genes with comparable expression levels to those in \code{sig} (this is redundant but possibly helpful).
#' @export

sig_score <- function(
    mat,
    sig,
    nbin = nrow(mat) %/% 110,
    n = 100,
    replace = FALSE,
    return_control_sets = FALSE
) {

    # We could use Matrix::rowMeans for either case, but it seems to be faster if we first choose the right function using an if statement.
    # if(any(c('dgTMatrix', 'dgCMatrix', 'dgRMatrix') %in% class(mat))) {
    #     gene_averages <- sort(Matrix::rowMeans(mat))
    # } else {
    #     gene_averages <- sort(base::rowMeans(mat))
    # }
    gene_averages <- sort(Matrix::rowMeans(mat))
    bins <- setNames(cut(seq_along(gene_averages), breaks = nbin, labels = FALSE, include.lowest = TRUE), names(gene_averages))

    if(return_control_sets) {

        # Define control gene sets for distribution of scores:
        sig_gene_controls <- sapply(
            sig,
            function(g) sample(names(bins)[bins == bins[g]], n, replace = replace),
            simplify = FALSE,
            USE.NAMES = TRUE
        )
        comparable_gene_sets <- lapply(1:n, function(i) sapply(sig_gene_controls, `[`, i))
        sig_scores <- rowMeans(sapply(sig, function(g) mat[g, ] - colMeans(mat[sig_gene_controls[[g]], ])))

        return(list(scores = sig_scores, controls = sig_gene_controls, comparable_gene_sets = comparable_gene_sets))

    } else {

        sig_scores <- rowMeans(
            sapply(sig, function(g) mat[g, ] - Matrix::colMeans(mat[sample(names(bins)[bins == bins[g]], n, replace = replace), ]))
        )

        return(sig_scores)

    }

}
