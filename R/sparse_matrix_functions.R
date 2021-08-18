#' Number of zeros
#'
#' Find the number of non-zero elements in a matrix per row, per column or in total. In the case of sparse matrices (class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}), this is much faster than using \code{apply} and/or \code{sum}. Note there already exists a function \code{nnzero} in the Matrix package that finds the number of non-zero elements in a sparse matrix, but \code{nnz} is faster.
#' @param mat A matrix. May be sparse (i.e. of class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}).
#' @return An integer or named integer vector.
#' @export
#' @examples
#' mat <- matrix(rbinom(900, 1, 0.1), 30, 30)
#' col_nnz(mat)
#' rownames(mat) <- replicate(30, paste(sample(letters, 10), collapse = ''))
#' row_nnz(mat)
#' nnz(mat)
#' mat <- as(mat, 'dgCMatrix')
#' row_nnz(mat)
#' nnz(mat)
col_nnz <- function(mat) {
    if(any(c('dgTMatrix', 'dgRMatrix') %in% class(mat))) {
        setNames(tabulate(mat@j + 1, nbins = mat@Dim[2]), colnames(mat))
    } else if('dgCMatrix' %in% class(mat)) {
        setNames(diff(mat@p), colnames(mat))
    } else {
        apply(mat, 2, function(x) sum(x != 0))
    }
}

#' @rdname col_nnz
#' @export
row_nnz <- function(mat) {
    if(any(c('dgTMatrix', 'dgCMatrix') %in% class(mat))) {
        setNames(tabulate(mat@i + 1, nbins = mat@Dim[1]), rownames(mat))
    } else if('dgRMatrix' %in% class(mat)) {
        setNames(diff(mat@p), rownames(mat))
    } else {
        apply(mat, 1, function(x) sum(x != 0))
    }
}

#' @rdname col_nnz
#' @export
nnz <- function(mat) {
    if(any(c('dgTMatrix', 'dgCMatrix', 'dgRMatrix') %in% class(mat))) {
        length(mat@x)
    } else {
        sum(mat != 0)
    }
}





#' Convert to fractions of row or column totals
#'
#' Divides the elements of a matrix by their corresponding row or column totals. This is useful, for example, when converting a gene expression count matrix to transcripts per million (TPM).
#' @param mat A matrix. May be sparse (i.e. of class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}).
#' @param MARGIN The subscript to apply over. Possible values are \code{1} and \code{2}. If \code{1}, each element will be divided by its corresponding row total; if \code{2}, its column total.
#' @return A matrix of the same size and class as \code{mat}.
#' @export
#' @examples
#' mat <- matrix(round(runif(100, max = 100)), 10, 10)
#' to_frac(mat)
#' to_frac(mat, MARGIN = 1)
#' mat[sample(1:100, 80)] <- 0
#' mat <- as(mat, 'dgCMatrix')
#' to_frac(mat)

to_frac <- function(mat, MARGIN = 2) {

    MARGIN <- match.arg(as.character(MARGIN), c('1', '2'))

    if(MARGIN == '1') {

        if(any(c('dgTMatrix', 'dgCMatrix') %in% class(mat))) {
            out <- mat
            # In the following, we add names to mat@x so that we can preserve the order of x:
            out@x <- unname(
                unlist(
                    unname(tapply(setNames(mat@x, as.character(1:length(mat@x))), mat@i, function(y) y/sum(y), simplify = FALSE))
                )[as.character(1:length(mat@x))]
            )
            out
        } else if(class(mat) == 'dgRMatrix') { # Then mat doesn't have an i attribute
            out <- mat
            # We don't need to use names in the following, because p is ordered with respect to x.
            # Will there be problems if there are all-zero rows at the bottom of the matrix?
            out@x <- unlist(unname(tapply(mat@x, cut(1:length(mat@x), mat@p), function(y) y/sum(y), simplify = FALSE)))
            out
        } else {
            apply(mat, 2, function(x) x/sum(x))
        }

    } else { # Then MARGIN == '2'

        if(any(c('dgTMatrix', 'dgRMatrix') %in% class(mat))) {
            out <- mat
            out@x <- unname(
                unlist(
                    unname(tapply(setNames(mat@x, as.character(1:length(mat@x))), mat@j, function(y) y/sum(y), simplify = FALSE))
                )[as.character(1:length(mat@x))]
            )
            out
        } else if(class(mat) == 'dgCMatrix') { # Then mat doesn't have a j attribute
            out <- mat
            # Will there be problems in the following if there are all-zero columns at the right of the matrix?
            out@x <- unlist(unname(tapply(mat@x, cut(1:length(mat@x), mat@p), function(y) y/sum(y), simplify = FALSE)))
            out
        } else {
            t(apply(mat, 1, function(x) x/sum(x)))
        }

    }

}





# NOTE ON THE FOLLOWING FUNCTION: There are log1p() and expm1() functions that do this and seem to preserve sparse matrix class. Check if this is so
# and check which is faster. I think these two functions do not allow choosing a base (so base is always e), but it would be good to mention them in
# the description if they work.

#' Log transformation of a matrix
#'
#' Log transforms a matrix after adding 1, or the reverse, exponentiates the matrix and then subtracts 1. In the case of sparse matrices (class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}), only the non-zero elements are transformed. This is faster than manually adding or subtracting 1, as it retains the sparse matrix class.
#' @param mat A matrix. May be sparse (i.e. of class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}).
#' @param base A positive number: the base with respect to which the logarithms are computed. Default: 2.
#' @param reverse If \code{FALSE}, log transform \code{mat}; if \code{TRUE}, reverse log transform, i.e. exponentiate \code{mat}.
#' @return A matrix of the same size and class as \code{mat}.
#' @export
#' @examples
#' mat <- matrix(round(runif(100, max = 100)), 10, 10)
#' log_transform(mat)
#' log_transform(log_transform(mat), reverse = TRUE)
#' mat[sample(1:100, 80)] <- 0
#' mat <- as(mat, 'dgCMatrix')
#' log_transform(mat)

log_transform <- function(mat, base = 2, reverse = FALSE) {
    if(any(c('dgTMatrix', 'dgCMatrix', 'dgRMatrix') %in% class(mat))) {
        out <- mat
        if(reverse) {
            out@x <- base^out@x - 1
        } else {
            out@x <- log(out@x + 1, base = base)
        }
        return(out)
    } else {
        if(reverse) {
            return(base^mat - 1)
        } else {
            return(log(mat + 1, base = base))
        }
    }
}





#' Shuffle matrix elements by row or column
#'
#' Re-orders the elements in each row or column of the matrix. Useful for defining null distributions. In the case of large sparse matrices (class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}), this is faster than using \code{apply}, although for small matrices it may be slower.
#' @param mat A matrix. May be sparse (i.e. of class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}).
#' @param MARGIN The subscript to apply over. Possible values are \code{1} and \code{2}. If \code{1}, the rows of the matrix are re-ordered; if \code{2}, the columns.
#' @return A matrix of the same size and class as \code{mat}.
#' @export
#' @examples
#' mat <- matrix(round(runif(100, max = 100)), 10, 10)
#' shuffle(mat)
#' shuffle(mat, MARGIN = 1)
#' mat[sample(1:100, 80)] <- 0
#' mat <- as(mat, 'dgCMatrix')
#' shuffle(mat)

shuffle <- function(mat, MARGIN = 2) {

    # To do: add an option for MARGIN = NULL, which would mean we just shuffle all values randomly irrespective of rows or columns.

    MARGIN <- match.arg(as.character(MARGIN), c('1', '2'))

    if('dgCMatrix' %in% class(mat)) {
        out <- mat
        if(MARGIN == '1') {
            out_attrs <- data.table(ind = 1:length(mat@i), rn = mat@i, val = mat@x)[,
                .(ind = ind, val = switch((.N == 1) + 1, sample(val), val), cn = sample(1L:ncol(mat) - 1L, .N)),
                by = rn
            ][order(cn, rn)]
            out@i <- out_attrs$rn
            out@p <- out_attrs[, c(0L, cumsum(tabulate(cn + 1L)))]
            out@x <- out_attrs$val
        } else {
            out_attrs <- data.table(ind = 1:ncol(mat), p0 = mat@p[-length(mat@p)] + 1, p1 = mat@p[-1], pdiff = diff(mat@p))[
                p0 <= p1,
                .(rn = sort(sample(1L:nrow(mat) - 1L, pdiff)), val = switch((p0 == p1) + 1, sample(mat@x[p0:p1]), mat@x[p0])),
                by = ind
            ]
            out@i <- out_attrs$rn
            out@x <- out_attrs$val
        }
        return(out)
    } else if('dgRMatrix' %in% class(mat)) {
        out <- mat
        if(MARGIN == '1') {
            out_attrs <- data.table(ind = 1:nrow(mat), p0 = mat@p[-length(mat@p)] + 1, p1 = mat@p[-1], pdiff = diff(mat@p))[
                p0 <= p1,
                .(cn = sort(sample(1L:ncol(mat) - 1L, pdiff)), val = switch((p0 == p1) + 1, sample(mat@x[p0:p1]), mat@x[p0])),
                by = ind
            ]
            out@j <- out_attrs$cn
            out@x <- out_attrs$val
        } else {
            out_attrs <- data.table(ind = 1:length(mat@j), cn = mat@j, val = mat@x)[,
                .(ind = ind, val = switch((.N == 1) + 1, sample(val), val), rn = sample(1L:nrow(mat) - 1L, .N)),
                by = cn
            ][order(rn, cn)]
            out@j <- out_attrs$cn
            out@p <- out_attrs[, c(0L, cumsum(tabulate(rn + 1L)))]
            out@x <- out_attrs$val
        }
        return(out)
    } else if('dgTMatrix' %in% class(mat)) {
        out <- mat
        if(MARGIN == '1') {
            out@j <- data.table(ind = 1:length(mat@i), rn = mat@i)[, .(ind = ind, cn = sample(1L:ncol(mat) - 1L, .N)), by = rn][order(ind), cn]
        } else {
            out@i <- data.table(ind = 1:length(mat@j), cn = mat@j)[, .(ind = ind, rn = sample(1L:nrow(mat) - 1L, .N)), by = cn][order(ind), rn]
        }
        return(out)
    } else {
        if(MARGIN == '1') {
            return(set_colnames(t(apply(mat, 1, sample)), colnames(mat)))
        } else {
            return(set_rownames(apply(mat, 2, sample), rownames(mat)))
        }
        # and if margin is NULL: matrix(sample(mat), nrow(mat), ncol(mat)), then set rownames and colnames to same as mat
    }

}
