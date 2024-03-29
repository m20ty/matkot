% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_matrix_functions.R
\name{shuffle}
\alias{shuffle}
\title{Shuffle matrix elements by row or column}
\usage{
shuffle(mat, MARGIN = 2)
}
\arguments{
\item{mat}{A matrix. May be sparse (i.e. of class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}).}

\item{MARGIN}{The subscript to apply over. Possible values are \code{1} and \code{2}. If \code{1}, the rows of the matrix are re-ordered; if \code{2}, the columns.}
}
\value{
A matrix of the same size and class as \code{mat}.
}
\description{
Re-orders the elements in each row or column of the matrix. Useful for defining null distributions. In the case of large sparse matrices (class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}), this is faster than using \code{apply}, although for small matrices it may be slower.
}
\examples{
mat <- matrix(round(runif(100, max = 100)), 10, 10)
shuffle(mat)
shuffle(mat, MARGIN = 1)
mat[sample(1:100, 80)] <- 0
mat <- as(mat, 'dgCMatrix')
shuffle(mat)
}
