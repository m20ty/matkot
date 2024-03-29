% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse_matrix_functions.R
\name{col_nnz}
\alias{col_nnz}
\alias{row_nnz}
\alias{nnz}
\title{Number of zeros}
\usage{
col_nnz(mat)

row_nnz(mat)

nnz(mat)
}
\arguments{
\item{mat}{A matrix. May be sparse (i.e. of class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}).}
}
\value{
An integer or named integer vector.
}
\description{
Find the number of non-zero elements in a matrix per row, per column or in total. In the case of sparse matrices (class \code{dgTMatrix}, \code{dgCMatrix} or \code{dgRMatrix}), this is much faster than using \code{apply} and/or \code{sum}. Note there already exists a function \code{nnzero} in the Matrix package that finds the number of non-zero elements in a sparse matrix, but \code{nnz} is faster.
}
\examples{
mat <- matrix(rbinom(900, 1, 0.1), 30, 30)
col_nnz(mat)
rownames(mat) <- replicate(30, paste(sample(letters, 10), collapse = ''))
row_nnz(mat)
nnz(mat)
mat <- as(mat, 'dgCMatrix')
row_nnz(mat)
nnz(mat)
}
