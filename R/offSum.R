#' Off-Diagonal Sum Of Squares.
#'
#' Calculate the sum of squares of the upper off-diagonal elements of two matrices.
#'
#' @param mat1 a matrix
#' @param mat2 a matrix of same dimensions as \code{mat1}
#' @param weights the weights to use in calculating the sum
#' @return the sum of squares
#' @export

offSum <- function(mat1, mat2, weights) {
  tmp <- (mat1 - mat2) / weights
  tmp <- t(t(tmp) / weights)
  return(sum((tmp[row(tmp) <= (col(tmp) - 1)])^2))
}
