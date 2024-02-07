#' Find Optimal Permuted Matrix.
#'
#' Find \eqn{A_{perm}}, the optimal permuted matrix of \eqn{A}, such that
#' \eqn{|A_{perm} - B|_F} is minimized.
#'
#' @param A a matrix
#' @param B a matrix of same dimensions as \eqn{A}
#' @return the optimal permutation of \eqn{A}
#' @export

permA <- function(A, B) {
  if (sum(dim(A) == dim(B)) == 2) {
    pureInd <- pureRowInd(B)
    optPerm <- getOptPerm(A[pureInd, ], B[pureInd, ])
    permutatedA <- A %*% optPerm[[1]]
    A <- t(t(permutatedA) * optPerm[[2]])
  }
  return(A)
}
