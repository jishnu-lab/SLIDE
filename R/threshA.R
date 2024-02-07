#' Threshold \eqn{A}.
#'
#' Threshold the given matrix \code{A} based on the given \eqn{\mu}. If \code{scale} is true,
#' then normalize each row of \code{A} such that the \eqn{L_1} norm of each row is not larger than 1.
#'
#' @param A an estimated matrix
#' @param mu threshold, a positive value
#' @param scale boolean for whether to normalize
#' @return thresholded version of matrix \code{A}
#' @export

threshA <- function(A, mu, scale = FALSE) {
  scaledA <- A
  for (i in 1:nrow(A)) {
    colInd <- abs(A[i, ]) <= mu
    scaledA[i,colInd] = 0
    if (scale && sum(abs(scaledA[i, ])) > 1) {
      scaledA[i, ] <- scaledA[i, ] / sum(abs(scaledA[i, ]))
    }
  }
  return(scaledA)
}
