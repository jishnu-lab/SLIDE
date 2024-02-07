#' Calculate Norms.
#'
#' Calculate the scaled matrix \eqn{L_1} and Frobenius norms.
#'
#' @param A a matrix
#' @return a vector containing the scaled \eqn{L_1} and Frobenius norms
#' @export

calARates <- function(A) {
  frob <- norm(A, "F") / sqrt(ncol(A) * nrow(A))
  l1 <- sum(abs(A)) / nrow(A) / ncol(A)
  return(c(l1 = l1, l2 = frob))
}
