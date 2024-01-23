#' Estimate \eqn{\Omega}.
#'
#' For a given \eqn{\lambda} and \eqn{C}, find \eqn{C^{-1}}.
#'
#' @param lambda \eqn{\lambda}
#' @param C a square, symmetric matrix
#' @return \eqn{\Omega = C^{-1}}
#' @export

estOmega <- function(lambda, C) {
  K <- nrow(C)
  omega <- matrix(0, K, K)
  for (i in 1:K) {
    omega[, i] <- solveRow(col_ind = i, C = C, lambda = lambda)
  }
  return(omega)
}
