#' Dantzig Estimation.
#'
#' Estimate a vector using a Dantzig-type estimator.
#'
#' @param C_hat a matrix
#' @param target_vec a vector
#' @param lambda \eqn{\lambda}, a positive constant
#' @return a row of \eqn{A_J} estimated with a Dantzig-type estimator
#' @export

dantzig <- function(C_hat, target_vec, lambda) {
  K <- length(target_vec) ## number of clusters (columns in Y_hat)
  c_vec <- rep(1, 2 * K) ## just 1s because we want to minimize x itself
  b_vec <- c(lambda + target_vec, lambda - target_vec, rep(0, 2 * K))
  new_C_hat <- matrix(0, K, 2 * K)
  for (i in 1:K) {
    new_C_hat[i, ] <- c(C_hat[i, ], -C_hat[i, ]) ## two copies of C_hat, but second copy is opposite signed
  }
  A_mat <- rbind(new_C_hat, -new_C_hat, diag(-1, nrow = 2 * K))
  ## minimize c_vec'x subject to A_mat %*% x ≤ b_vec and x ≥ 0
  LP_sol <- linprog::solveLP(c_vec, b_vec, A_mat, lpSolve = T, maxiter = 1000000)$solution
  if (is.null(LP_sol)) {
    return (NULL)
  }
  AJ_row <- LP_sol[1:K] - LP_sol[(K + 1):(2 * K)] ## split up positive and negative portions of beta
  return (AJ_row)
}


