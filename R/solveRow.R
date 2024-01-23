#' ???
#'
#' ???
#'
#' @param col_ind indices of columns
#' @param C matrix?
#' @param lambda \eqn{\labmda}
#' @return value
#' @export

solveRow <- function(col_ind, C, lambda) {
  #### set up linear program
  K <- nrow(C)
  c_vec <- c(1, rep(0, 2 * K))
  A_mat <- -c_vec
  A_mat <- rbind(A_mat, c(-1, rep(1, 2 * K)))
  tmp_constr <- C %x% t(c(1, -1))
  A_mat <- rbind(A_mat, cbind(-1 * lambda, rbind(tmp_constr, -tmp_constr)))
  tmp_vec <- rep(0, K)
  tmp_vec[col_ind] <- 1
  b_vec <- c(0, 0, tmp_vec, -tmp_vec)

  ## solving (cvec)'(x) subject to Amat %*% x =  bvec and x >= 0
  lp_result <- linprog::solveLP(c_vec, b_vec, A_mat, lpSolve = T)$solution
  while (length(lp_result) == 0) {
    cat("The penalty lambda =", lambda, "is too small and increased by 0.01...\n")
    lambda <- lambda + 0.01
    A_mat[-(1:2), 1] <- lambda
    lp_result <- linprog::solveLP(c_vec, b_vec, A_mat, lpSolve = T)$solution[-1]
  }
  ind <- seq(2, 2 * K, 2)
  return(lp_result[ind] - lp_result[ind + 1])
}
