#' Find \eqn{A_I} From Pure Nodes.
#'
#' Recover the estimated submatrix \eqn{A_I} by given the pure node group.
#'
#' @param pure_list a list of group indices of the pure nodes with sign
#' @param p the number of features
#' @return a matrix of dimensions \eqn{p \times K}
#' @export

recoverAI <- function(pure_list, p) {
  K <- length(pure_list)
  A <- matrix(0, p, K)
  for (i in 1:K) {
    group_i <- pure_list[[i]]
    A[group_i[[1]], i] <- 1
    group_i2 <- group_i[[2]]
    if (length(group_i2) != 0) {
      A[group_i2, i] <- -1
    }
  }
  return(A)
}
