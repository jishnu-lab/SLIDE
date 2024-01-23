#' Estimate \eqn{Z}.
#'
#' Additional prediction function using the support of \eqn{\beta}.
#'
#' @param x data matrix, standardized
#' @param er_res the results of a run of \code{plainER()} or \code{priorER()}
#' @return estimates for \eqn{Z}
#' @export

predZ <- function(x, er_res) {
  A_hat <- er_res$A
  C_hat <- er_res$C
  Gamma_hat <- er_res$Gamma
  Gamma_hat <- ifelse(Gamma_hat == 0, 1e-10, Gamma_hat)
  Gamma_hat_inv <- diag(Gamma_hat ** (-1))
  G_hat <- crossprod(A_hat, Gamma_hat_inv) %*% A_hat + solve(C_hat)
  Z_hat <- x %*% Gamma_hat_inv %*% A_hat %*% MASS::ginv(G_hat)
  return (Z_hat)
}
