#' compASV.
#'
#' I am not sure what this does.
#'
#' @param sigma_alt a matrix computed by \code{\link{estSigmaAlt}}
#' @param BI estimated matrix:
#' \eqn{\hat{A}_{\hat{I}\cdot}(\hat{A}_{\hat{I}\cdot}^\top \hat{A}_{\hat{I}\cdot}^\top)^{-1}}
#' @param Theta_hat estimated matrix:
#' \eqn{(\hat{\Sigma}_{\cdot \hat{I}}-\hat{\Gamma}_{\cdot \hat{I}})BI}
#' @param Gamma_hat estimated matrix
#' @param beta_hat matrix of the estimated values of \eqn{\beta}
#' @param Omega_hat estimated matrix
#' @param I_hat vector of pure node indices
#' @param I_hat_list list of pure node indices by cluster
#' @return ???
#' @export

compASV <- function(sigma_alt, BI, Theta_hat, Gamma_hat, beta_hat, Omega_hat, I_hat, I_hat_list) {
  D_tau_bar <- t(BI) %*% diag(Gamma_hat[I_hat]) %*% BI
  V1 <- as.numeric(sigma_alt + t(beta_hat) %*% D_tau_bar %*% beta_hat)
  Q <- t(solve(crossprod(Theta_hat), t(Theta_hat)))
  V2 <- Omega_hat + t(Q) %*% diag(Gamma_hat) %*% Q

  K_hat <- length(I_hat_list)

  D <- c()
  ms <- c()
  for (a in 1:K_hat) {
    group_a <- unlist(I_hat_list[[a]])
    m_a <- length(group_a)
    ms[a] <- m_a

    D1 <- (2 * m_a - 1) * (D_tau_bar[a,a] ^ 2) * (m_a ^ 2) - sum(Gamma_hat[group_a] ** 2)
    D[a] <- D1 * (beta_hat[a] ^ 2) / m_a / ((m_a - 1) ^ 2)
  }

  V3 <- t(Q[I_hat, ]) %*% diag(rep(D, ms)) %*% Q[I_hat, ]
  return(diag(V1 * V2 + V3))
}
