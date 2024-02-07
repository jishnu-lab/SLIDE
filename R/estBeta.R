#' Estimate \eqn{\beta}.
#'
#' Estimate the values of the coefficients \eqn{\beta} as well as confidence intervals for the estimates
#' if requested.
#'
#' @param y vector of responses of dimension \eqn{p}
#' @param x data matrix of dimensions \eqn{n \times p}
#' @param sigma sample correlation matrix of the data, \eqn{X}
#' @param A_hat estimated matrix of dimensions \eqn{p \times K}
#' @param C_hat estimated matrix
#' @param Gamma_hat estimated matrix
#' @param I_hat estimated matrix
#' @param I_hat_list estimated matrix
#' @param alpha_level value for confidence intervals
#' @return a list including the estimates for \eqn{\beta} and the confidence intervals (if requested)
#' @export

estBeta <- function(y, x, sigma, A_hat, C_hat, Gamma_hat, I_hat, I_hat_list, alpha_level = 0.05) {
  n <- nrow(x); p <- ncol(x)
  R <- matrix(0, nrow = p, ncol = ncol(A_hat))

  BI <- t(solve(crossprod(A_hat[I_hat, ]), t(A_hat[I_hat, ])))
  #### R = Theta_hat (supplement 2.2)
  R[I_hat, ] = (sigma[I_hat, I_hat] - diag(Gamma_hat[I_hat])) %*% BI
  R[-I_hat, ] = sigma[-I_hat, I_hat] %*% BI

  #### (supplement 2.3)
  beta_est <- solve(crossprod(R), t(R) %*% crossprod(x, y) / n) ## no support

  sigma_alt <- estSigmaAlt(y = y,
                           h_hat = t(BI) %*% crossprod(x[, I_hat], y) / n,
                           beta_hat = beta_est,
                           C_hat = C_hat)
  Omega_hat <- solve(C_hat)

  ### (1 - alpha / 2) CI of beta (for each individual component)
  beta_var <- compASV(sigma_alt = sigma_alt,
                      BI = BI,
                      Theta_hat = R,
                      Gamma_hat = Gamma_hat,
                      beta_hat = beta_est,
                      Omega_hat = Omega_hat,
                      I_hat = I_hat,
                      I_hat_list = I_hat_list)
  beta_var[beta_var < 0] = 0
  alpha_level <- alpha_level / ncol(A_hat) ## Bonferroni correction
  CIs_beta <- cbind(lower = beta_est - stats::qnorm(1 - alpha_level / 2) * sqrt(beta_var) / sqrt(n),
                    upper = beta_est + stats::qnorm(1 - alpha_level / 2) * sqrt(beta_var) / sqrt(n))

  return(list(beta_hat = beta_est, conf_int = CIs_beta, beta_var = beta_var))
}
