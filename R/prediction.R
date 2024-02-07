#' Essential Regression Prediction.
#'
#' Perform prediction on the data using the ER predictor.
#'
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param A_hat a matrix of dimensions \eqn{p \times K}
#' @param Gamma_hat a matrix of dimensions
#' @param I_hat a matrix of dimensions
#' @return a list including \eqn{\hat{\theta}}, the predicted values, and \eqn{\Theta}
#' @export

prediction <- function(y, x, sigma, A_hat, Gamma_hat, I_hat) {
  K <- ncol(A_hat) #### number of clusters
  p <- ncol(x); n <- nrow(x)
  #### R = K x p
  R <- matrix(0, nrow = K, ncol = p)

  #### solve for X: t(A_hat[I_hat, ]) %*% A_hat[I_hat, ] %*% X = t(A_hat[I_hat, ])
  #### A_hat = p x K, A_hat[I_hat, ] = |I_hat| x K
  #### BI = (K x p) %*% (p x K) %*% X = (K x p) -> (K x p)
  BI <- solve(t(A_hat[I_hat, ]) %*% A_hat[I_hat, ], t(A_hat[I_hat, ]))
  #### R = Theta_hat (supplement 2.3)
  #### sigma = p x p
  #### sigma[I_hat, I_hat] = |I_hat| x |I_hat|
  #### Gamma_hat = p x ?
  #### R = K x p
  R[, I_hat] = BI %*% (sigma[I_hat, I_hat] - diag(Gamma_hat[I_hat]))
  R[, -I_hat] = BI %*% sigma[I_hat, -I_hat]
  #### Q = n x K
  #### Q = (theta_hat^T %*% X^T)
  Q <- x %*% t(R)
  #### theta_hat = theta_hat_ER (supplement 2.2.1)
  theta_hat <- try(t(R) %*% solve(crossprod(Q), crossprod(Q, y)), silent = T)
  if (class(theta_hat)[1] == "try-error") {
    theta_hat <- t(R) %*% MASS::ginv(crossprod(Q)) %*% crossprod(Q, y)
  }

  #### pred_val is Y*_hat_ER (supplement 2.6)
  pred_val <- x %*% theta_hat
  return(list(er_predictor = theta_hat, pred_vals = pred_val, theta_hat = t(R)))
}
