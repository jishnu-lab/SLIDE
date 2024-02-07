#' estSigmaAlt.
#'
#' I am not sure what this does.
#'
#' @param y a response vector of dimension \eqn{p}
#' @param h_hat estimated matrix:
#' \eqn{\frac{1}{n}(\hat{A}_{\hat{I}}^\top \hat{A_{\hat{I}}^\top)^{-1}\hat{A}_{\hat{I}}^\top \mathbf{X}_{\hat{I}}^\top \mathbf{Y}})}
#' @param beta_hat estimated values of \eqn{\beta}
#' @param C_hat estimated matrix
#' @return ???
#' @export

estSigmaAlt <- function(y, h_hat, beta_hat, C_hat) {
  n <- length(y)
  sigma_alt <- crossprod(y) / n - 2 * t(beta_hat) %*% h_hat + t(beta_hat) %*% C_hat %*% beta_hat
  return (ifelse(sigma_alt < 0, 0, sigma_alt))
}
