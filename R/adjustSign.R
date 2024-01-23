#' Adjust Matrix Signs.
#'
#' Perform a sign operation on matrix \eqn{\Sigma} according to the sign of \eqn{A_I}.
#' Rows corresponding to nodes not in \eqn{A_I} are set to 0.
#'
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param AI a matrix of dimensions \eqn{p \times K}
#' @return a matrix of dimensions \eqn{p \times p}
#' @export

adjustSign <- function(sigma, AI) {
  signed_sigma <- matrix(0, nrow(AI), nrow(AI))
  for (i in 1:nrow(AI)) {
    index <- which(AI[i, ] != 0)
    if (length(index) != 0) {
      signed_sigma[i, ] = sign(AI[i, index]) * sigma[i, ]
    }
  }
  return(signed_sigma)
}
