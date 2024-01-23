#' Estimate \eqn{\Sigma_{\hat{TJ}}}.
#'
#' Estimate the matrix, \eqn{\Sigma_{\hat{TJ}}}.
#'
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param AI a matrix of dimensions \eqn{p \times K}
#' @param pure_vec a vector of indices of the pure nodes
#' @return an estimate for \eqn{\Sigma_{\hat{TJ}}} of dimensions \eqn{K \times |J|}
#' @export

estSigmaTJ <- function(sigma, AI, pure_vec) {
  #### adjust sign of entries in sigma according to AI
  signed_sigma <- adjustSign(sigma = sigma, AI = AI)
  #### subset columns of signed_sigma to just nonpure variables
  sigma_J <- matrix(signed_sigma[, -pure_vec], nrow = nrow(sigma))
  sigma_TJ <- matrix(0, ncol(AI), nrow(AI) - length(pure_vec))
  for (i in 1:ncol(AI)) {
    group_i <- which(AI[, i] != 0) #### get pure nodes in cluster i
    #### subset sigma_J to just rows of pure variables in cluster i
    sigma_iJ <- as.matrix(sigma_J[group_i, ])
    #### set row i of sigma_TJ to the mean correlations of the
    #### nonpure variables across the pure variables in cluster i
    sigma_TJ[i, ] <- apply(sigma_iJ, 2, mean)
  }
  return(sigma_TJ)
}
