#' Estimate \eqn{A_I}.
#'
#' Use the given \eqn{\delta} to calculate the fitted \eqn{A_I}, pure variable
#' indices in list form and vector form. Also return estimated \eqn{Y} and \eqn{C} for
#' the following Dantzig estimation.
#'
#' @param sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param delta \eqn{\delta}, a numerical constant
#' @param se_est standard deviations of the columns of the data matrix, \code{x}
#' @return a list including: \eqn{A_I}, a matrix of dimensions \eqn{p \times K},
#' a vector of the indices of estimated pure variables, and a list of the indices of
#' the estimated pure variables
#' @export


estAI <- function(sigma, delta, se_est) {
  #### get absolute value of covariance matrix
  abs_sigma <- abs(sigma)

  #### set entries on main diagonal to 0
  diag(abs_sigma) <- 0

  #### calculate the maximal absolute value for each row of sigma
  result_max <- findRowMax(abs_sigma = abs_sigma)
  max_vals <- result_max$max_vals #### maximal abs values
  max_inds <- result_max$max_inds #### first index where max abs values are achieved

  #### estimate list of pure node indices for given sigma and delta
  result_pure <- findPureNode(abs_sigma = abs_sigma,
                              delta = delta,
                              max_vals = max_vals,
                              max_inds = max_inds,
                              se_est = se_est)
  pure_list <- result_pure$pure_list
  pure_vec <- as.vector(unlist(result_pure$pure_vec))

  #### Estimate the sign subpartition of pure node sets. If there is an element
  #### of a list is empty, then a empty list will be put in that position
  signed_pure_list <- findSignPureNode(pure_list = pure_list, sigma = sigma)

  #### Recover the estimated submatrix A_I given the pure node group.
  AI <- recoverAI(pure_list = signed_pure_list, p = nrow(abs_sigma))

  return(list(AI = AI, pure_vec = pure_vec, pure_list = signed_pure_list))
}
