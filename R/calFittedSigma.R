#' Calculate Fitted Sigma.
#'
#' Calculate the fitted value of \eqn{A_I \cdot C \cdot A_I^\top} for given
#' \eqn{\Sigma} and \eqn{\delta}.
#'
#' @param sigma a correlation matrix of dimensions \eqn{p \times p}
#' @param delta a threshold parameter
#' @param max_vals the calculated maximal values of \eqn{\Sigma} by row
#' @param max_inds vector of column indices at which the values in \code{max_vals} are achieved in \code{| sigma |}
#' @param se_est the estimated standard errors
#' @return a list containing a vector of the indices of the estimated pure variables and
#' the fitted value of \eqn{A_I \cdot C \cdot A_I^\top}. returns -1 if only one pure node is identified
#' @export

calFittedSigma <- function(sigma, delta, max_vals, max_inds, se_est) {
  #### find pure nodes
  pure_nodes <- findPureNode(abs_sigma = abs(sigma),
                             delta = delta,
                             max_vals = max_vals,
                             max_inds = max_inds,
                             se_est = se_est)
  pure_list <- pure_nodes$pure_list

  if (singleton(pure_list = pure_list)) {
    return(list(pure_vec = NULL, fit_sigma = -1))
  }

  signed_pure_list <- findSignPureNode(pure_list = pure_list, sigma = sigma)
  AI <- recoverAI(pure_list = signed_pure_list, p = length(se_est))
  C <- estC(sigma = sigma, AI = AI)

  if (length(pure_list) == 1) {
    fit_sigma <- -1
  } else {
    sub_AI <- AI[pure_nodes$pure_vec, ]
    fit_sigma <- sub_AI %*% C %*% t(sub_AI)
  }
  return(list(pure_vec = pure_nodes$pure_vec, fit_sigma = fit_sigma))
}
