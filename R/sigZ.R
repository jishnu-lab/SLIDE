#' Get Significant \eqn{Z}s. NOT USED ANYMORE
#'
#' Function to find significant \eqn{\beta}s and \eqn{Z}s from Essential Regression results.
#'
#' @importFrom magrittr '%>%'
#' @param er_res a list output by the function \link{plainER}
#' @param y a vector of responses
#' @param alpha a numeric constant for the p-value threshold
#' @return a list including the p-values for the \eqn{\beta}s, the significant cluster indices,
#' and the features found in the significant clusters
#' @export

sigZ <- function(er_res, y, alpha_level = 0.05) {
  #### get significant clusters
  p_vals <- 2 * pnorm(abs(er_res$beta) / sqrt(er_res$beta_var / length(y)), lower.tail = F)
  Z_ind <- which(p_vals <= alpha_level)
  A_hat <- er_res$A
  sig_vars <- which(rowSums(abs(A_hat[,Z_ind, drop = F])) != 0)
  return(list("p_vals" = p_vals,
              "sig_clusters" = Z_ind,
              "sig_vars" = sig_vars))
}
