#' Find \eqn{\Sigma} Maximum Indices.
#'
#' Calculate indices of each row such that the absolute values of these indices
#' are within \eqn{2\delta} of the maximal absolute value \code{max_val} of this row.
#' This is an implementation of step 4 of Algorithm 1 in Bing et al. (2020).
#'
#' @param i the row index
#' @param max_val the maximal absolute value of row \eqn{i} of the covariance/correlation matrix
#' @param max_ind the first index in row \eqn{i} at which \code{max_val} is achieved
#' @param row_i a row of \code{abs_sigma}
#' @param delta \eqn{\delta}, a numeric constant
#' @param se_est vector of standard deviations of features (columns of \code{x})
#' @return a vector of indices
#' @export

findRowMaxInd <- function(i, max_val, max_ind, row_i, delta, se_est) {
  ## lbd <- delta * sd(feat i) * sd(feat max_val) + delta * sd(feat i) * sd(each feat)
  ## lbd is a vector of values for each column abs_sigma (each feature)
  lbd <- delta * se_est[i] * se_est[max_ind] + delta * se_est[i] * se_est

  ## which features satisfy max_val ≤ lbd + value in row_i
  ## ALG 1.4 - find set of indices within 2delta of max
  indices <- which(max_val <= lbd + row_i)
  return(indices)
}
