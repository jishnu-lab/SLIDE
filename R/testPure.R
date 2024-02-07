#' Test Node Purity.
#'
#' For given row, check if it is a pure node by iteratively checking the nodes
#' in \code{Si}. Return TRUE if the given row corresponds to a pure variable.
#' This is an implementation of steps 6-9 of Algorithm 1 in Bing et al. (2020).
#'
#' @param sigma_row a vector of dimension \eqn{p}
#' @param row_ind the index of the row
#' @param s_i a vector of indices that are within \eqn{2\delta} of the maximum value of \code{sigma_row}
#' @param max_vals a vector of the largest absolute values of each of the rows in \code{sigma}
#' @param max_inds a vector of the first index at which each value in Ms is achieved
#' @param delta \eqn{\delta}, a numerical constant
#' @param se_est a vector of estimates of the standard deviations of the rows of the data matrix, \eqn{x}
#' @return TRUE or FALSE
#' @export

testPure <- function(sigma_row, row_ind, s_i, max_vals, max_inds, delta, se_est) {
  for (i in 1:length(s_i)) {
    #### go row by row through list of indices where abs val is within 2*delta of max
    j <- s_i[i] #### j is some index of a row in sigma that has large enough abs cov/corr
    #### some sort of cutoff value
    delta_j <- (se_est[row_ind] + se_est[max_inds[j]]) * se_est[j] * delta
    #### check if abs diff between jth entry in row i and max value in row j is greater than cutoff
    #### sigma_row[j] = sigma_{ij} (ith by jth entry)
    #### maxes[j] = max(sigma_{j.}) (max of jth row)
    if (abs(sigma_row[j] - max_vals[j]) > delta_j) {
      return(FALSE)
    }
  }
  return(TRUE)
}
