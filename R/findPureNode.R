#' Estimate Pure Nodes.
#'
#' Estimate list of pure node indices for given \eqn{\Sigma} and \eqn{\delta}.
#' This code is an implementation of Algorithm 1 from Bing et al. (2020).
#'
#' @param abs_sigma a sample correlation matrix of dimensions \eqn{p \times p} with entries
#' in absolute value and diagonal of 0
#' @param delta \eqn{\delta}, a numerical constant
#' @param max_vals the largest absolute values of each row of \code{abs_sigma}
#' @param max_inds vector of column indices at which the values in \code{max_vals} are achieved in \code{abs_sigma}
#' @param se_est standard deviations of features (columns of \code{x})
#' @return A list including the list of estimated pure node indices and a vector
#' of the estimated pure node indices.
#' @export

findPureNode <- function(abs_sigma, delta, max_vals, max_inds, se_est) {
  G <- list() #### groups

  for (i in 1:nrow(abs_sigma)) { #### loop through rows
    row_i <- abs_sigma[i, ]
    #### Calculate indices of ith row such that the absolute values of these indices
    #### are within 2 * delta from the maximal absolute value of this row.
    #### these values correspond to columns of abs_sigma
    s_i <- findRowMaxInd(i, max_vals[i], max_inds[i], row_i, delta, se_est) ## ALG 1.4

    if (length(s_i) != 0) {
      #### For given row, check if it is a pure node by iteratively checking the other
      #### indices in s_i. Return TRUE if the given row corresponds to a pure variable.
      pure_flag <- testPure(sigma_row = row_i,
                            row_ind = i,
                            s_i = s_i,
                            max_vals = max_vals,
                            max_inds = max_inds,
                            delta = delta,
                            se_est = se_est) ## ALG 1.7

      if (pure_flag) { #### if pure
        G <- mergeUnion(G, c(s_i, i)) ## no merge()
      }
    }
  }
  return(list(pure_list = G, pure_vec = unlist(G)))
}
