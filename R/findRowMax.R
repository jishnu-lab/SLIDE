#' Find \eqn{\Sigma} Maxima.
#'
#' Calculate the maximal absolute value for each row of the given matrix.
#'
#' @param abs_sigma a matrix of dimensions \eqn{p \times p} with entries in absolute value and
#' diagonal of 0s
#' @return a list including the index of the (first) maximal absolute values of the rows of
#' \eqn{\Sigma} as well as the values
#' @export

findRowMax <- function(abs_sigma) {
  p <- nrow(abs_sigma)
  max_vals <- max_inds <- rep(0, p)
  for (i in 1:p) {
    row_i <- abs_sigma[i,]
    max_inds[i] <- which.max(row_i)
    max_vals[i] <- row_i[max_inds[i]]
  }
  return(list(max_inds = max_inds, max_vals = max_vals))
}
