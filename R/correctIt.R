#' Correct the response for the marginal when finding interactions.
#'
#' The new response vector is the residuals after regressing the true response values on the values of
#' one marginal \eqn{Z}. This is used for determining which interaction terms that include this marginal
#' variable are truly significant.
#' 
#' @param y a vector or data frame with one column of numeric values; the true response
#' @param z a vector or data frame with one column of numeric values; the values of a marginal variable
#' @return a vector of numeric values
#' @export

correctIt <- function(z, y) {
  y_hat <- stats::residuals(stats::lm(y ~ z))
  return(y_hat)
}