#' FDR Thresholding of \eqn{\hat{Sigma}}.
#'
#' Threshold \eqn{\hat{\Sigma}} according to false discovery rate threshold specification \code{thresh}.
#' For each entry in \eqn{\hat{\Sigma}}, calculate the associated \eqn{t}-test statistic with degrees of freedom \eqn{n-2}:
#' \deqn{t = \frac{r\sqrt{n-2}}{\sqrt{1-r^2}}}
#' Then find the two-tailed p-value for this according to the \eqn{t}-distribution, correct for
#' multiple testing, and threshold \eqn{\hat{\Sigma}} by setting any entry with an associated p-value
#' less than or equal to \code{thresh} to 0.
#'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param thresh a numeric constant
#' @return a matrix of dimensions \eqn{p \times p}
#' @export

threshSigma <- function(x, sigma, thresh) {
  #### get dimensions of data matrix
  n <- nrow(x)
  p <- ncol(x)

  #### find p-values of entries in sigma
  sigma_pvals <- matrix(0, nrow = nrow(sigma) , ncol = ncol(sigma))

  #### find p-values for the entries in sigma
  res <- apply(sigma, c(1, 2), function(x){corrToP(x, n = n)$p})

  #### perform multiple testing correction: Benjamini-Hochberg
  res_adjust <- as.data.frame(matrix(stats::p.adjust(res, method = 'BH'), nrow = nrow(res)))

  #### adjust sigma according to the p-values
  res_filtered <- apply(res_adjust, c(1, 2), function(x){ifelse (x <= thresh, 1, 0)})
  sigma <- sigma * res_filtered
  return (list("thresh_sigma" = sigma,
               "kept_entries" = res_filtered))
}


