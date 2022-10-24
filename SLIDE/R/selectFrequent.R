#' Find marginal variables using Method 2.
#'
#' Find significant marginal variables:
#' 2 - Significant marginal variable is the single most freqently appearing variable returned from \code{niter} iterations of
#'     second order knockoffs
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param niter a numeric constant; the number of times to run knockoffs
#' @param fdr a numeric constant; the target false discovery rate for knockoffs
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a vector of names corresponding to the selected columns of \code{z} that are determined to be significant
#' @export

selectFrequent <- function(z, y, niter = 1000, fdr = 0.1, parallel = TRUE) {
  ## run second order knockoffs
  results <- secondKO(z = z,
                      y = y,
                      statistic = knockoff::stat.glmnet_lambdasmax,
                      fdr = fdr,
                      niter = niter,
                      parallel = parallel)
  
  ## find most frequently selected variable
  ii <- which.max(results$tab_data)
  selected_vars <- names(results$tab_data)[ii]
  
  return (selected_vars)
}