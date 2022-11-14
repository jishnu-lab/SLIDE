#' Find marginal variables using Method 1.
#'
#' Find significant marginal variables:
#' 1 - Significant marginal variables are those appearing in the longest list returned from \code{niter} iterations of
#'     Gaussian knockoffs.
#'     
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param niter a numeric constant; the number of times to run knockoffs
#' @param fdr a numeric constant; the target false discovery rate for knockoffs
#' @return a vector of names corresponding to the selected columns of \code{z} that are determined to be significant
#' @export

selectLongest <- function(z, y, niter = 1000, fdr = 0.1) {
  ## make data as matrices - z must have names
  z <- as.matrix(z)
  y <- as.matrix(y)
  
  ## calculate the covariance matrix and mean vector of z
  sigma <- t(z) %*% z / dim(z)[1]
  sigma <- makePosDef(sigma) # make covariance matrix positive definite
  mu <- rep(0, dim(z)[2])
  
  ## knockoffs function
  gauss_kos <- function(x) {
    knockoff::create.gaussian(X = x,
                              mu = mu, 
                              Sigma = sigma)
  }
  
  ## initialize results list
  selected_list <- list()
  knockoff_list <- list()
  
  ## run iterations in parallel
  selected_list <- foreach::foreach(i = 1:niter) %dopar% {
    ## run gaussian knockoffs
    result <- knockoff::knockoff.filter(X = z, 
                                        y = y, 
                                        knockoffs = gauss_kos, 
                                        statistic = knockoff::stat.glmnet_lambdasmax,
                                        offset = 0,
                                        fdr = fdr)
    print(result$selected)
    result$selected
  }
  
  ## find the longest list returned by iterations
  list_len <- lapply(selected_list, function(x){ length(x) })
  mm <- which.max(unlist(list_len))
  selected_vars  <- selected_list[[mm]]
  
  return (names(selected_vars))
}
