#' Find marginal variables.
#'
#' Find significant marginal variables. There are currently six supported methods:
#' top_betas - Significant marginal variables are selected as the variables corresponding to the absolute largest coefficients from Essential Regression.
#'     Selected variables are in the top \eqn{\frac{\text{spec}}{2}} proportion of the positive or negative largest coefficients.
#' marginals - Significant marginal variables are selected as the provided vector of indices. This requires passing in a vector of indices as \code{marginals}.
#' 1 - Significant marginal variables are those appearing in the longest list returned from \code{niter} iterations of
#'     Gaussian knockoffs
#' 2 - Significant marginal variable is the single most freqently appearing variable returned from \code{niter} iterations of
#'     second order knockoffs
#' 3 - Significant marginal variables are those appearing in the best list returned from \code{niter} iterations of
#'     second order knockoffs. In this case, "best" refers to the shortest list that has the maximal overlap with the list of
#'     most frequently selected variables. The frequency list is comprised of those variables appearing at least \code{spec}
#'     proportion of the interactions. If the number of selected significant marginal variables is longer than the number of samples
#'     (number of rows of \code{z}), then \code{spec} is increased by 0.1, and the selected process is repeated.
#' 4 - Significant marginal variables are selected in the same way as in Method 3 but without reselection in the case of too long a 
#'     list.
#'     
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param method an integer; the selected method to use
#' @param betas a vector of numeric values; the coefficients corresponding to the variables in \code{z} (found using Essential Regression)
#' @param top_prop a numeric constant between 0 and 1; the proportion of \code{betas} to consider significant
#' @param marginals a vector of numeric values; the indices of the variables to consider significant
#' @param elbow a boolean flag; whether to select using \code{spec} or by identifying the "joint" of the histogram elbow
#' @param spec a numeric constant between 0 and 1; the proportion of iterations that a variable must be selected in to be considered significant
#' @param fdr a numeric constant between 0 and 1; the target false discovery rate for knockoffs
#' @param niter an integer; the number of times to run knockoffs
#' @param f_size an integer; the target size for the subset (number of columns in each subset)
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a vector of names corresponding to the selected columns of \code{z} that are determined to be significant
#' @export

marginalSLIDE <- function(z, y, method = 4, betas = NULL, top_prop = NULL, marginals = NULL, elbow = FALSE,
                          spec = 0.5, fdr = 0.1, niter = 1000, f_size = 100, parallel = TRUE) {
  #### check that data is satisfactory
  if (is.null(z) || is.null(y)) {
    stop("z and y must both be provided \n")
  } else if (fdr < 0 || fdr > 1) {
    stop("fdr must be between 0 and 1 \n")
  }
  
  ## get model size maximum
  n <- nrow(z)
  ## boolean flag for model size monitoring
  greater_than_n <- TRUE
  while (greater_than_n) {
    if (!is.null(top_prop)) {
      cat("      selecting marginal variables using method 5 . . . \n")
      marginal_vars <- getTopBetas(betas = betas, top_prop = top_prop)
    } else if (!is.null(marginals)) {
      cat("      selecting marginal variables using method 6 . . . \n")
      cat("      randomly picking marginals . . . \n")
      marginal_vars <- marginals
    } else if (method == 1) { ## select marginal Zs by longest result from niter replicates of gaussian KOs
      cat("      selecting marginal variables using method 1 . . . \n")
      marginal_vars <- selectLongest(z = z, 
                                     y = y, 
                                     niter = niter, 
                                     fdr = fdr)
    } else if (method == 2) { ## select marginal Z as the most frequently appearing Z from niter replicates of second order KOs
      cat("      selecting marginal variables using method 2 . . . \n")
      marginal_vars <- selectFrequent(z = z,
                                      y = y,
                                      niter = niter,
                                      fdr = fdr,
                                      parallel = parallel)
    } else if (method == 3) { ## select marginal Z by frequency/shortest list with subsetting
      cat("      selecting marginal variables using method 3 . . . \n")
      marginal_vars <- selectShortFreq(z = z,
                                       y = y,
                                       niter = niter,
                                       f_size = f_size,
                                       spec = spec,
                                       fdr = fdr,
                                       parallel = parallel)
    } else if (method == 4) { ## select marginal Z by frequency/shortest list with subsetting
     # cat("      selecting marginal variables using method 4 . . . \n")
      marginal_vars <- selectShortFreq(z = z,
                                       y = y,
                                       niter = niter,
                                       elbow = elbow,
                                       f_size = f_size,
                                       spec = spec,
                                       fdr = fdr,
                                       parallel = parallel)
    } else {
      stop("Invalid method provided. Please provide a marginal variable selection method. \n")
    }
    
    ## exit loop if model is small enough
    if (length(marginal_vars) < n) {
      greater_than_n <- FALSE
    } else {
      cat("         too many marginal terms . . . \n")
      if (elbow) { ## selecting with elbow requires spec to be decreased for more stringent model
        cat("             decreasing spec . . . \n")
        spec <- spec - 0.05
      } else { ## selecting with other methods requires spec to be increased for more stringent model
        cat("             increasing spec . . . \n")
        spec <- spec + 0.1
      }
    }
  }
  
  cat(paste0("         final marginal spec: ", spec, "\n"))
  
  return (marginal_vars)
}


