#' Run aggregated second order knockoffs.
#'
#' Iteratively run knockoffs and report variable selection frequencies.
#' 
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param niter a numeric constant; the number of times to run knockoffs
#' @param fdr a numeric constant; the target false discovery rate for knockoffs
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a list including the names of the selected columns of \code{z} that are determined to be significant and a table of their selection frequencies
#' @export

secondKO <- function(z, y, niter = 100, fdr = 0.2, parallel = TRUE) {
  ## Z and Y sanity check
  if (is.null(z) || is.null(y)) {
    stop("Z and Y must not be empty.")
  }
  
  if (parallel == T) { ## parallel computing
    selected_list <- foreach::foreach(i = 1:niter) %dopar% {
      result <- knockoff::knockoff.filter(X = z,
                                          y = as.matrix(y), 
                                          knockoffs = knockoff::create.second_order, 
                                          statistic = knockoff::stat.glmnet_lambdasmax,
                                          offset = 0,
                                          fdr = fdr)
      names(result$selected)
    }
  } else { ## sequential computing
    selected_list <- list()
    selected_list <- foreach::foreach(i = 1:niter) %do% {
      result <- knockoff::knockoff.filter(X = z, 
                                          y = y, 
                                          knockoffs = knockoff::create.second_order, 
                                          statistic = knockoff::stat.glmnet_lambdasmax,
                                          offset = 0,
                                          fdr = fdr)
    }
  }
  
  ## create frequency table
  tab_data <- table(unlist(selected_list)) 
  return (list(selected_list = selected_list,
               tab_data = as.matrix(t(tab_data))))
}
