#' Run aggregated second order knockoffs with splitting.
#'
#' Split data into subsets and run knockoffs on each subset separately. Aggregate results and do one final run of knockoffs to get final selected variables.
#' 
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param spec a numeric constant; the maximal difference in frequency proportions
#' @param fdr a numeric constant; the target false discovery rate for knockoffs
#' @param f_size the subset size (number of columns); must be less than or equal to \eqn{p}
#' @param niter a numeric constant; the number of times to run knockoffs
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a vector of names corresponding to the selected columns of \code{z} that are determined to be significant
#' @export

splitSecond <- function(z, y, spec = 0.1, fdr = 0.05, f_size = 100, niter = 20, parallel = T) {
  ## record sample size
  n <- nrow(z)
  ## boolean flag to monitor model size
  greater_than_n <- TRUE
  
  ## split up columns of Zs for aggregated knockoffs
  n_splits <- ceiling(ncol(z) / f_size) # number of subsets
  feature_split <- ceiling(ncol(z) / n_splits) # number of features per subset
  feature_start <- seq(1, ncol(z), feature_split) # index of first variable in subset
  feature_stop  <- pmin(feature_start + feature_split - 1, ncol(z)) # index of last variable in subset
  
  ## initial results
  screen_var <- NULL
  screen_tab <- NULL
  for (i in 1:length(feature_start)) { # iterate through subsets
    ## run second order knockoffs
    sec_ko_res <- secondKO(z = z[, c(feature_start[i]:feature_stop[i])],
                           y = y,
                           fdr = fdr,
                           niter = niter,
                           parallel = parallel)
    tab_data <- sec_ko_res$tab_data
    selected_list <- sec_ko_res$selected_list
    freq_vars <- findFrequent(spec = spec,
                              niter = niter,
                              tab_data = tab_data)
    ## find variables selected in the optimal run of knockoffs
    selected_vars <- findOptIter(n = Inf, # set Inf so that all are kept
                                 freq_vars = freq_vars,
                                 selected_list = selected_list)
    screen_var <- c(screen_var, selected_vars)
    screen_tab <- c(screen_tab, tab_data)
  }
  
  ## aggregation run
  if (is.null(screen_var)) { # if none were ever selected in any subsets, return NULL
    final_var <- NULL
  } else {
    if (length(screen_var) == 1) { # if only one variable is selected, just return it
      final_var = screen_var
    } else { # if more than one, aggregate all the results from the splits and run knockoffs again
      final_var <- secondKO(z = z[, screen_var],
                            y = y,
                            fdr = fdr,
                            niter = niter,
                            parallel = parallel)
      tab_data <- final_var$tab_data
      selected_list <- final_var$selected_list
      freq_vars <- findFrequent(spec = spec,
                                niter = niter,
                                tab_data = tab_data)
      ## find which variables were selected frequently enough
      final_var <- findOptIter(n = n,
                               freq_vars = freq_vars,
                               selected_list = selected_list)
    }
  }
  return (final_var)
}
