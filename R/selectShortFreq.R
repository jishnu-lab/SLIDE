#' Find significant variables according to length and frequency.
#'
#' Significant marginal variables are those appearing in the best list returned from \code{niter} iterations of
#' second order knockoffs. In this case, "best" refers to the shortest list that has the maximal overlap with the list of
#' most frequently selected variables. The frequency list is comprised of those variables appearing at least \code{spec}
#' proportion of the interactions. 
#' 
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param spec a numeric constant; the proportion of iterations that a variable must be selected in to be considered significant
#' @param fdr a numeric constant; the target false discovery rate for knockoffs
#' @param elbow a boolean flag; whether to select significant variables according to the "joint" of the elbow in a histogram or not
#' @param niter a numeric constant; the number of times to run knockoffs
#' @param f_size a numeric constant; the target size for the subset (number of columns in each subset)
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a vector of names corresponding to the selected columns of \code{z} that are determined to be significant
#' @export

selectShortFreq <- function(z, y, spec = 0.3, fdr = 0.1, elbow = FALSE, niter = 1000, f_size = 100, parallel = TRUE) {
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
    ## extract results
    tab_data <- sec_ko_res$tab_data
    selected_list <- sec_ko_res$selected_list
    ## find which variables are selected frequently enough
    if (elbow) { ## select by elbow plot
      freq_vars <- findFrequent(spec = 0.1,
                                niter = niter,
                                tab_data = tab_data)
    } else { ## select by frequency
      freq_inds <- which(tab_data >= niter * spec)
      ## extract these variables' names
      freq_vars <- colnames(tab_data)[freq_inds]
    }
    ## find variables selected in the optimal run of knockoffs
    if (length(freq_vars) > 0) {
      selected_vars <- findOptIter(freq_vars = freq_vars,
                                   selected_list = selected_list)
      screen_var <- c(screen_var, selected_vars)
      screen_tab <- c(screen_tab, tab_data)
    }
  }
  
  ## aggregation run, if needed
  if (n_splits > 1) {
    #cat("               running aggregated results . . . \n")
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
        ## find which variables are selected frequently enough
        if (elbow) { ## select by elbow plot
          ## find frequent variables
          freq_vars <- findFrequent(spec = spec,
                                    niter = niter,
                                    tab_data = tab_data)
        } else { ## select by frequency
          freq_inds <- which(tab_data >= niter * spec)
          ## extract these variables' names
          freq_vars <- colnames(tab_data)[freq_inds]
        }
        ## find which variables were selected frequently enough
        if (!is.null(freq_vars)) {
          final_var <- findOptIter(freq_vars = freq_vars,
                                   selected_list = selected_list)
        } else {
          final_var <- NULL
        }
      }
    }
  } else { ## if no aggregation - all in one run
    # cat("               no splitting . . . skipping aggregation \n")
   final_var <- screen_var 
  }
  return (final_var)
}
