#' Run SLIDE.
#'
#' Run SLIDE.
#'
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param method an integer between 1 and 6; the selected method to use
#' @param do_interacts a boolean flag; whether to select interaction terms or not
#' @param betas a vector of numeric values; the coefficients corresponding to the variables in \code{z} (found using Essential Regression)
#' @param top_prop a numeric constant between 0 and 1; the proportion of \code{betas} to consider significant
#' @param marginals a vector of integers; the indices of the variables to consider significant
#' @param spec a numeric constant; the proportion of iterations that a variable must be selected in to be considered significant
#' @param fdr a numeric constant between 0 and 1; the target false discovery rate for knockoffs
#' @param niter an integer; the number of times to run knockoffs
#' @param elbow a boolean flag; whether to select significant variables according to the "joint" of the elbow in a histogram or not
#' @param f_size an integer; the target size for the subset (number of columns in each subset)
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @param ncore an integer; the number of cores to use for parallel computing
#' @export

SLIDE <- function(z, y, method = 4, do_interacts = TRUE, betas = NULL, top_prop = NULL, marginals = NULL, 
                  spec = 0.3, fdr = 0.1, niter = 1000, elbow = FALSE, f_size = 100, parallel = TRUE, ncore = 10) {
  #### HOUSEKEEPING ############################################################
  ## record number of samples
  require(dplyr)
  n <- nrow(z)
  y <- as.matrix(y)
  
  ## boolean flag for monitoring model size
  greater_than_n <- TRUE
  
 
  
  
  
  ## do parallel computing if specified
  if (parallel == T) {
    cl <- parallel::makePSOCKcluster(ncore)
    doParallel::registerDoParallel(cl)
  }
  
  #### MARGINAL SELECTION ######################################################
  ## select marginal variables
  colnames(z) <- paste0("z",1:ncol(z))
  marginal_vars <- marginalSLIDE(z = z,
                                 y = y,
                                 method = method,
                                 niter = niter,
                                 f_size = f_size,
                                 parallel = parallel,
                                 fdr = fdr,
                                 elbow = elbow,
                                 spec = spec,
                                 marginals = marginals,
                                 betas = betas,
                                 top_prop = top_prop)
  
  
  
  
  
  
  ## if no marginal variables are selected, skip making interaction terms
  if (is.null(marginal_vars)) {
    cat("    no interaction terms . . . no marginals \n")
    if (parallel == T) { ## stop parallel computing
      parallel::stopCluster(cl)
    }
    return (list("marginal_vars" = marginal_vars,
                 "interaction_vars" = NULL,
                 "interactions" = NULL))
  }
  
  ## if too many marginal variables are selected, skip making interaction terms 
  m <- n - length(marginal_vars)
  if (m <= 0) {
    cat("    no interaction terms . . . too many marginals \n")
    if (parallel == T) { ## stop parallel computing
      parallel::stopCluster(cl)
    }
    return(list("marginal_vars" = marginal_vars,
                "interaction_vars" = NULL,
                "interactions" = NULL))
  }
  
  #### INTERACTION SELECTION ###################################################
  if (!do_interacts) { ## don't do interactions
    if (parallel == T) { ## stop parallel computing
      parallel::stopCluster(cl)
    }
    return (list("marginal_vars" = marginal_vars,
                 "interaction_vars" = NULL,
                 "interactions" = NULL))
  }
  
  cat("      starting interaction selection . . . \n")
  ## do interaction term selection
  print("Before doing interaction SLIDE")
  print(marginal_vars)
  marginal_vars<- as.numeric(sub("z","",marginal_vars))
  interactions <- interactionSLIDE(z = z,
                                   y = y,
                                   m = m,
                                   method = method,
                                   elbow = elbow,
                                   marginals = marginal_vars,
                                   spec = spec,
                                   niter = niter,
                                   f_size = f_size,
                                   parallel = parallel, 
                                   fdr = fdr)
  
  print("printig the yhat of each maginals:")
  print(interactions$upsilons)
  
  #### METHOD 4 ADDITIONAL WORK ################################################
  if (method == 4 && !is.null(interactions$upsilon)) {
    cat("      running knockoffs on marginal/interaction submodels . . . \n")
    final_upsilon <- testUpsilon(upsilon = interactions$upsilon,
                                 y = y,
                                 niter = niter,
                                 elbow = elbow,
                                 spec = spec,
                                 fdr = fdr, 
                                 f_size = f_size, 
                                 parallel = parallel)
    print("upsilon colnames:")
    colnames(final_upsilon)
    
    final_mods <- interactions$upsilon_mods[colnames(final_upsilon)]
  } else {
    final_upsilon <- NULL
    final_mods <- NULL
  }
  
  if (parallel == T) { ## stop parallel computing
    parallel::stopCluster(cl)
  }
  
  
  
  
  #### RETURN ##################################################################
  return (list("marginal_vars" = marginal_vars, # significant marginal variables
               "interaction_vars" = interactions$interaction_vars, # significant interaction variables
               "interaction_vals" = interactions$interaction_vals, # significant interaction variable values
               "upsilon" = final_upsilon,
               "upsilon_mods" = final_mods)) # method 4 only
  
}
