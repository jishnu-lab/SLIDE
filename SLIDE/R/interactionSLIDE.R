#' Find significant interaction terms.
#'
#' Iterate through marginal variables and construct and test their interaction terms.
#' Concatenate all marginal variable results.
#'
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param m an integer; the target model size
#' @param method an integer; the selected method to use
#' @param elbow a boolean flag; whether to select using \code{spec} or by identifying the "joint" of the histogram elbow
#' @param spec a numeric constant between 0 and 1; the proportion of iterations that a variable must be selected in to be considered significant
#' @param fdr a numeric constant between 0 and 1; the target false discovery rate for knockoffs
#' @param niter an integer; the number of times to run knockoffs
#' @param f_size an integer; the target size for the subset (number of columns in each subset)
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a list containing the names of the significant interaction terms, the values of those terms, and the values of \eqn{\upsilon}
#' (NULL if not method 4)
#' @export

interactionSLIDE <- function(z, y, m, marginals, method = 4, elbow = FALSE, spec = 0.3, 
                             fdr = 0.1, niter = 1000, f_size = 100, parallel = TRUE) {
  ## set flag to monitor model size
  greater_than_m <- TRUE
  
  ## loop until model is small enough
  while (greater_than_m) {
    cat("starting interactions........ \n")
    ## initialize vectors of results
    sig_interacts <- NULL # a vector; significant interaction variables
    interact_terms <- NULL # a matrix; values of sig_interacts
    upsilons <- matrix(nrow = nrow(z), ncol = 0) # a matrix; values of upsilons
    upsilon_mods <- NULL
    used_margs <- NULL # a vector; monitor which marginals have already been used
    
    ## loop through marginals and do interactions separately
    for (i in 1:length(marginals)) {
      ## get one marginal variable
      marg_var <- marginals[i]
      ## add it to the list of already used maorginal variables
      used_margs <- c(used_margs, marg_var)
      
      ## construct interaction terms with union logic
      cat(paste0("          interaction terms with ", marg_var), "\n")
      int_list <- interUnion(marginal_vars = marg_var,
                             z = z)
      
      ## determine if the reverse interaction term has already been tested
      ## for example, if we already tested Var1.Var2, we do not want to consider Var2.Var1
      already_made <- paste0(marg_var, ".", used_margs)
      interactions <- int_list$interaction # extract interaction term values
      which_already_made <- which((colnames(interactions) %in% already_made) == TRUE)
      if (length(which_already_made) > 0) {
        ## remove columns if the interaction term was tested previously
        interactions <- interactions[, -which_already_made]
      }
      
      ## correct y for the marginal variable
      corrected_y <- correctIt(z = as.matrix(z[, marg_var]),
                               y = y)
      
      if (method == 1) { # Method 1: longest iteration
        ## same routine as with marginals, use same code
        interaction_var <- selectLongest(z = interactions,
                                         y = corrected_y,
                                         fdr = fdr,
                                         niter = niter)
      } else { ## Methods 2, 3, 4: shortest, best iteration
        ## same routine as with marginals, use same code
        interaction_var <- selectShortFreq(z = interactions,
                                           y = corrected_y,
                                           fdr = fdr,
                                           elbow = elbow,
                                           niter = niter,
                                           spec = spec,
                                           parallel = parallel)
      }
      
      ## set flag, model size is small enough
      if (length(interaction_var) <= m) {
        greater_than_m <- FALSE
      } else {
        cat("         too many interaction terms . . . \n")
        if (elbow) { # elbow needs to decrease spec to be more stringent (smaller model size)
          cat("             decreasing spec . . . \n")
          spec <- spec - 0.05
        } else { # other methods need to increase spec to be more stringent (smaller model)
          cat("             increasing spec . . . \n")
          spec <- spec + 0.1
        }
        greater_than_m <- TRUE
        break
      }
      
      ## Method 4 requires additional work to construct the temporary variables via linear regression
      if (method == 4) {
        ## create upsilon (fitted y values when regression on marginal and its interaction terms)
        upsilon <- getUpsilon(z = z,
                              y = y,
                              interaction_var = interaction_var,
                              interactions = interactions,
                              marg_var = marg_var)
        ## append results
        upsilons <- cbind(upsilons, upsilon$upsilon)
        upsilon_mods[[paste0("ups", marg_var)]] <- upsilon$upsilon_mod
      }
      
      ## if interaction terms were selected, then concatenate them with past results
      if (!is.null(interaction_var)) {
        ## add interaction term names
        sig_interacts <- c(sig_interacts, interaction_var)
        ## add interaction term values
        int_var_df <- int_list$interaction[, interaction_var]
        ## rename interaction term value column names
        old_names <- colnames(interact_terms)
        interact_terms <- cbind(interact_terms, int_var_df)
        colnames(interact_terms) <- c(old_names, interaction_var)
      }
    }
  }
  
  # cat(paste0("         final interaction spec: ", spec, "\n"))
  # 
  # debugFileName <- paste0("interactionSLIDE_debug", rnorm(n = 1), ".rds")
  # 
  # cat(paste0("The SLIDE funtion debugg file name is ", debugFileName))
  #  
  # 
  #  saveRDS(list(
  #     z=z,
  #     y=y, 
  #     m=m, 
  #     marginals=marginals,
  #     method = method, 
  #     elbow = FALSE, 
  #     spec = spec, 
  #     fdr = fdr, 
  #     niter = niter,
  #     f_size = f_size, 
  #     parallel = TRUE,
  #     upsilon = upsilons,
  #     upsilon_mods = upsilon_mods), 
  #     file = paste0("/ix/djishnu/Javad/Poholek_scale_v_noscale/w_scale_220922/DebugData/SLIDE/interactionSlide/", debugFileName),
  #     )
  # 
  
  
  
  return (list(interaction_vars = sig_interacts,
               interaction_vals = interact_terms,
               upsilon = upsilons,
               upsilon_mods = upsilon_mods))
}