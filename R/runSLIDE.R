#' Run SLIDE for cross-validation.
#'
#' Creates a regression model from the results of SLIDE and
#' predicts response values from the validation data. 
#' 
#' @param train_y a vector of numeric values; the response (training set)
#' @param valid_y a vector of numeric values; the response (validation set)
#' @param train_z a matrix or data frame of numeric values; the latent factor data (training set)
#' @param valid_z a matrix or data frame of numeric values; the latent factor data (validation set)
#' @param method an integer (1, 2, 3, 4); the selection method for determining marginal/interaction significance
#' @param spec a numeric constant; the proportion of iterations that a variable must be selected in to be considered significant
#' @param niter an integer; the number of times to run knockoffs
#' @param fdr a numeric constant between 0 and 1; the target false discovery rate for knockoffs
#' @param f_size an integer; the target size for the subset (number of columns in each subset)
#' @param betas a vector of numeric values; the coefficients corresponding to the variables in \code{z} (found using Essential Regression)
#' @param top_prop a numeric constant between 0 and 1; the proportion of \code{betas} to consider significant
#' @param marginals a vector of integers; the indices of the variables to consider significant
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @param ncore an integer; the number of cores to use for parallel computing
#' @param y_factor a boolean flag; whether the response is categorical or not
#' @return a list containing the training model and the predicted values of the response
#' using \code{valid_x}
#' @export

runSLIDE <- function(train_y, valid_y, train_z, valid_z, method, spec, niter, fdr, 
                     f_size, betas, top_prop, marginals, parallel, ncore, y_factor) {
  ## run SLIDE
  res <- SLIDE(z = train_z,
               y = train_y,
               method = method,
               niter = niter,
               spec = spec,
               marginals = marginals,
               parallel = parallel,
               f_size = f_size,
               betas = betas,
               top_prop = top_prop,
               ncore = ncore,
               elbow = FALSE,
               do_interacts = TRUE,
               fdr = fdr)
  ## if no sig margs, select 5 random ones (like in lasso)
  if (length(res$marginal_vars) == 0) {
    cat("SLIDE selects no features - Randomly selecting 5 features. . . \n")
    sig_margs <- sample(1:ncol(train_z), 5)
    res <- SLIDE(z = train_z,
                 y = train_y,
                 method = method,
                 niter = niter,
                 spec = spec,
                 marginals = sig_margs,
                 parallel = parallel,
                 f_size = f_size,
                 betas = betas,
                 top_prop = top_prop,
                 ncore = ncore,
                 do_interacts = TRUE,
                 fdr = fdr)
  }
  ## get training Zs (margs + interacts)
  if (!is.null(res$upsilon)) {  ## interactions found
    zs_train <- res$upsilon
    ## get validation interactions
    zs_ints_valid <- matrix(nrow = nrow(valid_z), ncol = 0)
    
    
    
    
    saveRDS(list(
      train_z = train_z,
      train_y = train_y,
      valid_z = valid_z,
      res = res
    ), file = paste0("/ix/djishnu/Javad/Poholek_scale_v_noscale/w_scale_220922/DebugData/", debugFileName))
    
    ## loop through marginals and create upsilons by regressing on marg + its interactions
    ## and extracting the fitted values from the regression model
    
    print(zs_train)
    print("/n")
    print(dim(zs_train))
    zs_train <- as.data.frame(zs_train)
    
    
    for (t in 1:ncol(zs_train)) {
      ## get name of marginal and upsilon
      marg_var <- colnames(zs_train)[t]
      only_ind <- gsub("upsZ", "", marg_var)
      marg_var <- gsub("ups", "", marg_var)
      ## make interaction terms
      print("possible error:")
      print(colnames(valid_z))
      
      inter_union <- interUnion(marginal_vars = marg_var, z = valid_z)
      ## get interaction variable names that were chosen as significant
      selected_interactions <- res$interaction_vars
      ## remove interacting variables (get just marginals)
      only_marg <- gsub("\\..*","", selected_interactions)
      ## get indices of interaction terms with the target marginal
      sel <- which(gsub("Z","",only_marg) == gsub("Z","",marg_var))
      ## extract the values of these chosen interaction terms
      selected_interactions <- inter_union$interaction[, selected_interactions[sel]]
      ## make a data frame from the target marginal and interactions
      print(marg_var)
      # if(!grepl("Z",marg_var)){marg_var <- paste0("Z",marg_var)}
      # print(colnames(valid_z))
      
      
      
      ## Sometime variables are coming with Z and sometimes not so this section do double checking!
      
      if(!grepl("Z",colnames(res$upsilon)[1])){
      colnames(valid_z) <- gsub("Z","",colnames(valid_z))
      marg_var          <- gsub("Z","",marg_var)}else{
        
        if(!grepl("Z",colnames(valid_z)[1]) & !grepl("Z",marg_var) ){colnames(valid_z) <- paste0("Z",colnames(valid_z)); 
        print("new columns are ") 
        marg_var <- paste0("Z",marg_var)}
        
        if(!grepl("Z",colnames(valid_z)[1]) & grepl("Z",marg_var) ){colnames(valid_z) <- paste0("Z",colnames(valid_z)); 
        print("new columns are ")
        colnames(valid_z)
        }
        
        if(grepl("Z",colnames(valid_z)[1]) & !grepl("Z",marg_var) ){ marg_var <- paste0("Z",marg_var); 
        sprintf("new marginal is %s",marg_var) 
        }
        
        
      }
      
      valid_df <- data.frame(valid_z[, marg_var], selected_interactions)
      colnames(valid_df) <- c(marg_var, res$interaction_vars[sel])
      
      ## predict using trained model
      upsilon_name <- paste0("ups", marg_var)
      print(res$upsilon_mods[[upsilon_name]])
      fit_vals <- stats::predict(object = res$upsilon_mods[[upsilon_name]], newdata = valid_df)
      ## extract fitted values and rename
      valid_upsilon <- fit_vals %>%
        as.data.frame()
      colnames(valid_upsilon) <- upsilon_name
      zs_ints_valid <- cbind(zs_ints_valid, valid_upsilon)
    }
    
    ## extract the significant upsilons
    selected_upsilon <- colnames(res$upsilon)
    zs_valid <- zs_ints_valid[, selected_upsilon] %>%
      as.data.frame()
    colnames(zs_valid) <- selected_upsilon
  } else { ## no interactions
    cat("NO INTERACTIONS...USING MARGINALS \n")
    ## just use the marginal variables and make into a data frame
    zs_train <- cbind(train_z[, res$marginal_vars]) %>%
      as.data.frame()
    colnames(zs_train) <- c(res$marginal_vars)
    zs_valid <- cbind(valid_z[, res$marginal_vars]) %>%
      as.data.frame()
    colnames(zs_valid) <- c(res$marginal_vars)
  }
  ## make model
  if (y_factor) {
    ## logistic regression if y is categorical
    slide_mod <- stats::glm(train_y ~ ., data = zs_train, family = "binomial")
    ## predict values for validation set
    valid_prob <- predict(slide_mod, zs_valid, type = "response")
    ## get contrasts
    contrast <- stats::contrasts(train_y)
    ## get labels from class probabilities - validation set
    valid_pred <- ifelse(valid_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
    
    ## predict values for training set
    train_prob <- slide_mod$fitted.values
    ## get labels from class probabilities - training set
    train_pred <- ifelse(train_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
  } else {
    ## simple linear regression if y is continuous
    slide_mod <- stats::lm(train_y ~ ., data = zs_train)
    ## predict values for validation set
    valid_pred <- predict(slide_mod, zs_valid, type = "response")
    ## predict values for training set
    train_pred <- slide_mod$fitted.values
  }
  
  ## add SLIDE results to replicate report
  cat("  SLIDE model size: ", ncol(zs_train), "\n")
  
  
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))
}
