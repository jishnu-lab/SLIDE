#' Calculate prediction error.
#'
#' The meat of \code{runSims()}. This function splits the data into a training set
#' and a validation set, trains the methods, and then calculates the error (MSE) using the
#' validation set.
#'
#' @param x a data frame of numeric values; the data
#' @param y a vector or data frame with one column; the response
#' @param spec a numeric constant; the proportion of iterations that a variable must be selected in to be considered significant
#' @param niter an integer; the number of times to run knockoffs
#' @param fdr a numeric constant between 0 and 1; the target false discovery rate for knockoffs
#' @param f_size an integer; the target size for the subset (number of columns in each subset)
#' @param ncore an integer; the number of cores to use for parallel computing
#' @param out_file a string path/file name; where to save the results
#' @return a data frame with 5 rows and 2 columns (method name and method MSE)
#' @export



predError <- function(x, y, spec = 0.3, niter = 1000, fdr = 0.1, f_size = 100, ncore = 64, out_file, k = NULL) {
  ## split into training and validation sets --- use 75% of data for training
  num_valid <- ceiling(nrow(x) * 0.25) ## get size of validation set
  valid_ind <- sample(seq(1, nrow(x)), num_valid) ## sample validation set indices
  ## get validation set
  valid_x <- x[valid_ind, ]
  valid_y <- y[valid_ind]
  ## get training set
  train_x <- x[-valid_ind, ]
  train_y <- y[-valid_ind]

  ## standardize sets
  stds <- standCV(
    valid_x = valid_x,
    valid_y = valid_y,
    train_x = train_x,
    train_y = train_y
  )
  valid_x_std <- stds$valid_x
  valid_y_std <- stds$valid_y
  train_x_std <- stds$train_x
  train_y_std <- stds$train_y

  ## run EssReg
  cat("        EssReg . . . . \n")
  res_essreg <- runEssReg(
    train_y = train_y_std,
    train_x_raw = train_x,
    train_x_std = train_x_std,
    valid_x = valid_x_std,
    delta = 0.1,
    lambda = 1,
    thresh_fdr = 0.2,
    rep_cv = 100,
    alpha_level = 0.05
  )
  res_z_train <- EssReg::predZ(x = train_x_std, er_res = res_essreg$model)
  res_z_valid <- EssReg::predZ(x = valid_x_std, er_res = res_essreg$model)
  colnames(res_z_train) <- colnames(res_z_valid) <- paste0("Z", 1:ncol(res_z_valid))

  # run SLIDE
  cat("        SLIDE . . . . \n")
  res_slide <- runSLIDE(
    train_y = train_y_std,
    valid_y = valid_y_std,
    train_z = res_z_train,
    valid_z = res_z_valid,
    method = 4,
    spec = spec,
    niter = niter,
    fdr = fdr,
    f_size = f_size,
    betas = NULL,
    top_prop = NULL,
    marginals = NULL,
    parallel = T,
    ncore = ncore,
    y_factor = FALSE
  )

  ## run LASSO
  cat("        LASSO . . . . \n")
  res_lasso <- runLASSO(
    k = 0,
    train_y = train_y_std,
    train_x = train_x_std,
    valid_x = valid_x_std,
    lasso_fam = "gaussian",
    y_factor = FALSE
  )

  ## run PLSR
  cat("        PLSR . . . . \n")
  random_num <- rnorm(n=1)
  
  saveRDS(file = paste0("/ix/djishnu/Javad/simulationResult/dbug","/",random_num,".rds"), list(
    train_y_std = train_y_std,
    train_x_std = train_x_std,
    valid_x_std = valid_x_std
  ))


  res_plsr <- runPLSR(
    train_y = train_y_std,
    train_x = train_x_std,
    valid_x = valid_x_std, 
    n_comp = 100
  )

  ## run PCR
  cat("        PCR . . . . \n")
  res_pcr <- runPCR(
    train_y = train_y_std,
    train_x = train_x_std,
    valid_x = valid_x_std
  )

  #save(res_essreg, res_slide, res_lasso, res_plsr, res_pcr, file = out_file)

  ## evaluate model performances (MSE)
  ess_reg_mse <- mean((valid_y_std - res_essreg$valid_pred)^2)
  slide_mse <- mean((valid_y_std - res_slide$valid_pred)^2)
  lasso_mse <- mean((valid_y_std - res_lasso$valid_pred)^2)
  plsr_mse <- mean((valid_y_std - res_plsr$valid_pred)^2)
  pcr_mse <- mean((valid_y_std - res_pcr$valid_pred)^2)
  ## combine into data frame
  mse_results <- data.frame(
    "method" = c("ER", "SLIDE", "LASSO", "PLSR", "PCR"),
    "MSE" = c(ess_reg_mse, slide_mse, lasso_mse, plsr_mse, pcr_mse)
  )
  return(mse_results)
}
