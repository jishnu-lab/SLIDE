#' Benchmarking methods with cross-validation.
#'
#' Runs one replicate of a cross-validation regime for method performance comparison. Currently supports
#' plainER, SLIDE, LASSO regression, PHATE regression. Additional methods compared include principal components
#' regression and partial least squares regression for continuous response and principal components logistic regression
#' and partial least squares discriminant analysis for categorical response.
#'
#' @param rep an integer; the replicate number
#' @param k an integer; the number of folds for cross-validation
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param x a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param z a matrix or data frame; latent factors identified using \code{x}
#' @param y_order a vector; the levels of \code{y} if it is ordinal
#' @param std a boolean flag; whether standardization is performed outside (FALSE) or within cross-validation splitting (TRUE)
#' @param eval_type a string either "corr" or "auc"; the method by which to evaluate performance
#' @param niter an integer; the number of times to run knockoffs
#' @param spec a numeric constant; the proportion of iterations that a variable must be selected in to be considered significant
#' @param fdr a numeric constant between 0 and 1; the target false discovery rate for knockoffs
#' @param f_size an integer; the target size for the subset (number of columns in each subset)
#' @param method an integer (1, 2, 3, 4); the selection method for determining marginal/interaction significance
#' @param delta \eqn{\delta}, a numerical value used for thresholding
#' @param lambda \eqn{\lambda}, a numerical constant used in thresholding
#' @param thresh_fdr a numerical constant used for thresholding the correlation matrix to
#' control the false discovery rate, default is 0.2
#' @param rep_cv number of replicates for inner cross-validation of \code{delta}
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @param ncore an integer; the number of cores to use for parallel computing
#' @param out_path a string path; where to save output
#' @return a matrix of model performances
#' @export

benchCV <- function(rep, k = 5, y, x, z, y_order = NULL, std = TRUE, eval_type, niter = 100, spec = 0.1,
                    fdr = 0.1, f_size = 100, method = 4, delta, lambda = 0.1, thresh_fdr = 0.2,
                    rep_cv = 50, alpha_level = 0.05, parallel = T, ncore = 10, out_path) {
  sprintf("the dimension of xs  are %s", dim(x)[1])

  raw_x       <- x
  ## if standardizing outside cross-validation
  if (!std) {
    cat("Standardizing data outside cross-validation . . . \n")
    x         <- scale(raw_x, T, T)
  }

  ## setting up settings for evaluation and model building
  if (eval_type == "auc") {
    if (length(unique(as.vector(y))) > 2) { ## check y is binary
      stop("y must be binary - please re-factor")
    }
    lasso_fam <- "binomial"
    y_factor  <- T
  } else { ## if evaluating with correlation, treat y as continuous (regardless of truth)
    lasso_fam <- "gaussian"
    y_factor  <- F
  }

  ## divide into folds
  group_inds <- getFolds(
    k = k,
    x = x,
    y = y
  )

  ## methods list - these will be performed in the method comparison
  methods   <- c("plainER", "plainER_y", "SLIDE", "lasso","plsr","pcr", "lasso_y", "phate")
  # methods <- c("lasso", "lasso_y", "phate")

  ## add methods depending upon response type
  if (y_factor) {
    methods <- c(methods, "plsda", "pclr")
  } else {
    methods <- c(methods, "plsr", "pcr")
  }

  ## initialization
  results <- NULL

  ## set up timing of entire replicate
  tictoc::tic.clearlog()
  tictoc::tic("entire replicate")

  for (i in 1:k) { ## loop through folds
    cat("######## FOLD ", i, "######## \n")
    tictoc::tic(paste0("FOLD ", i)) ## timing for fold

    ## separate data into training and validation sets
    valid_ind   <- group_inds[[i]] ## validation indices
    train_ind   <- seq(1, nrow(x))[-valid_ind] ## training indices

    ## add validation set indices to replicate report
    cat("validation set indices: ", paste0(valid_ind, " ", collapse = ""), "\n")

    train_y     <- y[-valid_ind] ## training y's
    valid_y     <- y[valid_ind] ## validation y's

    train_x     <- x[-valid_ind, ] ## training x's --- already standardized if std_cv = F
    train_x_raw <- raw_x[-valid_ind, ] ## training x's
    valid_x     <- matrix(x[valid_ind, ], ncol = ncol(x)) ## validation x's
    valid_x_raw <- matrix(raw_x[valid_ind, ], ncol = ncol(raw_x)) ## validation x's

    train_z     <- z[-valid_ind, ] ## training z's
    valid_z     <- z[valid_ind, ] ## validation z's

    ## standardize the training and validation sets
    stands <- EssReg::standCV(
      train_y = train_y,
      train_x = train_x_raw,
      valid_y = valid_y,
      valid_x = valid_x_raw
    )

    if (std) { ## use inner standardized x's
      cat("Standardizing data within cross-validation . . . \n")
      train_x_std <- stands$train_x
      valid_x_std <- stands$valid_x
      train_y_std <- stands$train_y
      valid_y_std <- stands$valid_y
    } else { ## use outer standardized x's
      train_x_std <- train_x
      valid_x_std <- valid_x
      train_y_std <- train_y
      valid_y_std <- valid_y
    }

    ## rename columns
    colnames(train_x_std) <- colnames(valid_x_std) <- colnames(train_x_raw) <- colnames(valid_x_raw) <- colnames(x)

    perm_ind <- sample(1:nrow(train_x_std))
    train_y_perm <- train_y[perm_ind] ## not standardized to maintain labels
    train_y_std_perm <- train_y_std[perm_ind]

    ## get labels if factor
    if (y_factor) {
      train_y_labs <- factor(train_y, levels = y_order)
      train_y_labs_perm <- factor(train_y_perm, levels = y_order)
      valid_y_labs <- factor(valid_y, levels = y_order)
    }

    for (j in 1:length(methods)) { ## loop through methods
      method_j <- methods[j]
      cat("CV for ", method_j, ". . . \n")
      tictoc::tic(paste0("fold ", i, " ", method_j)) ## timing for each method

      ## determine which y to train with
      if (grepl(x = method_j, pattern = "_y", fixed = TRUE)) { ## if doing y permutation
        if (y_factor) { ## if y is a factor, use permuted y labels
         cat("        using permuted y labels \n")
         use_y_train <- train_y_labs_perm
      } else { ## if y is not a factor, use permuted y values
         cat("        using permuted y values \n")
         use_y_train <- train_y_std_perm
      }
        use_y_train_nonstd <- train_y_perm ## need non-standardized continuous y for plainER
       } else { ## if not doing y permutation
         if (y_factor) { ## if y is a factor, use true y labels
          cat("        using true y labels \n")
          use_y_train <- train_y_labs
        } else { ## if y is not a factor, use true y values
          cat("        using true y values \n")
          use_y_train <- train_y_std
        }
        use_y_train_nonstd <- train_y ## need non-standardized continuous y for plainER
      }
      
      
      
      
      
      if (grepl(x = method_j, pattern = "plainER_y", fixed = TRUE)) { ## plain essential regression, predict with all Zs
        cat("        using permuted y values on ER \n")
        result <- runEssReg(
          train_y = use_y_train_nonstd[perm_ind],
          train_x_raw = train_x_raw,
          train_x_std = train_x_std,
          valid_x = valid_x_std,
          delta = delta,
          lambda = lambda,
          thresh_fdr = thresh_fdr,
          rep_cv = rep_cv,
          alpha_level = alpha_level
        )}
      
        if (grepl(x = method_j, pattern = "lasso_y", fixed = TRUE)) { ## LASSO REGRESSION
          
          cat("using permuted y values on lasso \n")
          print(paste0(k, "\n"))
          print(paste0(dim(use_y_train)[1], "\n"))
          print(paste0(dim(train_x_std)[1], "\n"))
          
          result <- runLASSO(
            k = k,
            train_y = use_y_train[perm_ind],
            train_x = train_x_std,
            valid_x = valid_x_std,
            lasso_fam = lasso_fam,
            y_factor = y_factor
          )}
      

      if (grepl(x = method_j, pattern = "plainER", fixed = TRUE)) { ## plain essential regression, predict with all Zs
        result <- runEssReg(
          train_y = use_y_train,
          train_x_raw = train_x_raw,
          train_x_std = train_x_std,
          valid_x = valid_x_std,
          delta = delta,
          lambda = lambda,
          thresh_fdr = thresh_fdr,
          rep_cv = rep_cv,
          alpha_level = alpha_level
        )
      } else if (grepl(x = method_j, pattern = "SLIDE", fixed = TRUE)) { ## run SLIDE
        result <- runSLIDEBeta(
          train_y = use_y_train,
          train_z = train_z,
          valid_z = valid_z,
          method = method,
          niter = niter,
          spec = spec,
          marginals = NULL,
          parallel = parallel,
          f_size = f_size,
          y_factor = y_factor,
          betas = NULL,
          top_prop = NULL,
          ncore = ncore,
          fdr = fdr
        )
      } else if (grepl(x = method_j, pattern = "plsr", fixed = TRUE)) { ## PLSR
        cat("Running PLSR ... \n")
        result <- runPLSR(
          train_y = use_y_train,
          train_x = train_x_std,
          valid_x = valid_x_std
        )
      } else if (grepl(x = method_j, pattern = "pcr", fixed = TRUE)) { ## PCR
        cat("Running PCR ... \n")
        result <- runPCR(
          train_y = use_y_train,
          train_x = train_x_std,
          valid_x = valid_x_std
        )
      } else if (grepl(x = method_j, pattern = "plsda", fixed = TRUE)) { ## PLS-DA
        cat("Running PLSDA ... \n")
        result <- runPLSDA(
          train_y = use_y_train,
          train_x = train_x_std,
          valid_x = valid_x_std
        )
      } else if (grepl(x = method_j, pattern = "pclr", fixed = TRUE)) { ## PRINCIPAL COMPONENTS LOGISTIC REGRESSION
        result <- runPCLR(
          train_y = use_y_train,
          train_x = train_x_std,
          valid_x = valid_x_std
        )
      } else if (grepl(x = method_j, pattern = "phate", fixed = TRUE)) { ## PHATE
        cat("Running Phate ... /n")
        result <- runPHATE(
          train_y = use_y_train,
          train_x = train_x_std,
          valid_x = valid_x_std,
          y_factor = y_factor
        )
      } else if (grepl(x = method_j, pattern = "lasso", fixed = TRUE)) { ## LASSO REGRESSION

        cat("Running lasso ... \n")
        print(paste0(k, "\n"))
        print(paste0(dim(use_y_train)[1], "\n"))
        print(paste0(dim(train_x_std)[1], "\n"))

        result <- runLASSO(
          k = k,
          train_y = use_y_train,
          train_x = train_x_std,
          valid_x = valid_x_std,
          lasso_fam = lasso_fam,
          y_factor = y_factor
        )
      } else {
        stop("Invalid method specified. Please try again. \n")
      }

      res                  <- result$model
      valid_pred           <- result$valid_pred
      train_pred           <- result$train_pred

      if (eval_type == "auc") { ## if using area under roc curve to evaluate model fit
        valid_pred         <- as.data.frame(valid_pred)
        train_pred         <- as.data.frame(train_pred)
        fold_res           <- cbind(method_j, valid_pred, valid_y_labs, valid_ind, "valid") ## organize validation results
        train_res          <- cbind(method_j, train_pred, use_y_train, train_ind, "train") ## organize training results
        colnames(fold_res) <- colnames(train_res) <- c("method", "pred_vals", "true_vals", "index", "set")
        fold_res           <- rbind(fold_res, train_res) ## join together all results
      } else { ## if using correlation to evaluate model fit
        fold_res           <- cbind(method_j, valid_pred, valid_y, valid_ind, "valid") ## organize validation results
        train_res          <- cbind(method_j, train_pred, use_y_train, train_ind, "train") ## organize training results
        colnames(fold_res) <- colnames(train_res) <- c("method", "pred_vals", "true_vals", "index", "set")
        fold_res           <- rbind(fold_res, train_res) ## join together all results
      }

      ## save fold results
      saveRDS(list(
        train_y = use_y_train,
        valid_y = valid_y,
        train_x = train_x_std,
        valid_x = valid_x_std,
        train_z = train_z,
        valid_z = valid_z,
        results = fold_res,
        model = res,
        valid_pred = valid_pred,
        train_pred = train_pred
      ),
      file = paste0(out_path, "replicate", rep, "/fold", i, "_", method_j, ".rds")
      )

      results <- rbind(results, fold_res)
      tictoc::toc()
      cat("\n")
    }
    tictoc::toc()
    cat("\n")
  }

  ## set results data frame column names
  results <- as.data.frame(results)
  saveRDS(results, file = paste0(out_path, "replicate", rep, "/predicted_values.rds"))

  ## evaluate the performance of each method
  final_results <- evalCV(
    cv_results = results,
    methods = methods,
    eval_type = eval_type
  )

  ## end timing for replicate
  tictoc::toc()
  cat("\n")

  ## cast correlation/AUC to numeric and save results
  final_results <- as.data.frame(final_results)
  final_results[, 2] <- as.numeric(final_results[, 2])
  final_results[, 3] <- as.numeric(final_results[, 3])
  saveRDS(final_results, file = paste0(out_path, "replicate", rep, "/model_evaluations.rds"))

  return(final_results)
}
