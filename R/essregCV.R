#' Cross-Validation For Essential Regression.
#'
#' Perform k-fold cross-validation for essential regression to select.
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param k an integer for number of folds to use in cross-validation
#' @param y a response vector of dimension \eqn{n}
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param delta \eqn{\delta}, a numerical constant used for thresholding
#' @param std_cv, a boolean flag of whether perform scaling x in CV and scaling Y
#' @param std_y a boolean flag of wether z score Y or not in the PlainER function
#' @param thresh_fdr a numerical constant used for thresholding the correlation matrix to
#' control the false discovery rate, default is 0.2
#' @param permute a boolean flag indicating whether to permute the response (\eqn{y}) during training
#' @param eval_type a string indicating what metric to use for evaluating model performance
#' (can be "auc" - area under ROC curve, or "corr" - for spearman correlation)
#' @param lambda \eqn{\lambda}, a numerical constant used in thresholding
#' @param rep_cv number of replicates for cross-validation
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @param rep an integer indicating the replicate number
#' @param out_path path for saving output
#' @param benchmark a boolean flag that decides if running lasso, pcr, pclr, plsr and plsda.
#' @param run_lasso a boolean flag that decides if running lasso for comparison
#' @return An object of class \sQuote{data.frame}
#' @export

essregCV <- function(k = 5, y, x, delta, std_cv, std_y, thresh_fdr = 0.2, lambda = 0.1,
                     rep_cv = 50, alpha_level = 0.05, permute = T,y_levels = NULL,
                     eval_type, out_path, rep, benchmark = F, run_lasso = F) {

  #get raw_x and get scaled x ###########################
  raw_x <- x
  x <- scale(x, T, T)
  if (std_cv != std_y ) {stop("std_cv and std_y should be the same...\n")}

  if (eval_type == "auc") {
    lasso_fam <- "binomial"
    y_factor <- T
  } else { ## if evaluating with correlation, treat y as continuous (regardless of truth)
    lasso_fam <- "gaussian"
    y_factor <- F
  }

  ## create output directory
  new_dir <- paste0(out_path, "replicate", rep, "/")
  dir.create(file.path(new_dir), showWarnings = F, recursive = T)

  #################################################################
  ##                  Splitting Data Into Folds                  ##
  #################################################################
  zero_in <- TRUE
  while (zero_in) {

    if (y_factor) {
      group_inds <- caret::createFolds(factor(y), k = k, list = TRUE)
    } else {
      group_inds <- caret::createFolds(y, k = k, list = TRUE)
    }


    y_groups_val <- NULL
    y_groups_train <- NULL
    for (i in 1:length(group_inds)) {
      y_groups_val[[length(y_groups_val) + 1]] <- y[unlist(group_inds[[i]])]
      y_groups_train[[length(y_groups_train) + 1]] <- y[-unlist(group_inds[[i]])]
    }
    # check if we are doing LOOCV
    len <- lapply(y_groups_val, function(x){length(x)})
    if (! 1 %in% len) {
      group_vars_val <- sapply(y_groups_val, sd) ## get standard deviation of responses (validation set)
    }else{
      group_vars_val <- unlist(y_groups_val) ## don't calculate sd if LOOCV
    }
    group_vars_train <- sapply(y_groups_train, sd) ## get standard deviation of responses (training set)
    zero_in_val <- 0 %in% group_vars_val ## see if any are 0s from the sd calculation
    zero_in_train <- 0 %in% group_vars_train
    if (! 1 %in% len) {
      if (zero_in_val || zero_in_train) {
        zero_in <- T
      } else {
        zero_in <- F
      }
    } else{
      # whether we are doing LOOCV or not, let's check the training data.
      if (zero_in_train) {
        zero_in <- T
      } else{
        zero_in <- F
      }
    }

    if (eval_type == "auc") { ## if AUC, go ahead with first split
      if (!zero_in_train) {
        zero_in <- FALSE
      }
    }
  }

  #################################################################
  ##                       Running Methods                       ##
  #################################################################
  methods <- c("plainER")
  if (benchmark){
    if (y_factor) {
      methods <- c(methods, "lasso", "plsda", "pclr")
    } else {
      methods <- c(methods, "lasso", "plsr", "pcr")
    }
  } else if (run_lasso) {
    methods <- c(methods, "lasso")
  }

  if (permute) {
    methods <- c(methods, paste0(methods, "_y"))
  }

  results <- NULL
  ##---------------------------------------------------------------
  ##                  starting cross validation                  --
  ##---------------------------------------------------------------
  for (i in 1:k) { ## loop through folds
    cat("FOLD ", i, ". . . . \n")
    valid_ind <- group_inds[[i]] ## validation indices
    cat("validation indices", valid_ind, "\n")
    train_y_raw <- y[-valid_ind] ## training y's
    valid_y_raw <- y[valid_ind] ## validation y's
    train_x_raw <- raw_x[-valid_ind, ]
    valid_x_raw <- matrix(raw_x[valid_ind, ], ncol = ncol(x))

    # if we are doing z-scoring X within CV and z-scoring Y
    if (std_cv) {
      stands <- standCV(train_y = train_y_raw,
                        train_x = train_x_raw,
                        valid_y = valid_y_raw,
                        valid_x = valid_x_raw)

      train_x_std <- stands$train_x
      train_y <- stands$train_y
      valid_x_std <- stands$valid_x
      valid_y <- stands$valid_y

      # if we are z-scoring X outside of CV and not z-scoring Y
    } else {
      train_x_std <- x[-valid_ind, ] ## training x's
      valid_x_std <- matrix(x[valid_ind, ], ncol = ncol(x)) ## validation x's
      train_y <- train_y_raw
      valid_y <- valid_y_raw
    }


    ## rename columns
    colnames(train_x_std) <- colnames(valid_x_std) <- colnames(x)

    ## permute y's
    # try sampling differently
    # perm_ind = unlist(caret::createFolds(train_y, k = length(train_y) - 1))
    perm_ind <- sample(1:nrow(train_x_std))
    # note that if std_cv == FALSE, train_y_perm == train_y_perm_raw
    train_y_perm <- train_y[perm_ind]
    train_y_perm_raw <- train_y_raw[perm_ind]

    ## get labels if factor
    if (y_factor) {
      train_y_labs <- factor(train_y_raw, levels = y_levels)
      train_y_labs_perm <- factor(train_y_perm_raw, levels = y_levels)
      valid_y_labs <- factor(valid_y_raw, levels = y_levels)
    }

    ##----------------------------------------------------------------
    ##                   loop through all methods                   --
    ##----------------------------------------------------------------
    for (j in 1:length(methods)) { ## loop through methods
      method_j <- methods[j]
      cat("CV for ", method_j, ". . . \n")

      ## determine which y to train with
      if (grepl(x = method_j, pattern = "_y", fixed = TRUE)) { ## if doing y permutation
        if (y_factor) { ## if y is a factor, use permuted y labels
          cat("        using permuted y labels \n")
          use_y_train <- train_y_labs_perm
        } else { ## if y is not a factor, use permuted y values
          cat("        using permuted y values \n")
          use_y_train <- train_y_perm
        }
        use_y_train_ER <- train_y_perm_raw ## need non-standardized continuous y for plainER
      } else { ## if not doing y permutation
        if (y_factor) { ## if y is a factor, use true y labels
          cat("        using true y labels \n")
          use_y_train <- train_y_labs
        } else { ## if y is not a factor, use true y values
          cat("        using true y values \n")
          use_y_train <- train_y
        }
        use_y_train_ER <- train_y_raw ## need non-standardized continuous y for plainER
      }


      ##-------------
      ##  plainER   -
      ##-------------
      if (grepl(x = method_j, pattern = "plainER", fixed = TRUE)) { ## plain essential regression, predict with all Zs

        res <- plainER(y = use_y_train_ER,
                       x = train_x_raw,
                       x_std = train_x_std,
                       std_y = std_y,
                       sigma = NULL,
                       delta = delta,
                       lambda = lambda,
                       thresh_fdr = thresh_fdr,
                       rep_cv = rep_cv,
                       alpha_level = alpha_level)

        if (is.null(res)) {
          return (NULL)
        }
        pred_all_betas <- res$pred$er_predictor
        pred_vals <- valid_x_std %*% pred_all_betas
        #pred_vals <- t((t(beta_valid) - centers_y) / scales_y)

        ##----------
        ##  plsr   -
        ##----------
      } else if (grepl(x = method_j, pattern = "plsr", fixed = TRUE)) { ## PLSR
        res <- pls::plsr(use_y_train ~ train_x_std, validation = "CV", segments = 5)
        n_comp <- pls::selectNcomp(res, method = "randomization")
        n_comp <- ifelse(n_comp == 0, 1, n_comp)
        pred_vals <- predict(res, comps = n_comp, newdata = valid_x_std, type = "response")
        #pred_vals <- t((t(pred_vals) - centers_y) / scales_y)

        ##---------
        ##  pcr   -
        ##---------
      } else if (grepl(x = method_j, pattern = "pcr", fixed = TRUE)) { ## PCR
        res <- pls::pcr(use_y_train ~ train_x_std, validation = "CV", segments = 5)
        n_comp <- pls::selectNcomp(res, method = "randomization")
        n_comp <- ifelse(n_comp == 0, 1, n_comp)
        pred_vals <- predict(res, comps = n_comp, newdata = valid_x_std, type = "response")
        #pred_vals <- t((t(pred_vals) - centers_y) / scales_y)

        ##-----------
        ##  plsda   -
        ##-----------
      } else if (grepl(x = method_j, pattern = "plsda", fixed = TRUE)) { ## PLS-DA
        ## convert to syntactically valid names
        use_y_train_syn <- make.names(use_y_train)
        use_y_valid_syn <- make.names(valid_y_labs)
        ctrl <- caret::trainControl(method = "repeatedcv",
                                    classProbs = T,
                                    savePredictions = T,
                                    summaryFunction = caret::twoClassSummary)
        res <- caret::train(y = use_y_train_syn,
                            x = train_x_std,
                            trControl = ctrl,
                            metric = "ROC",
                            method = "pls")

        pred_vals <- predict(res, newdata = valid_x_std, type = "prob")

        ##----------
        ##  pclr   -
        ##----------
      } else if (grepl(x = method_j, pattern = "pclr", fixed = TRUE)) { ## pc logistic regression
        # res <- randomForest::randomForest(use_y_train ~ .,
        #                                   data = train_x_std)
        # pred_vals <- predict(res, newdata = valid_x_std, type = "prob")
        #pred_vals <- t((t(pred_vals) - centers_y) / scales_y)
        princ_comps <- stats::prcomp(x = train_x_std, center = FALSE, scale = FALSE)
        # cumul_props <- summary(princ_comps)$importance[3, ]
        # num_pcs <- min(which(cumul_props > 0.9))
        # num_pcs <- ifelse(num_pcs >= nrow(train_x_std), nrow(train_x_std) - 1, num_pcs) ## if too many, then just use max number to avoid rank deficiency
        num_pcs <- nrow(train_x_std) - 1 ## use maximum number of PCs
        train_pcs <- princ_comps$x[, 1:min(num_pcs, ncol(princ_comps$x))] ## get PCs for regression
        valid_pcs <- scale(valid_x_std, princ_comps$center, princ_comps$scale) %*% princ_comps$rotation ## project validation set to PC space
        res <- stats::glm(as.numeric(as.character(use_y_train)) ~ ., data = as.data.frame(train_pcs), family = "binomial") ## make model
        pred_vals <- predict(res, newdata = as.data.frame(valid_pcs), type = "response") ## predict validation set values
      } else if (grepl(x = method_j, pattern = "lasso", fixed = TRUE)) { ## lasso for comparison

          #if ((nrow(train_x_std) / 10) < 3) { ## sample size too small
          if (k == nrow(x)) { #LOOCV
            cat("\n using LOOCV for lasso\n")

            # check if we're going to have a problem with cv.glmnet
            if (min(table(as.factor(use_y_train))) <= 2) {
              # we will have an error, so don't do cross val in glmnet
              cat("\ncan't do cv.glmnet, trying glmnet instead\n")
              res = glmnet::glmnet(x = train_x_std, y = use_y_train,
                                   family = lasso_fam, alpha = 1,
                                   lambda = 1,
                                   standardize = F)

              # in lieu of refactoring this entire thing, we are just going
              # to add our lambda value to the res list
              res$lambda.min = res$lambda
            } else {
              res <- glmnet::cv.glmnet(train_x_std,
                                       y = as.factor(use_y_train),
                                       alpha = 1,
                                       nfolds = nrow(train_x_std),
                                       standardize = F,
                                       # type.measure = "class",
                                       grouped = F,
                                       family = lasso_fam)
            }
          } else {
            res <- glmnet::cv.glmnet(train_x_std,
                                     use_y_train,
                                     alpha = 1,
                                     nfolds = 5,
                                     standardize = F,
                                     grouped = F,
                                     family = lasso_fam)
          }

          beta_hat <- coef(res, s = res$lambda.min)[-1]
          sub_beta_hat <- which(beta_hat != 0)

          ## if lasso selects no variable, randomly pick 5 features instead
          if (length(sub_beta_hat) == 0) {
            cat("Lasso selects no features - Randomly selecting 5 features. . . \n")

            sub_beta_hat <- ifelse(ncol(train_x_std) > 5, sample(1:ncol(train_x_std), 5), sample(1:ncol(train_x_std)))

            # sub_beta_hat <- sample(1:ncol(train_x_std), 5)
          }
          ## predict and standardize - use linear model rather than glmnet object
          lasso_train <- as.data.frame(train_x_std[, sub_beta_hat])
          lasso_valid <- as.data.frame(valid_x_std[, sub_beta_hat])
          if (dim(lasso_valid)[[2]] ==1 ){
            lasso_valid <- as.data.frame(t(valid_x_std[, sub_beta_hat]))
          }
          colnames(lasso_train) <- colnames(lasso_valid) <- paste0("X", sub_beta_hat)
          if (y_factor) {
            lasso_lm <- stats::glm(use_y_train ~ ., data = lasso_train, family = lasso_fam)
          } else {
            lasso_lm <- stats::lm(use_y_train ~ ., data = lasso_train)
          }
          pred_vals <- predict(lasso_lm, lasso_valid, type = "response")
          # lasso_pred_vals <- t((t(pred_vals) - centers_y) / scales_y)
        }

        if (eval_type == "auc") { ## if using area under roc curve to evaluate model fit
          pred_vals <- as.data.frame(pred_vals)
          if (ncol(pred_vals) > 1) {
            pred_vals <- pred_vals[, 1]
          }
          fold_res <- cbind(method_j, pred_vals, as.numeric(as.character(valid_y_labs)))
          #print(valid_y_labs)
          colnames(fold_res) <- c("method", "pred_vals", "true_vals")
          results <- rbind(results, fold_res)
          #print(results)
        } else { ## if using correlation to evaluate model fit
          fold_res <- cbind(method_j, pred_vals, valid_y)
          colnames(fold_res) <- c("method", "pred_vals", "true_vals")
          results <- rbind(results, fold_res)
        }
      }
    }

    ## set results data frame column names
    results <- as.data.frame(results)
    final_results <- NULL
    if (eval_type == "auc") {
      for (i in 1:length(methods)) {
        method_res <- results %>% dplyr::filter(method == methods[i])
        predicted <- as.numeric(method_res$pred_vals)
        true <- as.numeric(method_res$true_vals)
        print("before prediction")
        predicted[is.na(predicted)] <- median(predicted, na.rm = TRUE)
        method_roc <- ROCR::prediction(predicted, true)
        print("after prediction")
        method_auc <- ROCR::performance(method_roc, "auc")
        method_auc <- method_auc@y.values[[1]]
        if (method_auc < 0.5) { ## if classifier auc is < 0.5, reverse it to be > 0.5
          method_auc <- 1 - method_auc
        }
        method_res <- c("method" = methods[i],
                        "auc" = as.numeric(method_auc))
        final_results <- rbind(final_results, method_res)
      }
    } else {
      for (i in 1:length(methods)) {
        method_res <- results %>%
          dplyr::filter(method == methods[i])
        predicted <- as.numeric(method_res$pred_vals)
        true <- as.numeric(method_res$true_vals)
        method_corr <- cor(predicted, true, method = "spearman")
        method_mse <- sum((predicted - true)^2) / length(predicted)
        method_res <- c("method" = methods[i],
                        "corr" = as.numeric(method_corr),
                        "mse" = as.numeric(method_mse))
        final_results <- rbind(final_results, method_res)
      }
    }

    final_results <- as.data.frame(final_results)
    final_results[, 2] <- as.numeric(final_results[, 2])

    combined_res <- NULL
    combined_res$each_fold <- results
    combined_res$final_corr <- final_results
    saveRDS(combined_res, file = paste0(new_dir, "results.rds"))
    return (final_results)
}
