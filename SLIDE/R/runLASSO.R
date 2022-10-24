#' Run LASSO regression for cross-validation.
#'
#' Runs LASSO to create a regression model and predicts response values 
#' for the validation set.
#' 
#' @param k the number of folds for cross-validation
#' @param train_y a vector of numeric values; the response (training set)
#' @param train_x a matrix or data frame of numeric values; the data (training set)
#' @param valid_x a matrix or data frame of numeric values; the data (validation set)
#' @param lasso_fam a string indicating whether to perform linear or logistic regression
#' @param y_factor a boolean flag; whether the response is categorical or not
#' @return a list containing the training model and the predicted values of the response
#' using \code{valid_x} and \code{train_x}
#' @export

runLASSO <- function(k, train_y, train_x, valid_x, lasso_fam, y_factor) {
  ## train a lasso model
  if ((nrow(train_x) / 10) < 3) { ## sample size too small
    res <- glmnet::cv.glmnet(train_x,
                             train_y,
                             alpha = 1,
                             nfolds = 5,
                             standardize = F,
                             grouped = F,
                             family = lasso_fam)
  } else {
    res <- glmnet::cv.glmnet(train_x,
                             train_y,
                             alpha = 1,
                             nfolds = 10,
                             standardize = F,
                             grouped = F,
                             family = lasso_fam)
  }
  
  ## get nonzero coefficients
  beta_hat <- coef(res, s = res$lambda.min)[-1]
  sub_beta_hat <- which(beta_hat != 0)
  
  ## if lasso selects no variable, randomly pick 5 features instead
  if (length(sub_beta_hat) == 0) {
    cat("Lasso selects no features - Randomly selecting 5 features. . . \n")
    sub_beta_hat <- sample(1:ncol(train_x), 5)
  }
  
  ## predict - use linear model rather than glmnet object
  lasso_train <- as.data.frame(train_x[, sub_beta_hat])
  # if (k == nrow(x)) { ## if LOOCV
  #   lasso_valid <- as.data.frame(matrix(valid_x[, sub_beta_hat], nrow = 1))
  # } else {
  #   lasso_valid <- as.data.frame(valid_x[, sub_beta_hat])
  # }
  lasso_valid <- as.data.frame(valid_x[, sub_beta_hat])
  ## rename columns of data frames
  colnames(lasso_train) <- colnames(lasso_valid) <- paste0("X", sub_beta_hat)
  
  ## use logistic regression if y is categorical
  ## use simple linear regression if y is continuous
  if (y_factor) {
    lasso_lm <- stats::glm(train_y ~ ., data = lasso_train, family = "binomial")
    ## predict values for validation set
    valid_prob <- predict(lasso_lm, lasso_valid, type = "response")
    ## get contrasts
    contrast <- stats::contrasts(train_y)
    ## get labels from class probabilities - validation set
    valid_pred <- ifelse(valid_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
    
    ## predict values for training set
    train_prob <- lasso_lm$fitted.values
    ## get labels from class probabilities - training set
    train_pred <- ifelse(train_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
  } else {
    lasso_lm <- stats::lm(train_y ~ ., data = lasso_train)
    ## predict values for validation set
    valid_pred <- predict(lasso_lm, lasso_valid, type = "response")
    ## predict values for training set
    train_pred <- lasso_lm$fitted.values
  }

  cat("  LASSO model size:", ncol(lasso_train), "\n")
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))
}