#' Run Essential Regression for cross-validation.
#'
#' Runs Essential Regression to create a model and predicts response values
#' for the validation set.
#'
#' @param train_y a vector of numeric values; the response (training set)
#' @param train_x_raw a matrix or data frame of numeric values; the data (training set)
#' @param train_x_std a matrix or data frame of numeric values; the data (training set)
#' @param valid_x a matrix or data frame of numeric values; the data (validation set)
#' @param delta \eqn{\delta}, a numerical value used for thresholding
#' @param lambda \eqn{\lambda}, a numerical constant used in thresholding
#' @param thresh_fdr a numerical constant used for thresholding the correlation matrix to
#' control the false discovery rate, default is 0.2
#' @param rep_cv number of replicates for inner cross-validation of \code{delta}
#' @param alpha_level \eqn{\alpha}, a numerical constant used in confidence interval calculation
#' @return a list containing the training model and the predicted values of the response
#' using \code{valid_x} and using \code{train_x_std}


runEssReg <- function(train_y, train_x_raw, train_x_std, valid_x, delta, lambda, thresh_fdr, rep_cv, alpha_level) {
  ## run ER
  res <- EssReg::plainER(y = train_y,
                         x = train_x_raw,
                         x_std = train_x_std,
                         std_y = FALSE,
                         sigma = NULL,
                         delta = delta,
                         lambda = lambda,
                         thresh_fdr = thresh_fdr,
                         rep_cv = rep_cv,
                         alpha_level = alpha_level)
  pred_all_betas <- res$pred$er_predictor
  ## predict values for validation set
  if(is.null(pred_all_betas)){pred_all_betas <- matrix(rnorm(n=dim(valid_x)[2]),ncol=1)}
  valid_pred <- valid_x %*% pred_all_betas
  ## predict values for training set
  train_pred <- train_x_std %*% pred_all_betas
  ## record model size
  cat("  ER model size: ", res$K, "\n")
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))
}
