#' Run principal components regression for cross-validation.
#'
#' Creates a linear regression model on the principal components of the data 
#' set and predicts response values for the validation set.
#' 
#' @param train_y a vector of numeric values; the response (training set)
#' @param train_x a matrix or data frame of numeric values; the data (training set)
#' @param valid_x a matrix or data frame of numeric values; the data (validation set)
#' @return a list containing the training model and the predicted values of the response
#' using \code{valid_x} and \code{train_x}
#' @export

runPCR <- function(train_y, train_x, valid_x) {
  ## fit pcr model
  res <- pls::pcr(train_y ~ train_x, validation = "CV", segments = 5)
  ## get maximum number of components
  #n_comp <- min(nrow(train_x)-1,  floor(0.05*ncol(train_x)))
  n_comp <- res$ncomp
  ## predict values for validation set
  valid_pred <- predict(res, comps = n_comp, newdata = valid_x, type = "response")
  ## predict values for training set
  train_pred <- predict(res, comps = n_comp, type = "response")
  
  cat("  PCR model size: ", n_comp, "\n")
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))
}