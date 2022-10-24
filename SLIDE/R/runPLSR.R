#' Run partial least-squares regression for cross-validation.
#'
#' Creates a partial least-squares regression model and
#' predicts response values from the validation data.
#'
#' @param train_y a vector of numeric values; the response (training set)
#' @param train_x a matrix or data frame of numeric values; the data (training set)
#' @param valid_x a matrix or data frame of numeric values; the data (validation set)
#' @return a list containing the training model and the predicted values of the response
#' using \code{valid_x} and \code{train_x}
#' @export



runPLSR <- function(train_y, train_x, valid_x, n_comp = NULL) {
  

  res <- pls::plsr(train_y ~ train_x, validation = "CV", segments = 5)
  ## get max number of components
  n_comp <- res$ncomp
  ## predict values for validation set
  valid_pred <- predict(res, comps = n_comp, newdata = valid_x, type = "response") ## use maximum number of PCs
  ## predict values for training set
  train_pred <- predict(res, comps = n_comp, type = "response") ## use maximum number of PCs
  
  cat("  PLSR model size: ", n_comp, "\n")
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))

  cat("  PLSR model size: ", n_comp, "\n")
  
}
