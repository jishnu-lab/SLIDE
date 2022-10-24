#' Run partial least-squares discriminant analysis for cross-validation.
#'
#' Creates a partial least-squares discriminant analysis classifier and
#' predicts response values from the validation data. 
#' 
#' @param train_y a vector of numeric values; the response (training set)
#' @param train_x a matrix or data frame of numeric values; the data (training set)
#' @param valid_x a matrix or data frame of numeric values; the data (validation set)
#' @return a list containing the training model and the predicted class probabilities of the response
#' using \code{valid_x} and \code{train_x}
#' @export

runPLSDA <- function(train_y, train_x, valid_x) {
  ## convert to syntactically valid names
  train_y <- make.names(train_y)
  
  ## create cv control
  ctrl <- caret::trainControl(method = "repeatedcv",
                              classProbs = T,
                              savePredictions = T,
                              summaryFunction = caret::twoClassSummary)
  
  ## do plsda
  res <- caret::train(y = train_y,
                      x = train_x,
                      trControl = ctrl,
                      metric = "ROC",
                      method = "pls")
  
  ## using auc to evaluate, predict class probabilities
  valid_pred <- predict(res, newdata = valid_x, type = "raw")
  ## predict for training set
  train_pred <- predict(res, type = "raw")
  
  cat("  PLSDA model size: ", res$bestTune[[1]], "\n")
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))
}