#' Run principal components logistic regression for cross-validation.
#'
#' Creates a logistic regression model on the principal components of the data 
#' set and predicts response values for the validation set.
#' 
#' @param train_y a vector of numeric values; the response (training set)
#' @param train_x a matrix or data frame of numeric values; the data (training set)
#' @param valid_x a matrix or data frame of numeric values; the data (validation set)
#' @return a list containing the training model and the predicted class probabilities of the response
#' using \code{valid_x} and \code{train_x}
#' @export

runPCLR <- function(train_y, train_x, valid_x) {
  ## get principal components - should already be scaled
  princ_comps <- stats::prcomp(x = train_x, center = FALSE, scale = FALSE)
  ## get number of components - used maximum number
  num_pcs <- min(nrow(train_x)-1,  floor(0.05*ncol(train_x)))
  ## get PCs for regression
  train_pcs <- princ_comps$x[, 1:num_pcs] %>%
    as.data.frame()
  ## project validation set to PC space
  valid_pcs <- valid_x %*% princ_comps$rotation %>%
    as.data.frame()
  ## make regression model
  res <- stats::glm(train_y ~ ., data = train_pcs, family = "binomial")
  
  ## predict validation set class probabilities
  valid_prob <- predict(res, newdata = valid_pcs, type = "response")
  ## get contrasts
  contrast <- stats::contrasts(train_y)
  ## get labels from class probabilities - validation set
  valid_pred <- ifelse(valid_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
  
  ## predict values for training set
  train_prob <- predict(res, type = "response")
  ## get labels from class probabilities - training set
  train_pred <- ifelse(train_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
  
  cat("  PCLR model size: ", ncol(train_pcs), "\n")
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))
}