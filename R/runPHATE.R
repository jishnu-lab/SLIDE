#' Run PHATE regression for cross-validation.
#'
#' Creates a regression model on the embeddings of the data 
#' set found from PHATE and predicts response values for the validation set.
#' 
#' @param train_y a vector of numeric values; the response (training set)
#' @param train_x a matrix or data frame of numeric values; the data (training set)
#' @param valid_x a matrix or data frame of numeric values; the data (validation set)
#' @param y_factor a boolean flag; whether the response is categorical or not
#' @return a list containing the training model and the predicted values of the response
#' using \code{valid_x} and \code{train_x}
#' @export

runPHATE <- function(train_y, train_x, valid_x, y_factor) {
  ## make phate data frame
  phate_df <- rbind(train_x, valid_x)
  ## use max number of components - 1
  phate_mod <- phateR::phate(phate_df, n_pca = as.integer(nrow(phate_df) - 1)) 
  phate_embed <- phate_mod$embedding %>%
    as.data.frame()
  
  ## get training embeddings
  phate_train <- phate_embed[1:nrow(train_x), ]
  ## get validation embeddings
  phate_valid <- phate_embed[(nrow(train_x) + 1):nrow(phate_embed), ]
  
  ## if y is continuous, use simple linear regression
  ## if else, use logistic regression
  if (y_factor) {
    res <- stats::glm(train_y ~ ., data = phate_train, family = "binomial")
    ## predict validation set class probabilities
    valid_prob <- predict(res, newdata = phate_valid, type = "response")
    ## get contrasts
    contrast <- stats::contrasts(train_y)
    ## get labels from class probabilities - validation set
    valid_pred <- ifelse(valid_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
    
    ## predict values for training set
    train_prob <- predict(res, type = "response")
    ## get labels from class probabilities - training set
    train_pred <- ifelse(train_prob > 0.5, rownames(contrast)[2], rownames(contrast)[1])
  } else {
    res <- stats::lm(train_y ~ ., data = phate_train)
    ## predict values for validation set
    valid_pred <- predict(res, phate_valid, type = "response")
    ## predict values for training set
    train_pred <- res$fitted.values
  }
  
  cat("  PHATE model size: ", ncol(phate_train), "\n")
  return (list("model" = res,
               "valid_pred" = valid_pred,
               "train_pred" = train_pred))
}