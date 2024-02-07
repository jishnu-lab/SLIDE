#' Sample Cross Validation.
#'
#' Calculate SLIDE performance using the sample CV framework for a quick evaluation. This is meant to help user to choose which model to use for the more rigorous cross validation. 
#'
#' @param y the response matrix.
#' @param z_matrix the z_matrix, the output of CalcZMatrix function
#' @param SLIDE_res the output of GetTopFeatures function.
#' @param fraction the fraction of samples to choose for sampling.
#' @param condition use auc or corr(spearman) to evaluate performance.
#' @param sampleCV_iter a integer for number of iterations.
#' @param logistic a boolean flag on wether running logistic or lm.
#' @param out_path the output path.
#' @return the mean performance.
#' @export


sampleCV <- function(y, z_matrix, SLIDE_res, fraction = 2/3, condition, sampleCV_iter = 20,  logistic = FALSE, out_path){
  
  # if the number of SLIDE chosen lfs are bigger than 2/3(default number) of the sample numbers, return na
  lf_idx = union(SLIDE_res$marginal_vals, union(SLIDE_res$interaction$p1, SLIDE_res$interaction$p2))
  if (length(lf_idx) > length(sample(nrow(z_matrix), ceiling(fraction * nrow(z_matrix))))) {
    perf = "NA"
    return(perf)
  } 
  
  # get the significant lf and interactors
  sigK <- SLIDE_res$marginal_vals
  sigIn <- as.vector(SLIDE_res$SLIDE_res$interaction_vars)
  
  # get the z_matrix for SLIDE selected LFs and the interaction terms
  IntData <- pairwiseInteractions(sigK,z_matrix)
  Dataint <- IntData$interaction[, sigIn]
  # create the data for lm or logistic regression
  lm_data <- data.frame(y = y, z_matrix[, sigK], Dataint)
  colnames(lm_data)[1] = 'y'
  
  perfs = c()
  for (iter in 1:sampleCV_iter){
    
    # randomly choose idx each round
    sample_idx <- sample(nrow(z_matrix), ceiling(fraction * nrow(z_matrix)))
    train_lm_data <- lm_data[sample_idx, ]
    test_lm_data <- lm_data[-sample_idx, ]
    # get rid of the y in the test data
    y_true <- test_lm_data$y
    test_lm_data <- test_lm_data[, -1]
    if(nrow(train_lm_data) + nrow(test_lm_data) != nrow(z_matrix)){stop("The data split in sampleCV is wrong...")}
    if(length(y_true) != nrow(test_lm_data)){stop("The dimensions of true y and test data doesn't match in samplesCV.")}
    
    # build the model and get y_hat and y_true
    if (logistic == FALSE){ # running linear regression
      model <-  lm(y ~ ., data = train_lm_data)
      y_hat = predict(model, newdata = test_lm_data)
    } else{ # running logistic 
      model <- glm(subset_y ~. , data = train_lm_data, family = binomial())
      y_hat = predict(model, newdata = test_lm_data, type = 'response')
      y_hat = ifelse(y_hat>0.5, 1, 0)
    }
    
    # calculate performance
    if (condition == 'corr'){
      perf = cor(y_hat, y_true, method = 'spearman')
    }else if (condition == 'auc'){
      perf = pROC::auc(response=as.matrix(y_true), predictor=as.matrix(y_hat))
    }else{
      stop("Error: Condition not found in sampleCV. ")
    }
    
    perfs = append(perfs, perf)
  }
  model$samplesCV$performances = perfs
  saveRDS(model, paste0(out_path, "sampleCV_model.RDS"))
  return(mean(perfs))
  
}