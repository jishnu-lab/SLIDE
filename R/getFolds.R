#' Split data for cross-validation.
#'
#' Randomly split data into \code{k} folds for cross-validation.
#' 
#' @param k an integer; the number of folds
#' @param x a data frame or matrix; the data
#' @param y a data frame or vector; the reponse
#' @return a list of vectors containing the indices of the rows of the data found in each fold
#' @export

getFolds <- function(k, x, y) {
  zero_in <- TRUE
  while (zero_in) { ## if any are 0, redo the fold divisions
    total_num <- nrow(x)
    num_group <- k
    remainder <- total_num %% num_group # get the remainder
    num_per_group <- total_num %/% num_group # get number of indices per fold
    ## partition data into folds
    partition <- rep(num_per_group, num_group) + c(rep(1, remainder), rep(0, num_group - remainder))
    ## extract the indices in each fold
    pre_vec <- sample(1:nrow(x))
    extract <- vector("list", length(partition))
    extract[[1]] <- pre_vec[1:partition[1]]
    for (i in 2:length(partition)) {
      temp_ind <- sum(partition[1:(i - 1)]) + 1
      extract[[i]] <- pre_vec[temp_ind:(temp_ind + partition[i] - 1)]
    }
    group_inds <- extract
    
    ## split the response vector into folds as well
    y_groups_val <- NULL
    y_groups_train <- NULL
    for (i in 1:length(group_inds)) {
      y_groups_val[[length(y_groups_val) + 1]] <- y[unlist(group_inds[[i]])]
      y_groups_train[[length(y_groups_train) + 1]] <- y[-unlist(group_inds[[i]])]
    }
    ## get standard deviation of responses (validation set)
    group_vars_val <- ifelse(k == nrow(x), 1, sapply(y_groups_val, sd)) 
    ## get standard deviation of responses (training set)
    group_vars_train <- sapply(y_groups_train, sd) 
    ## see if any are 0
    zero_in_val <- 0 %in% group_vars_val 
    ## see if any are 0
    zero_in_train <- 0 %in% group_vars_train 
    
    ## want a good mix of valid (don't want a fold with all the same value)
    if (zero_in_val || zero_in_train) {
      cat("no variation in one or more folds . . . resplitting \n")
      zero_in <- T
    } else {
      zero_in <- F
    }
  }
  
  return (group_inds)
}