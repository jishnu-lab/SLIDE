evalCV <- function(cv_results, methods, eval_type) {
  final_results <- NULL
  if (eval_type == "auc") { ## if using area under ROC curve for evaluations
    for (i in 1:length(methods)) { ## loop through all methods
      ## filter validation results for one method
      valid_res <- cv_results %>%
        dplyr::filter(method == methods[i], 
                      set == "valid") %>%
        dplyr::mutate(pred_vals = as.numeric(pred_vals),
                      true_vals = as.numeric(true_vals))
      ## filter training results for one method
      train_res <- cv_results %>%
        dplyr::filter(method == methods[i],
                      set == "train") %>%
        dplyr::group_by(index) %>%
        dplyr::mutate(pred_vals = as.numeric(pred_vals),
                      true_vals = as.numeric(true_vals))
      ## get mean predicted value for training set across folds
      freq_pred <- train_res %>%
        dplyr::group_by(index, pred_vals) %>% ## group by both index and predicted class
        dplyr::mutate(class_n = n()) ## get predicted class frequencies by index
      freq_pred <- freq_pred %>%
        dplyr::group_by(index) %>%
        dplyr::filter(class_n == max(class_n)) %>% ## only keep rows with the most frequently predicted class
        unique() ## only keep unique rows
      ## remove rows where neither class is more frequently predicted than the other
      freq_pred <- freq_pred %>%
        dplyr::group_by(index) %>%
        dplyr::mutate(index_n = n())
      
      ## observations where one class dominated
      single_freq <- freq_pred %>%
        dplyr::filter(index_n == 1)
      ## observations where classes were predicted equally frequently
      double_freq <- freq_pred %>%
        dplyr::filter(index_n == 2) %>%
        dplyr::select(-pred_vals) %>%
        unique()
      ## randomly select classes
      double_freq$pred_vals <- stats::rbinom(nrow(double_freq), 1, 0.5)
      double_freq <- double_freq %>%
        dplyr::relocate(pred_vals, .after = method)
      ## combine the two 
      train_res <- rbind(single_freq, double_freq)
      
      ## get spearman correlation between predicted and true response values
      ## training set
      train_pred <- as.numeric(train_res$pred_vals)
      train_true <- as.numeric(train_res$true_vals)
      ## validation set
      valid_pred <- as.numeric(valid_res$pred_vals)
      valid_true <- as.numeric(valid_res$true_vals)
      
      ## get area under ROC curve for predicted vs true response values
      train_roc <- ROCR::prediction(train_pred, train_true)
      valid_roc <- ROCR::prediction(valid_pred, valid_true)
      
      ## get AUCs
      train_auc <- ROCR::performance(train_roc, method = "auc")
      train_auc <- train_auc@y.values[[1]]
      valid_auc <- ROCR::performance(valid_roc, method = "auc")
      valid_auc <- valid_auc@y.values[[1]]
      
      ## if classifier auc is < 0.5, reverse it to be > 0.5
      if (train_auc < 0.5) { 
        train_auc <- 1 - train_auc
      }
      if (valid_auc < 0.5) {
        valid_auc <- 1 - valid_auc
      }
      
      ## append method results to final results table
      method_res <- c("method" = methods[i],
                      "train_auc" = as.numeric(train_auc),
                      "valid_auc" = as.numeric(train_auc))
      final_results <- rbind(final_results, method_res)
    }
  } else { ## if using spearman correlation for evaluations
    for (i in 1:length(methods)) { ## loop through all methods
      ## filter validation results for one method
      valid_res <- cv_results %>%
        dplyr::filter(method == methods[i], 
                      set == "valid") %>%
        dplyr::mutate(pred_vals = as.numeric(pred_vals),
                      true_vals = as.numeric(true_vals))
      ## filter training results for one method
      train_res <- cv_results %>%
        dplyr::filter(method == methods[i],
                      set == "train") %>%
        dplyr::group_by(index) %>%
        dplyr::mutate(pred_vals = as.numeric(pred_vals),
                      true_vals = as.numeric(true_vals))
      ## get mean predicted value for training set across folds
      mean_pred <- train_res %>%
        dplyr::group_by(index) %>%
        dplyr::summarise(mean_pred_vals = mean(pred_vals))
      ## create data frame with mean predicted value, true values, and indices
      train_res <- dplyr::left_join(train_res, mean_pred, by = c("index")) %>%
        dplyr::select(c("method", "mean_pred_vals", "true_vals", "index")) %>%
        unique()
      
      ## get spearman correlation between predicted and true response values
      ## validation set
      valid_pred <- as.numeric(valid_res$pred_vals)
      valid_true <- as.numeric(valid_res$true_vals)
      valid_corr <- cor(valid_pred, valid_true, method = "spearman")
      ## training set
      train_pred <- as.numeric(train_res$mean_pred_vals)
      train_true <- as.numeric(train_res$true_vals)
      train_corr <- cor(train_pred, train_true, method = "spearman")
      
      ## concatenate validation and training set results
      method_res <- c("method" = methods[i],
                      "valid_corr" = as.numeric(valid_corr),
                      "train_corr" = as.numeric(train_corr))
      ## training set
      final_results <- rbind(final_results, method_res)
    }
  }
  
  return (final_results)
}