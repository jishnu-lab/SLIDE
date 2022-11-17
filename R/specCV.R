specCV <- function(k = 5, y, x, z, out_path, rep, parallel = T, y_order = NULL, method = 4,
                   ncore = 10, eval_type, niter = 100, fdr = 0.1, f_size = 100, spec = 0.1) {
  ## create output directory
  new_dir <- paste0(out_path, "replicate", rep, "/")
  dir.create(file.path(new_dir), showWarnings = F, recursive = T)
  
  ## set up for auc/corr evaluations
  if (eval_type == "auc") {
    if (length(unique(as.vector(y))) > 2) { ## check y is binary
      stop("y must be binary - please re-factor")
    }
    lasso_fam <- "binomial"
    y_factor <- T
  } else { ## if evaluating with correlation, treat y as continuous (regardless of truth)
    lasso_fam <- "gaussian"
    y_factor <- F
  }
  
  ## divide into folds
  zero_in <- TRUE
  while (zero_in) { ## if any are 0, redo the fold divisions
    ## first part is partition()
    total_num <- nrow(x)
    num_group <- k
    remainder <- total_num %% num_group # get the remainder
    num_per_group <- total_num %/% num_group
    partition <- rep(num_per_group, num_group) + c(rep(1, remainder), rep(0, num_group - remainder))
    ## second part is extract()
    pre_vec <- sample(1:nrow(x))
    extract <- vector("list", length(partition))
    extract[[1]] <- pre_vec[1:partition[1]]
    for (i in 2:length(partition)) {
      temp_ind <- sum(partition[1:(i - 1)]) + 1
      extract[[i]] <- pre_vec[temp_ind:(temp_ind + partition[i] - 1)]
    }
    group_inds <- extract
    
    y_groups_val <- NULL
    y_groups_train <- NULL
    for (i in 1:length(group_inds)) {
      y_groups_val[[length(y_groups_val) + 1]] <- y[unlist(group_inds[[i]])]
      y_groups_train[[length(y_groups_train) + 1]] <- y[-unlist(group_inds[[i]])]
    }
    group_vars_val <- ifelse(k == nrow(x), 1, sapply(y_groups_val, sd)) ## get standard deviation of responses (validation set)
    group_vars_train <- sapply(y_groups_train, sd) ## get standard deviation of responses (training set)
    zero_in_val <- 0 %in% group_vars_val ## see if any are 0
    zero_in_train <- 0 %in% group_vars_train ## see if any are 0
    if (zero_in_val || zero_in_train) {
      zero_in <- T
    } else {
      zero_in <- F
    }
  }
  
  ## spec list
  specs <- c(0.0, 0.1, 0.3, 0.5, "elbow", "ER")
  
  ## initialization
  results <- NULL
  pick_num <- -1
  
  ## begin CV
  tictoc::tic.clearlog()
  tictoc::tic("entire replicate")
  for (i in 1:k) { ## loop through folds
    cat("######## FOLD ", i, "######## \n")
    tictoc::tic(paste0("FOLD ", i))
    valid_ind <- group_inds[[i]] ## validation indices
    ## add validation set indices to replicate report
    cat("validation set indices: ", paste0(valid_ind, " ", collapse = ""), "\n")
    
    train_y <- y[-valid_ind] ## training y's
    valid_y <- y[valid_ind] ## validation y's
    train_x <- x[-valid_ind, ] ## training x's
    valid_x <- matrix(x[valid_ind, ], ncol = ncol(x)) ## validation x's
    train_z <- z[-valid_ind, ] ## training z's
    valid_z <- z[valid_ind, ] ## validation z's
    
    ## standardize sets
    stands <- EssReg::standCV(train_y = train_y,
                              train_x = train_x,
                              valid_y = valid_y,
                              valid_x = valid_x)
    
    train_x_std <- stands$train_x
    train_y_std <- stands$train_y
    valid_x_std <- stands$valid_x
    valid_y_std <- stands$valid_y
    
    centers_y <- attr(train_y_std, "scaled:center")
    scales_y <- attr(train_y_std, "scaled:scale")
    
    ## rename columns
    colnames(train_x_std) <- colnames(valid_x_std) <- colnames(x)
    
    ## permute y's
    perm_ind <- sample(1:nrow(train_x_std))
    train_y_std_perm <- train_y_std[perm_ind]
    train_y_perm <- train_y[perm_ind]
    
    ## get labels if factor and determine which y to train with
    if (y_factor) {
      train_y_labs <- factor(train_y, levels = y_order)
      train_y_labs_perm <- factor(train_y_perm, levels = y_order)
      valid_y_labs <- factor(valid_y, levels = y_order)
      cat("        using true y labels \n")
      use_y_train <- train_y_labs
    } else {
      cat("        using true y values \n")
      use_y_train <- train_y
    }
    
    for (j in 1:length(specs)) { ## loop through different specs
      spec <- specs[j]
      if (spec == "elbow" || spec == "ER") {
        elbow <- TRUE
      } else {
        spec <- as.numeric(spec)
        elbow <- FALSE
      }
      
      if (spec == "ER") {
        ## run ER
        res <- EssReg::plainER(y = train_y,
                               x = train_x,
                               sigma = NULL,
                               delta = 0.06,
                               lambda = 1,
                               thresh_fdr = 0.2,
                               rep_cv = 100,
                               alpha_level = 0.05)
        ## predict values for validation set
        pred_all_betas <- res$pred$er_predictor
        pred_vals <- valid_x %*% pred_all_betas
        ## record model size
        cat("  ER model size: ", res$K, "\n")
      } else {
        cat("CV for spec ", spec, ". . . \n")
        ############# SLIDE ######################################################
        cat(paste0("        selecting Zs by method 4 with spec ", spec, "\n"))
        res <- SLIDE(z = train_z,
                     y = use_y_train,
                     method = method,
                     niter = niter,
                     spec = spec,
                     elbow = elbow,
                     marginals = NULL,
                     parallel = parallel,
                     f_size = f_size,
                     betas = NULL,
                     top_prop = NULL,
                     ncore = ncore,
                     do_interacts = TRUE,
                     fdr = fdr,
                     out_path = out_path)
        if (length(res$marginal_vars) == 0) { ## if no sig margs, select 5 random ones (like in lasso)
          cat("SLIDE selects no features - Randomly selecting 5 features. . . \n")
          sig_margs <- sample(1:ncol(train_z), 5)
          res <- SLIDE(z = train_z,
                       y = use_y_train,
                       method = method,
                       niter = niter,
                       spec = spec,
                       elbow = elbow,
                       marginals = sig_margs,
                       parallel = parallel,
                       f_size = f_size,
                       betas = NULL,
                       top_prop = NULL,
                       ncore = ncore,
                       do_interacts = TRUE,
                       fdr = fdr,
                       out_path = out_path)
        }
        ## save
       
        
        ## get training Zs (margs + interacts)
        if (!is.null(res$upsilon)) {  ## interactions found
          zs_train <- res$upsilon
          ## get validation interactions
          zs_ints_valid <- matrix(nrow = nrow(valid_x), ncol = 0)
          for (t in 1:ncol(zs_train)) {
            marg_var <- colnames(zs_train)[t]
            only_ind <- gsub("upsZ", "", marg_var)
            marg_var <- gsub("ups", "", marg_var)
            
            inter_union <- interUnion(marginal_vars = marg_var, z = valid_z)
            selected_interactions <- res$interaction_vars
            only_marg <- gsub("\\..*","", selected_interactions)
            sel <- which(only_marg == marg_var)
            selected_interactions <- inter_union$interaction[, selected_interactions[sel]]
            valid_df <- data.frame(valid_z[, marg_var], selected_interactions)
            colnames(valid_df) <- c(marg_var, res$interaction_vars[sel])
            
            upsilon_name <- paste0("ups", marg_var)
            fit_vals <- stats::predict(object = res$upsilon_mods[[upsilon_name]], newdata = valid_df)
            valid_upsilon <- fit_vals %>%
              as.data.frame()
            colnames(valid_upsilon) <- upsilon_name
            zs_ints_valid <- cbind(zs_ints_valid, valid_upsilon)
          }
          selected_upsilon <- colnames(res$upsilon)
          zs_valid <- zs_ints_valid[, selected_upsilon] %>%
            as.data.frame()
          colnames(zs_valid) <- selected_upsilon
        } else { ## no interactions
          cat("NO INTERACTIONS...USING MARGINALS \n")
          zs_train <- cbind(train_z[, res$marginal_vars]) %>%
            as.data.frame()
          colnames(zs_train) <- c(res$marginal_vars)
          zs_valid <- cbind(valid_z[, res$marginal_vars]) %>%
            as.data.frame()
          colnames(zs_valid) <- c(res$marginal_vars)
        }
        ## make model
        if (y_factor) {
          slide_mod <- stats::glm(use_y_train ~ ., data = zs_train, family = "binomial")
        } else {
          slide_mod <- stats::lm(use_y_train ~ ., data = zs_train)
        }
        
        ## do prediction
        pred_vals <- predict(slide_mod, zs_valid, type = "response")
        
        ## add SLIDE results to replicate report
        cat("  SLIDE model size: ", ncol(zs_train), "\n")
      }
      
      if (eval_type == "auc") { ## if using area under roc curve to evaluate model fit
        pred_vals <- as.data.frame(pred_vals)
        if (ncol(pred_vals) > 1) {
          pred_vals <- pred_vals[, 1]
        }
        fold_res <- cbind(spec, pred_vals, valid_y_labs, valid_ind)
        colnames(fold_res) <- c("spec", "pred_vals", "true_vals", "index")
      } else { ## if using correlation to evaluate model fit
        fold_res <- cbind(spec, pred_vals, valid_y, valid_ind)
        colnames(fold_res) <- c("spec", "pred_vals", "true_vals", "index")
      }
      results <- rbind(results, fold_res)
      tictoc::toc()
      cat("\n")
    }
    tictoc::toc()
    cat("\n")
  }
  
  ## set results data frame column names
  results <- as.data.frame(results)
  saveRDS(results, file = paste0(out_path, "replicate", rep, "/predicted_values.rds"))
  
  final_results <- NULL
  if (eval_type == "auc") {
    for (i in 1:length(specs)) {
      spec_res <- results %>%
        dplyr::filter(spec == specs[i])
      predicted <- spec_res$pred_vals
      true <- spec_res$true_vals
      spec_roc <- ROCR::prediction(predicted, true)
      spec_auc <- ROCR::performance(spec_roc, method = "auc")
      spec_auc <- spec_auc@y.values[[1]]
      if (spec_auc < 0.5) { ## if classifier auc is < 0.5, reverse it to be > 0.5
        spec_auc <- 1 - spec_auc
      }
      spec_res <- c("spec" = specs[i],
                    "auc" = as.numeric(spec_auc))
      final_results <- rbind(final_results, spec_res)
    }
  } else {
    for (i in 1:length(specs)) {
      spec_res <- results %>%
        dplyr::filter(spec == specs[i])
      predicted <- as.numeric(spec_res$pred_vals)
      true <- as.numeric(spec_res$true_vals)
      spec_corr <- cor(predicted, true, method = "spearman")
      spec_res <- c("spec" = specs[i],
                    "corr" = as.numeric(spec_corr))
      final_results <- rbind(final_results, spec_res)
    }
  }
  tictoc::toc()
  cat("\n")
  
  final_results <- as.data.frame(final_results)
  final_results[, 2] <- as.numeric(final_results[, 2])
  saveRDS(final_results, file = paste0(new_dir, "model_evaluations.rds"))
  
  return (final_results)
}
