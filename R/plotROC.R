plotROC <- function(er_input) {
  ## housekeeping ##############################################################
  pred_values <- list("plainER" = NULL,
                      "plainER_y" = NULL,
                      "HER" = NULL,
                      "HER_rand" = NULL,
                      "lasso" = NULL,
                      "lasso_y" = NULL,
                      "phate" = NULL,
                      "plsda" = NULL,
                      "pclr" = NULL)
  methods <- c("plainER", "plainER_y", "HER", "HER_rand", 
               "lasso", "lasso_y", "phate", "plsda", "pclr")
  
  ## create output directory
  new_dir <- paste0(out_path, "ROC_curves/")
  dir.create(file.path(new_dir), showWarnings = F, recursive = T)
  
  ## get predicted values from replicates of benchCV ###########################
  for (i in 1:er_input$nreps) { ## loop through replicate results
    rep_path <- paste0(er_input$out_path, "replicate", i, "/")
    method_res <- readRDS(paste0(rep_path, "predicted_values.rds"))
    for (j in 1:length(methods)) { ## loop through methods
      method_preds <- pred_values[[methods[j]]]
      method_filt <- method_res %>%
        dplyr::filter(method == methods[j]) %>%
        dplyr::mutate(pred_vals = as.numeric(pred_vals),
                      true_vals = as.numeric(true_vals),
                      index = as.factor(index)) %>%
        dplyr::select(-method)
      method_preds <- rbind(method_preds, method_filt)
      pred_values[[methods[j]]] <- method_preds
    }
  }
  
  ## make ROC curves ###########################################################
  roc_df <- NULL
  for (i in 1:length(methods)) { ## loop through methods
    method_preds <- pred_values[[methods[i]]]
    method_preds <- method_preds %>%
      dplyr::group_by(index) %>%
      dplyr::summarise(pred_vals = stats::median(pred_vals),
                       true_vals = stats::median(true_vals))

    predicted <- method_preds$pred_vals
    true <- method_preds$true_vals
    method_roc <- ROCR::prediction(predicted, true)
    method_auc <- ROCR::performance(method_roc)
    method_auc <- method_auc@y.values
    if (method_auc < 0.5) { ## if classifier auc is < 0.5, reverse it to be > 0.5
      method_auc <- 1 - method_auc
    }
    
    method_roc <- ggplot2::ggroc(method_roc, colour = 'steelblue', size = 2) +
      ggplot2::ggtitle(paste0("Median ROC Curve for ", method[i])) +
      ggplot2::theme_minimal() +
      ggplot2::annotate("text", x = 0.9, y = 0.9, label = paste0("AUC = ", method_auc))
    
    ggplot2::ggsave(paste0(new_dir, method[i], "_ROC.pdf"), bench_plot, width = 20, height = 15, units = "in")
  }
  
  
  
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  
  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)
  
  ## create violin plot of replicate correlations ##############################
  all_results <- list.files(path = er_input$out_path,
                            pattern = "model_evaluations.rds",
                            recursive = TRUE)
  
  bench_reps <- NULL
  for (i in 1:length(all_results)) {
    read_path <- paste0(er_input$out_path, all_results[[i]])
    rep_res <- readRDS(read_path)
    bench_reps <- rbind(bench_reps, rep_res)
  }
  
  ## set method levels
  if (er_input$eval_type == "corr") {
    method_order <- c("HER", "plainER", "lasso", "phate", "plsr", "pcr", "plainER_y", "lasso_y")
  } else {
    method_order <- c("HER", "plainER", "lasso", "phate", "plsda", "pclr", "plainER_y", "lasso_y")
  }
  
  viol_df <- bench_reps %>%
    dplyr::mutate(method = factor(method, levels = method_order)) %>%
    dplyr::mutate(perm = sub(".*_", "", method)) %>%
    dplyr::mutate(perm = ifelse(perm == method, "no_perm", paste0(perm, "_perm"))) %>%
    dplyr::mutate(method_perm = sub("*_.", "", method)) %>%
    dplyr::mutate(method = as.factor(method),
                  perm = as.factor(perm)) %>%
    dplyr::mutate(alpha = ifelse(perm == "no_perm", 1, 0.95))
  
  pdf_file <- paste0(er_input$out_path, "/benchmarks_violinplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
  
  if (er_input$eval_type == "corr") {
    bench_plot <- ggplot2::ggplot(data = viol_df,
                                  ggplot2::aes(x = method,
                                               y = corr,
                                               fill = method_perm,
                                               alpha = alpha)) +
      ggplot2::geom_violin() +
      ggplot2::labs(fill = "Method") +
      ggplot2::ylim(-1, 1) +
      ggplot2::scale_alpha(guide = 'none') +
      ggplot2::geom_vline(xintercept = 6.5, linetype = "dotted", size = 1.5)
  } else {
    bench_plot <- ggplot2::ggplot(data = viol_df,
                                  ggplot2::aes(x = method,
                                               y = auc,
                                               fill = method_perm,
                                               alpha = alpha)) +
      ggplot2::geom_violin() +
      ggplot2::labs(fill = "Method") +
      ggplot2::scale_alpha(guide = 'none') +
      ggplot2::ylim(0, 1) +
      ggplot2::geom_vline(xintercept = 6.5, linetype = "dotted", size = 1.5)
  }
  
  ggplot2::ggsave(pdf_file, bench_plot, width = 10, height = 7, units = "in")
}