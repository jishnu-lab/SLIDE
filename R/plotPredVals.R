plotPredVals <- function(yaml_path) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  
  ## get names of .rds files
  all_results <- list.files(path = er_input$out_path,
                            pattern = "predicted_values.rds",
                            recursive = TRUE)
  
  ## read in results from each replicate
  replicates <- NULL
  for (i in 1:length(all_results)) {
    read_path <- paste0(er_input$out_path, all_results[[i]])
    rep_res <- readRDS(read_path)
    replicates <- rbind(replicates, rep_res)
  }
  
  ## set method levels
  if (er_input$eval_type != "corr") {
    stop("This function is only for continuous y. \n")
  } else {
    method_order <- c("SLIDE", "plainER", "lasso", "phate", "plsr", "pcr", "plainER_y", "lasso_y")
  }
  
  ## loop through different methods
  for (i in 1:length(method_order)) {
    ## get one method
    meth <- method_order[i]
    ## create plot file name
    plot_name <- paste0(er_input$out_path, meth, "_predicted_values_plot.pdf")
    
    ## filter predicted values to just one method
    pred_vals <- replicates %>%
      dplyr::filter(method == meth) %>%
      dplyr::group_by(index) %>%
      dplyr::mutate(pred_vals = as.numeric(pred_vals)) %>%
      dplyr::mutate(true_vals = as.numeric(true_vals)) %>%
      dplyr::mutate(index = factor(index, levels = seq(1, unique(index))))
    
    ## create violin plot of predicted values
    pred_plot <- ggplot2::ggplot(data = pred_vals,
                                 ggplot2::aes(x = index,
                                              y = pred_vals,
                                              fill = index)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_point(ggplot2::aes(x = index,
                                       y = true_vals,
                                       fill = "#FF0000"),
                          pch = 21,
                          size = 2,
                          colour = "black") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "observation",
                    y = "predicted values",
                    title = paste0("Predicted Values Violin Plot - ", meth))
    
    ## save plot
    ggplot2::ggsave(plot_name, pred_plot, width = 10, height = 8, units = "in")
  }
}