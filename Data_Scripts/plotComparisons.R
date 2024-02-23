plotComparisons <- function(yaml_path) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized
  
  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)
  
  if (er_input$eval_type == "auc") {
    if (length(er_input$y_order) > 2) {
      stop("y must be binary - please refactor")
    }
    y <- toCont(y, er_input$y_order)
    saveRDS(y, file = paste0(er_input$out_path, "pipeline3_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }
  
  ## do replicates to make violin plots
  foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %do% {
    cat("replicate", j, "\n")
    temp <- NULL
    while (is.null(temp)) {
      temp <- compareER(k = er_input$k,
                        x = x,
                        y = y,
                        delta = er_input$delta,
                        eval_type = er_input$eval_type,
                        lambda = er_input$lambda,
                        out_path = er_input$out_path,
                        rep_cv = er_input$rep_cv,
                        alpha_level = er_input$alpha_level,
                        thresh_fdr = er_input$thresh_fdr,
                        parallel = er_input$parallel,
                        y_order = er_input$y_order,
                        spec = er_input$spec,
                        niter = er_input$niter,
                        f_size = er_input$f_size,
                        fdr = er_input$fdr,
                        ncore = er_input$ncore,
                        rep = j)
    }
    temp
  } -> comp_reps
  saveRDS(bench_reps, paste0(er_input$out_path, "benchmarks.RDS"))
  
  ## create violin plot of replicate correlations ##################################
  viol_df <- comp_reps %>%
    dplyr::mutate(method = as.factor(method))
  
  pdf_file <- paste0(er_input$out_path, "/ER_comparisons_violinplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
  
  if (er_input$eval_type == "corr") {
    comp_plot <- ggplot2::ggplot(data = viol_df,
                                 ggplot2::aes(x = method,
                                              y = corr,
                                              fill = method)) +
      ggplot2::geom_violin() +
      ggplot2::labs(fill = "Method")
  } else {
    comp_plot <- ggplot2::ggplot(data = viol_df,
                                 ggplot2::aes(x = method,
                                              y = auc,
                                              fill = method)) +
      ggplot2::geom_violin() +
      ggplot2::labs(fill = "Method")
  }
  
  ggplot2::ggsave(pdf_file, comp_plot, width = 20, height = 15, units = "in")
}

