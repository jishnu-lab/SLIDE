#' Essential Regression Pipeline - Step 5
#'
#' Run Essential Regression Pipeline with K-Fold Cross-Validation
#'     Step 5 - Run final replicates of K-Fold CV using optimal \eqn{\delta} and optimal \eqn{\lambda}.
#'              Create a boxplot and run \code{plainER()} one last time for final results.
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param yaml_path the path to a .yaml file containing all necessary parameters/arguments
#' for Essential Regression
#' @return nothing is returned, saves boxplot of cross-validation results and final ER output
#' @export

pipelineER3 <- function(yaml_path) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1, check.names = F)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1, check.names = F)) ## not standardized

  ## remove any zero sd cols
  x = x[, which(apply(x, 2, sd) != 0)]
  x_std <- scale(x, T, T)

  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)  

  if (er_input$y_factor) {
    y <- toCont(y, er_input$y_levels)
    saveRDS(y, file = paste0(er_input$out_path, "pipeline3_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }

  # check with benchmark methods we're doing
  run_lasso = F
  if (!is.null(er_input$lasso) & er_input$lasso) {
    run_lasso = T
  }

  if (er_input$k <= 0) {
    er_input$k <- nrow(x) #LOOCV
  }

  lambda_rep = data.frame()
  for (j in 1:er_input$nreps) {
    temp <- NULL
    while (is.null(temp)) {
      temp <- essregCV(k = er_input$k,
                       x = x,
                       y = y,
                       std_y = er_input$std_y,
                       std_cv = er_input$std_cv,
                       delta = er_input$delta,
                       permute = er_input$permute,
                       eval_type = er_input$eval_type,
                       y_levels = er_input$y_levels,
                       lambda = er_input$lambda,
                       out_path = er_input$out_path,
                       rep_cv = er_input$rep_cv,
                       alpha_level = er_input$alpha_level,
                       thresh_fdr = er_input$thresh_fdr,
                       rep = j,
                       benchmark = er_input$benchmark,
                       run_lasso = run_lasso)
    }
    lambda_rep = rbind(lambda_rep, temp)
  }
  saveRDS(lambda_rep, paste0(er_input$out_path, "pipeline_step5.rds"))

  ## create boxplot of replicate correlations ##################################
  if (er_input$eval_type == "auc") {
    methods <- c("plainER", "plainER_y", "lasso", "lasso_y", "pclr", "pclr_y", "plsda", "plsda_y")
  } else {
    methods <- c("plainER", "plainER_y", "lasso", "lasso_y", "pcr", "pcr_y", "plsr", "plsr_y")
  }
  bp_df <- lambda_rep %>%
    dplyr::mutate(perm = sub(".*_", "", method)) %>%
    dplyr::mutate(perm = ifelse(perm == method, "no_perm", paste0(perm, "_perm"))) %>%
    dplyr::mutate(method_perm = sub("*_.", "", method)) %>%
    dplyr::mutate(method = factor(method, levels = methods),
                  perm = as.factor(perm)) %>%
    dplyr::mutate(alpha = ifelse(perm == "no_perm", 1, 0.9))

  pdf_file <- paste0(er_input$out_path, "/opt_delta_lambda_boxplot.pdf")
  dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)

  if (er_input$eval_type == "corr") {
    lambda_boxplot <- ggplot2::ggplot(data = bp_df,
                                      ggplot2::aes(x = method,
                                                   y = corr,
                                                   fill = method_perm,
                                                   alpha = alpha)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(fill = "Method") +
      ggplot2::scale_alpha(guide = 'none')
  } else {
    lambda_boxplot <- ggplot2::ggplot(data = bp_df,
                                      ggplot2::aes(x = method,
                                                   y = auc,
                                                   fill = method_perm,
                                                   alpha = alpha)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(fill = "Method") +
      ggplot2::scale_alpha(guide = 'none')
  }

  ggplot2::ggsave(pdf_file, lambda_boxplot, width = 20, height = 15, units = "in")

  ## Final plainER Run #########################################################
  er_output <- plainER(x = x,
                       y = y,
                       x_std = x_std,
                       std_y = er_input$std_y,
                       sigma = NULL,
                       delta = er_input$delta,
                       lambda = er_input$lambda,
                       rep_cv = er_input$rep_cv,
                       alpha_level = er_input$alpha_level,
                       thresh_fdr = er_input$thresh_fdr,
                       out_path = er_input$out_path)

  if (is.null(er_output)) {
    cat("plainER failed --- infeasible linear program \n")
    return ()
  }

  saveRDS(er_output, paste0(er_input$out_path, "final_delta_", er_input$delta, "_lambda_", er_input$lambda, ".rds"))
}

