#' Essential Regression Pipeline - Steps 3, 4
#'
#' Run Essential Regression Pipeline with K-Fold Cross-Validation
#'     Step 3 - Fine Grid Search for \eqn{\delta}
#'     Step 4 - K-Fold Cross-Validation to find \eqn{\lambda}
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param yaml_path the path to a .yaml file containing all necessary parameters/arguments
#' for Essential Regression
#' @param steps an integer or string indicating which steps of the pipeline to perform: "3" or "all"
#' @return nothing is returned, saves boxplot of cross-validation results for user to use
#' in selecting optimal \eqn{\lambda}
#' @export

pipelineER2 <- function(yaml_path, steps = "all") {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized

  x_std <- scale(x, T, T)

  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)

  if (er_input$y_factor) {
    y <- toCont(y, er_input$y_levels)
    saveRDS(y, file = paste0(er_input$out_path, "pipeline2_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }

  # check with benchmark methods we're doing
  run_lasso = F
  if (!is.null(er_input$lasso) & er_input$lasso) {
    run_lasso = T
  }

  ## Step 3: Fine Delta Search ###############################################
  if (steps == "all"){
    if (file.exists(paste0(er_input$out_path, "pipeline_step3.rds"))) {
      fine_delta_er <- readRDS(file = paste0(er_input$out_path, "pipeline_step3.rds"))
    } else {
      if (length(er_input$delta) == 1) {
        d_lbd <- er_input$delta - er_input$delta / 2
        d_ubd <- er_input$delta + er_input$delta / 2
        delta_grid <- seq(d_lbd, d_ubd, er_input$delta / 100)
      } else {
        delta_grid <- er_input$delta
      }

      fine_delta_er <- plainER(y = y,
                               x = x,
                               x_std = x_std,
                               std_y = er_input$std_y,
                               sigma = cor(x),
                               delta = delta_grid,
                               lambda = 0.5,
                               rep_cv = er_input$rep_cv,
                               alpha_level = er_input$alpha_level,
                               thresh_fdr = er_input$thresh_fdr,
                               out_path = er_input$out_path)
      if (is.null(fine_delta_er)) {
        cat("plainER failed --- infeasible linear program \n")
        return ()
      }
      saveRDS(fine_delta_er, file = paste0(er_input$out_path, "pipeline_step3.rds"))
    }

    best_delta <- fine_delta_er$opt_delta
  }
  ## Step 4: Lambda Search ###################################################
  if (steps == "4" | steps == "all") {
    corr_bp_data <- NULL
    if (length(er_input$delta) == 1) {
      best_delta <- er_input$delta
    } else {
      cat("Delta is given as multiple values...")
    }
    for (i in 1:length(er_input$lambda)) {

      lambda <- er_input$lambda[[i]]
      cat("LAMBDA = ", lambda, " . . . \n")
      out_path <- paste0(er_input$out_path, "lambda_", lambda, "/")
      if (file.exists(paste0(er_input$out_path, "essregCV_lambda_", lambda, ".rds"))) {
        lambda_rep <- readRDS(paste0(er_input$out_path, "essregCV_lambda_", lambda, ".rds"))
      } else {
        lambda_rep = data.frame()
        for (j in 1:er_input$nreps) {
          if (file.exists(file = paste0(out_path, "replicate", j, "/output_table.rds"))) {
            readRDS(paste0(out_path, "replicate", j, "/output_table.rds"))
          } else {
            result <- NULL
            while (is.null(result)) {
              result <- essregCV(k = er_input$k,
                                 x = x,
                                 y = y,
                                 std_y = er_input$std_y,
                                 std_cv = er_input$std_cv,
                                 delta = best_delta,
                                 permute = er_input$permute,
                                 eval_type = er_input$eval_type,
                                 y_levels = er_input$y_levels,
                                 lambda = lambda,
                                 out_path = paste0(er_input$out_path, "lambda_", lambda, "/"),
                                 rep_cv = er_input$rep_cv,
                                 alpha_level = er_input$alpha_level,
                                 thresh_fdr = er_input$thresh_fdr,
                                 rep = j,
                                 benchmark = er_input$benchmark,
                                 run_lasso = run_lasso)
            }
            lambda_rep = rbind(lambda_rep, result)
          }
        }
        # foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
        #   if (file.exists(file = paste0(out_path, "replicate", j, "/output_table.rds"))) {
        #     readRDS(paste0(out_path, "replicate", j, "/output_table.rds"))
        #   } else {
        #     result <- NULL
        #     while (is.null(result)) {
        #       result <- essregCV(k = er_input$k,
        #                          x = x,
        #                          y = y,
        #                          std_y = er_input$std_y,
        #                          std_cv = er_input$std_cv,
        #                          delta = best_delta,
        #                          permute = er_input$permute,
        #                          eval_type = er_input$eval_type,
        #                          y_levels = er_input$y_levels,
        #                          lambda = lambda,
        #                          out_path = paste0(er_input$out_path, "lambda_", lambda, "/"),
        #                          rep_cv = er_input$rep_cv,
        #                          alpha_level = er_input$alpha_level,
        #                          thresh_fdr = er_input$thresh_fdr,
        #                          rep = j,
        #                          benchmark = er_input$benchmark)
        #     }
        #     result
        #   }
        # } -> lambda_rep

        #saveRDS(lambda_rep, file = paste0(er_input$out_path, "essregCV_lambda_", lambda, ".rds"))
      }

      ## make CV plot
      if (er_input$eval_type == "auc") {
        methods <- c("plainER", "plainER_y", "lasso", "lasso_y", "pclr", "pclr_y", "plsda", "plsda_y")
      } else {
        methods <- c("plainER", "plainER_y", "lasso", "lasso_y", "pcr", "pcr_y", "plsr", "plsr_y")
      }
      final_res <- lambda_rep %>%
        dplyr::mutate(perm = sub(".*_", "", method)) %>%
        dplyr::mutate(perm = ifelse(perm == method, "no_perm", paste0(perm, "_perm"))) %>%
        dplyr::mutate(method_perm = sub("*_.", "", method)) %>%
        dplyr::mutate(method = factor(method, levels = methods),
                      perm = as.factor(perm)) %>%
        dplyr::mutate(alpha = ifelse(perm == "no_perm", 1, 0.9))

      pdf_file <- paste0(er_input$out_path, "lambda_", lambda, "_boxplot.pdf")
      dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)

      if (er_input$eval_type == "corr") {
        lambda_boxplot <- ggplot2::ggplot(data = final_res,
                                          ggplot2::aes(x = method,
                                                       y = corr,
                                                       fill = method_perm,
                                                       alpha = alpha)) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(fill = "Method") +
          ggplot2::scale_alpha(guide = 'none')
      } else {
        lambda_boxplot <- ggplot2::ggplot(data = final_res,
                                          ggplot2::aes(x = method,
                                                       y = auc,
                                                       fill = method_perm,
                                                       alpha = alpha)) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(fill = "Method") +
          ggplot2::scale_alpha(guide = 'none')
      }

      ggplot2::ggsave(pdf_file, lambda_boxplot, width = 20, height = 15, units = "in")
      corr_bp_data[[length(corr_bp_data) + 1]] <- list("lambda" = lambda,
                                                       "result" = lambda_rep)
    }
    saveRDS(corr_bp_data, file = paste0(er_input$out_path, "pipeline_step4.rds"))

    ## create boxplot of replicate correlations ################################
    final_res <- NULL
    for (i in 1:length(corr_bp_data)) {
      bp_data <- corr_bp_data[[i]]
      bp_lambda <- bp_data$lambda
      bp_df <- bp_data$result %>%
        dplyr::filter(method == "plainER" | method == "plainER_y") %>%
        dplyr::mutate(lambda = bp_lambda)
      final_res <- rbind(final_res, bp_df)
    }
    final_res <- final_res %>%
      dplyr::mutate(lambda = as.factor(lambda),
                    method = as.factor(method))
    pdf_file <- paste0(er_input$out_path, "lambda_selection_boxplot.pdf")
    dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)

    if (er_input$eval_type == "corr") {
      lambda_boxplot <- ggplot2::ggplot(data = final_res,
                                        ggplot2::aes(x = lambda,
                                                     y = corr,
                                                     fill = method)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(fill = "Method")
    } else {
      lambda_boxplot <- ggplot2::ggplot(data = final_res,
                                        ggplot2::aes(x = lambda,
                                                     y = auc,
                                                     fill = method)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(fill = "Method")
    }

    ggplot2::ggsave(pdf_file, lambda_boxplot, width = 20, height = 15, units = "in")
  }
}
