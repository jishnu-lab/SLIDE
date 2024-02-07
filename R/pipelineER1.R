#' Essential Regression Pipeline - Steps 1, 2
#'
#' Run Essential Regression Pipeline with K-Fold Cross-Validation
#'     Step 1 - Coarse Grid Search for \eqn{\delta}
#'     Step 2 - K-Fold Cross-Validation to find magnitude of \eqn{\delta}
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param yaml_path the path to a .yaml file containing all necessary parameters/arguments
#' for Essential Regression
#' @param steps an integer or string indicating which steps of the pipeline to perform: "1" or "all"
#' @return nothing is returned, saves boxplot of cross-validation results for user to use
#' in selecting optimal \eqn{\delta} if \code{steps} is "all"
#' @export

pipelineER1 <- function(yaml_path, steps = "all") {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1, check.names = F)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1, check.names = F)) ## not standardized
  
  x_std <- scale(x, T, T)

  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)

  if (er_input$y_factor) {
    y <- toCont(y, er_input$y_levels)
    saveRDS(y, file = paste0(er_input$out_path, "pipeline1_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
  }

  deltas <- list(seq(0.0001, 0.001, 0.0001),
                 seq(0.001, 0.01, 0.001),
                 seq(0.01, 0.1, 0.01),
                 seq(0.1, 1, 0.1))

  # check with benchmark methods we're doing
  run_lasso = F
  if (!is.null(er_input$lasso) & er_input$lasso) {
    run_lasso = T
  }

  ## Step 1: Coarse Delta Search #############################################
  if (file.exists(paste0(er_input$out_path, "pipeline_step1.rds"))) {
    coarse_res <- readRDS(paste0(er_input$out_path, "pipeline_step1.rds"))
  } else {
    coarse_res = list()
    for (i in 1:length(deltas)) {
      cat(i, "\n")
      result <- plainER(y = y,
                        x = x, # x here is NOT z_scored x
                        x_std = x_std,
                        std_y = er_input$std_y,
                        sigma = NULL,
                        delta = deltas[[i]],
                        lambda = 0.5,
                        rep_cv = er_input$rep_cv,
                        alpha_level = er_input$alpha_level,
                        thresh_fdr = er_input$thresh_fdr,
                        out_path = er_input$out_path)
      if (is.null(result)) {
        cat("plainER failed --- infeasible linear program \n")
        return ()
      }
      coarse_res[[length(coarse_res) + 1]] = result
    }
    saveRDS(coarse_res, file = paste0(er_input$out_path, "pipeline_step1.rds"))
  }

  ## Step 2: K-Fold Cross-Validation To Find Delta Magnitude #################
  if (steps == "all") {
    corr_bp_data <- NULL
    for (i in 1:length(coarse_res)) {
      mag_delta <- coarse_res[[i]]$opt_delta
      magnitude <- deltas[[i]][1]
      cat("DELTA = ", mag_delta, " . . . \n")
      if (file.exists(paste0(er_input$out_path, "essregCV_delta_", mag_delta, ".rds"))) {
        delta_rep <- readRDS(paste0(er_input$out_path, "essregCV_delta_", mag_delta, ".rds"))
      } else {
        foreach::foreach (j = 1:er_input$nreps, .combine = rbind) %dopar% {
          if (file.exists(file = paste0(er_input$out_path, "delta_", mag_delta, "/replicate", j, "/output_table.rds"))) {
            readRDS(paste0(er_input$out_path, "delta", mag_delta, "/replicate", j, "/output_table.rds"))
          } else {
            result <- NULL
            while (is.null(result)) {
              result <- essregCV(k = er_input$k,
                                 x = x,
                                 y = y,
                                 std_cv = er_input$std_cv,
                                 std_y = er_input$std_y,
                                 delta = mag_delta,
                                 eval_type = er_input$eval_type,
                                 permute = er_input$permute,
                                 y_levels = er_input$y_levels,
                                 lambda = 0.5,
                                 rep_cv = er_input$rep_cv,
                                 alpha_level = er_input$alpha_level,
                                 thresh_fdr = er_input$thresh_fdr,
                                 out_path = paste0(er_input$out_path, "step2_delta_", mag_delta, "/"),
                                 rep = j,
                                 benchmark = er_input$benchmark,
                                 run_lasso = run_lasso)
            }
            result
          }
        } -> delta_rep
        #saveRDS(delta_rep, file = paste0(er_input$out_path, "essregCV_delta_", mag_delta, ".rds"))
      }
      ## make CV plot
      if (er_input$eval_type == "auc") {
        methods <- c("plainER", "plainER_y", "lasso", "lasso_y", "pclr", "pclr_y", "plsda", "plsda_y")
      } else {
        methods <- c("plainER", "plainER_y", "lasso", "lasso_y", "pcr", "pcr_y", "plsr", "plsr_y")
      }
      final_res <- delta_rep %>%
        dplyr::mutate(perm = sub(".*_", "", method)) %>%
        dplyr::mutate(perm = ifelse(perm == method, "no_perm", paste0(perm, "_perm"))) %>%
        dplyr::mutate(method_perm = sub("*_.", "", method)) %>%
        dplyr::mutate(method = factor(method, levels = methods),
                      perm = as.factor(perm)) %>%
        dplyr::mutate(alpha = ifelse(perm == "no_perm", 1, 0.9))
      pdf_file <- paste0(er_input$out_path, "delta_", mag_delta, "_boxplot.pdf")
      dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
      if (er_input$eval_type == "corr") {
        delta_boxplot <- ggplot2::ggplot(data = final_res,
                                         ggplot2::aes(x = method,
                                                      y = corr,
                                                      fill = method_perm,
                                                      alpha = alpha)) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(fill = "Method") +
          ggplot2::scale_alpha(guide = 'none')
      } else {
        delta_boxplot <- ggplot2::ggplot(data = final_res,
                                         ggplot2::aes(x = method,
                                                      y = auc,
                                                      fill = method_perm,
                                                      alpha = alpha)) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(fill = "Method") +
          ggplot2::scale_alpha(guide = 'none')
      }
      ggplot2::ggsave(pdf_file, delta_boxplot, width = 20, height = 15, units = "in")
      corr_bp_data[[length(corr_bp_data) + 1]] <- list("delta" = mag_delta,
                                                       "result" = delta_rep)
    }
    saveRDS(corr_bp_data, file = paste0(er_input$out_path, "pipeline_step2.rds"))

    ## create boxplot of replicate performance #################################
    final_res <- NULL
    for (i in 1:length(corr_bp_data)) {
      bp_data <- corr_bp_data[[i]]
      bp_delta <- bp_data$delta
      bp_df <- bp_data$result %>%
        dplyr::filter(method == "plainER" | method == "plainER_y") %>%
        dplyr::mutate(delta = bp_delta)
      final_res <- rbind(final_res, bp_df)
    }
    final_res <- final_res %>%
      dplyr::mutate(delta = as.factor(delta),
                    method = as.factor(method))
    pdf_file <- paste0(er_input$out_path, "/delta_selection_boxplot.pdf")
    dir.create(file.path(dirname(pdf_file)), showWarnings = F, recursive = T)
    if (er_input$eval_type == "corr") {
      delta_boxplot <- ggplot2::ggplot(data = final_res,
                                       ggplot2::aes(x = delta,
                                                    y = corr,
                                                    fill = method)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(fill = "Method")
    } else {
      delta_boxplot <- ggplot2::ggplot(data = final_res,
                                       ggplot2::aes(x = delta,
                                                    y = auc,
                                                    fill = method)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(fill = "Method")
    }

    ggplot2::ggsave(pdf_file, delta_boxplot, width = 20, height = 15, units = "in")
  }
}
