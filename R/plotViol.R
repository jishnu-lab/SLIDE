#' Create a violin plot of results.
#'
#' Reads in results from running the benchmarking cross-validation regime and creates a
#' violin plot of the replicate results.
#' 
#' @param yaml_path a string path to the .yaml file that configured the cross-validation
#' @return none
#' @export

plotViol <- function(yaml_path) {
  require(dplyr)
  require(ggplot2)
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  
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
    method_order <- c("SLIDE", "plainER", "lasso", "phate", "plsr", "pcr", "plainER_y", "lasso_y")
  } else {
    method_order <- c("SLIDE", "plainER", "lasso", "phate", "plsda", "pclr", "plainER_y", "lasso_y")
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
    #### ORIGINAL VIOLINPLOT
    # bench_plot <- ggplot2::ggplot(data = viol_df,
    #                               ggplot2::aes(x = method,
    #                                            y = corr,
    #                                            fill = method_perm,
    #                                            alpha = alpha)) +
    #   ggplot2::geom_violin(trim = TRUE) +
    #   ggplot2::geom_point(stat = "summary", fun = "mean", alpha = 1) +
    #   ggplot2::geom_errorbar(stat = "summary", 
    #                          fun.data = "mean_se", 
    #                          width = 0.2, 
    #                          alpha = 1,
    #                          fun.args = list(mult = 1.96)) +
    #   ggplot2::labs(fill = "Method") +
    #   ggplot2::ylim(-1, 1) +
    #   ggplot2::scale_alpha(guide = 'none') +
    #   ggplot2::geom_vline(xintercept = 6.5, linetype = "dotted", size = 1.5)
    #### BOXPLOT
    # bench_plot <- ggplot2::ggplot(data = viol_df,
    #                               ggplot2::aes(x = method,
    #                                            y = corr,
    #                                            fill = method_perm,
    #                                            alpha = alpha)) +
    #   ggplot2::geom_boxplot() +
    #   ggplot2::labs(fill = "Method") +
    #   ggplot2::ylim(-1, 1) +
    #   ggplot2::scale_alpha(guide = 'none') +
    #   ggplot2::geom_vline(xintercept = 6.5, linetype = "dotted", size = 1.5)
    #### MODIFIED VIOLINPLOT
    # bench_plot <- ggplot2::ggplot(data = viol_df,
    #                               ggplot2::aes(x = method,
    #                                            y = corr,
    #                                            fill = method_perm,
    #                                            alpha = alpha)) +
    #   ggplot2::geom_violin(trim = TRUE, scale = 'width') +
    #   ggplot2::geom_boxplot(width = 0.1) +
    #   ggplot2::labs(fill = "Method") +
    #   ggplot2::ylim(-1, 1) +
    #   ggplot2::scale_alpha(guide = 'none') +
    #   ggplot2::geom_vline(xintercept = 6.5, linetype = "dotted", size = 1.5)
    #### MODIFIED VIOLINPLOT #2
    bench_plot <- ggplot2::ggplot(data = viol_df,
                                  ggplot2::aes(x = method,
                                               y = valid_corr,
                                               fill = method_perm,
                                               alpha = alpha)) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::labs(fill = "Method") +
      ggplot2::ylim(-0.1, 1) +
      ggplot2::scale_alpha(guide = 'none') +
      ggplot2::geom_vline(xintercept = 6.5, linetype = "dotted", size = 1.5)+
      ggplot2::ylab("The spearman correlation")+
      ggplot2::theme_light()+ 
      ggplot2::theme_bw()+
      ggplot2:: theme(legend.position = "none")   +
      ggplot2::theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    ggplot2::ggsave("/ix/djishnu/Javad/Poholek_scale_v_noscale/w_scale_220922/plots/benchplot.pdf",width = 10,height = 10)
    
    
  } else {
    bench_plot <- ggplot2::ggplot(data = viol_df,
                                  ggplot2::aes(x = method,
                                               y = auc,
                                               fill = method_perm,
                                               alpha = alpha)) +
      ggplot2::geom_point(stat = "summary", fun = "mean", alpha = 1) +
      ggplot2::geom_violin(trim = TRUE) +
      ggplot2::geom_errorbar(stat = "summary", 
                             fun.data = "mean_se", 
                             width = 0.2, 
                             alpha = 1,
                             fun.args = list(mult = 1.96)) +
      ggplot2::labs(fill = "Method") +
      ggplot2::scale_alpha(guide = 'none') +
      ggplot2::ylim(-0.25, 1) +
      ggplot2::geom_vline(xintercept = 6.5, linetype = "dotted", size = 1.5)+
      ggplot2::theme_light()
  }
  
  ggplot2::ggsave(pdf_file, bench_plot, width = 10, height = 10, units = "in")
}