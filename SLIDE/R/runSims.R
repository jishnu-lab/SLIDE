#' Run simulations.
#'
#' Given different parameter sets, run simulations and save the results for plotting.
#' Of \code{N}, \code{P}, \code{K}, \code{sK}, and \code{sI}, one should be a vector of values to
#' loop over and the rest should be single values. This way the varying values of a single parameter
#' can be plotted on the x-axis against the model MSEs on the y-axis for the simulation line plots.
#'
#' @param N a single numeric value or a vector; sample size
#' @param P a single numeric value or a vector; number of features
#' @param K a single numeric value or a vector; number of clusters
#' @param sK a single numeric value or a vector; number of significant clusters used to
#' generate \eqn{Y}
#' @param sI a single numeric value or a vector; number of significant interaction terms used
#' to generate \eqn{Y}
#' @param spec a numeric constant; the proportion of iterations that a variable must be selected in to be considered significant
#' @param niter an integer; the number of times to run knockoffs
#' @param fdr a numeric constant between 0 and 1; the target false discovery rate for knockoffs
#' @param f_size an integer; the target size for the subset (number of columns in each subset)
#' @param ncore an integer; the number of cores to use for parallel computing
#' @param out_path a string path; where to save output
#' @return none
#' @export

runSims <- function(N, P, K, sK, sI, spec, niter, fdr, f_size, ncore, out_path,seed) {
  ## create output directory 
  dir.create(file.path(out_path), showWarnings = F, recursive = T)
  
  ## output path for saving datasets
  data_path <- "Data/Simulations/"
  
  ## initialize results data frame
  results <- NULL
  
  ## find which thing we are varying on x-axis
  loop_var <- max(length(N), length(P), length(K), length(sK), length(sI))

  ## loop over different values on x-axis
  for (i in 1:loop_var) {
    ## set parameter values
    if (length(N) > 1) { ## vary sample size
      setting <- "n"
      n <- N[i]; p <- P; k <- K; sig_k <- sK; sig_inx <- sI
      set_val <- n
      plot_title <- paste0("Simulations Varying Sample Size \n (p = ", p, ", k = ", k,
                           ", sig. k = ", sig_k, ", sig. interactions  = ", sig_inx, ")")
    } else if (length(P) > 1) { ## vary number of Xs
      setting <- "p"
      n <- N; p <- P[i]; k <- K; sig_k <- sK; sig_inx <- sI
      set_val <- p
      plot_title <- paste0("Simulations Varying Number of Features \n (n = ", n, ", k = ", k,
                           ", sig. k = ", sig_k, ", sig. interactions  = ", sig_inx, ")")
    } else if (length(K) > 1) { ## vary number of clusters
      setting <- "k"
      if(sI>0){
      n <- N; p <- P; k <- K[i]; sig_k <- ceiling(K[i]/10); sig_inx <- ceiling(K[i]/5)} else{
      n <- N; p <- P; k <- K[i]; sig_k <- K[i]/5; sig_inx <- 0
      
      print(sig_k)
      print(sig_inx)
        
      }
      set_val <- k
      plot_title <- paste0("Simulations Varying Number of Clusters \n (n = ", n, ", p = ", p,
                           ", sig. k = ", sig_k, ", sig. interactions  = ", sig_inx, ")")
    } else if (length(sK) > 1) { ## vary number of significant clusters
      setting <- "sig_k"
      n <- N; p <- P; k <- K; sig_k <- sK[i]; sig_inx <- sI
      set_val <- sig_k
      plot_title <- paste0("Simulations Varying Number of Clusters \n (n = ", n, ", p = ", p,
                           ", k = ", k, ", sig. interactions  = ", sig_inx, ")")
    } else { ## vary number of significant interaction terms
      setting <- "sig_inx"
      n <- N; p <- P; k <- K; sig_k <- sK; sig_inx <- sI[i]
      set_val <- sig_inx
      plot_title <- paste0("Simulations Varying Number of Clusters \n (n = ", n, ", p = ", p,
                           ", k = ", k, "sig. k = ", sig_k, ")")
    }
    
    cat("Running Loop", i, "\n")
    ## generate data
    print("n")
    print(n)
    print("p")
    print(p)
    print("k")
    print(k)
    print("sig_k")
    print(sig_k)
    print("sig_inx")
    print(sig_inx)
    
    sim_data <- genSimData(n = n, p = p, k = k, sig_k = sig_k, sig_inx = sig_inx)
    ## run model training and evaluation
    eval_models <- predError(x = sim_data$x,
                             y = sim_data$y,
                             spec = spec,
                             niter = niter,
                             fdr = fdr,
                             f_size = f_size,
                             ncore = ncore, 
                             out_file = paste0(out_path, "models_", setting, set_val, ".rda"),
                             k=k)
    ## save results
    eval_models$n <- n
    eval_models$p <- p
    eval_models$k <- k
    eval_models$sig_k <- sig_k
    eval_models$sig_inx <- sig_inx
    
    results <- rbind(results, eval_models)
    
    save(sim_data, results, file = paste0(out_path, "results_n", n, "_p", p, "_k", k, "_sigk", sig_k, "_siginx", sig_inx, ".rda"))
  }
  
  
  saveRDS(results, file = paste0(out_path, "results_n", n, "_p", p, "_k", k, "_sigk", sig_k, "_siginx", sig_inx,"_",seed,".rds"))
  ## plot results
  plot_df <- results %>%
    dplyr::mutate(method = as.factor(method))
  
  line_plot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = !!sym(setting), y = MSE, color = method)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(title = plot_title)

  
  ggplot2::ggsave(paste0(out_path, "results_n", n, "_p", p, "_k", k, "_sigk", sig_k, "_siginx", sig_inx,seed,".pdf"), line_plot, width = 8, height = 6, units = "in")
}
