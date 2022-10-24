plotSpecSims <- function(path) {
  ## list files in output directory
  all_results <- list.files(path = path,
                            pattern = ".rda",
                            recursive = TRUE)
  
  ## read in .rda results
  mean_fdrs <- NULL
  specs <- NULL
  for (i in 1:length(all_results)) {
    read_path <- paste0(path, all_results[[i]])
    ## find corresponding spec
    components <- strsplit(read_path, "_")[[1]]
    find_spec <- grep("spec", components)[2]
    spec_component <- components[find_spec]
    spec <- gsub(".*spec", "", spec_component)
    ## load in one result
    rep_res <- load(read_path)
    ## calculate mean fdr and record corresponding spec
    mean_fdr <- mean(fdrs)
    mean_fdrs <- c(mean_fdrs, mean_fdr)
    specs <- c(specs, spec)
  }
  
  ## make data frame for plotting
  plot_df <- data.frame(as.numeric(specs), as.numeric(mean_fdrs))
  ## rename data frame columns
  colnames(plot_df) <- c("spec", "mean fdr")
  
  ## make line plot
  line_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = `spec`, y = `mean fdr`)) +
    ggplot2::geom_point() + 
    ggplot2::geom_line() +
    ggplot2::labs(title = "Simulation Plot: Spec vs. Mean FDR")
  
  ## save to directory specified by path parameter
  ggplot2::ggsave(paste0(path, "simulation_spec_plot.pdf"), line_plot, width = 10, height = 6, units = "in")
}