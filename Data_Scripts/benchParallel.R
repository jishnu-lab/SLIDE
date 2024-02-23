benchParallel <- function(yaml_path, replicate) {
  ## process arguments
  er_input <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized
  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized
  z <- as.matrix(utils::read.csv(er_input$z_path, row.names = 1)) ## standardized
  
  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)
  
  if (er_input$eval_type == "auc") {
    if (length(er_input$y_order) > 2) {
      stop("y must be binary - please refactor")
    }
    y <- toCont(y, er_input$y_order)
    saveRDS(y, file = paste0(er_input$out_path, "benchViol_y_mapping.rds"))
    orig_y <- y$cat_y
    y <- y$cont_y
    ## set method levels
    method_order <- c("HER", "plainER", "lasso", "phate", "plsda", "pclr")
  } else {
    ## set method levels
    method_order <- c("HER", "plainER", "lasso", "phate", "plsr", "pcr")
  }
  
  ## create output directory
  new_dir <- paste0(er_input$out_path, "replicate", replicate, "/")
  dir.create(file.path(new_dir), showWarnings = F, recursive = T)
  sink(file = paste0(new_dir, "replicate_output.txt")) ## make replicate output file
  
  bench_rep <- NULL
  while (is.null(bench_rep)) {
    bench_rep <- benchCV(k = er_input$k,
                         x = x,
                         y = y,
                         z = z,
                         delta = er_input$delta,
                         eval_type = er_input$eval_type,
                         lambda = er_input$lambda,
                         out_path = er_input$out_path,
                         rep_cv = er_input$rep_cv,
                         method = er_input$method,
                         alpha_level = er_input$alpha_level,
                         thresh_fdr = er_input$thresh_fdr,
                         parallel = er_input$parallel,
                         y_order = er_input$y_order,
                         spec = er_input$spec,
                         niter = er_input$niter,
                         f_size = er_input$f_size,
                         fdr = er_input$fdr,
                         ncore = er_input$ncore,
                         rep = replicate)
  }
  
  sink() ## close file
  saveRDS(bench_rep, paste0(er_input$out_path, "benchmarks.RDS"))
}