#' Benchmarking methods with cross-validation wrapper.
#'
#' A wrapper for \code{benchCV} that configures the run by reading in a provided
#' .yaml file and creates the results directories.
#'
#' @param yaml_path a string path to the .yaml file with run configuration
#' @param replicate an integer; the replicate number used for output directory naming
#' @return none
#' @export

paperBench <- function(yaml_path, replicate) {
  ## process arguments

  er_input <- yaml::yaml.load_file(as.character(yaml_path))

  y <- as.matrix(utils::read.csv(er_input$y_path, row.names = 1)) ## not standardized

  x <- as.matrix(utils::read.csv(er_input$x_path, row.names = 1)) ## not standardized

  z <- as.matrix(utils::read.csv(er_input$z_path, row.names = 1)) ## standardized
  z <- scale(z,T,T)
  colnames(z)<- paste0("Z",1:ncol(z))

  ## create results directory (if not made yet)
  dir.create(file.path(er_input$out_path), showWarnings = F, recursive = T)



  ## create replicate output directory
  new_dir <- paste0(er_input$out_path, "replicate", replicate, "/")
  dir.create(file.path(new_dir), showWarnings = F, recursive = T)

  ## set up output report text file
  cat("Start of the analysis... \n")
  #sink(file = paste0(new_dir, "replicate_output.txt")) ## make replicate output file


  #
  er_input$std_cv <- T
  er_input$std_y <- T
  er_input$benchmark <- T
  er_input$permute <- T


  benchCV2(k = er_input$k,
           x = x,
           y = y,
           z=z,
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
           rep = replicate,
           benchmark = F,
           niter=er_input$niter,
           spec=er_input$spec,
           fdr=er_input$fdr,
           f_size=er_input$f_size,
           parallel = er_input$parallel,
           ncore=20)

  print(replicate)

}
