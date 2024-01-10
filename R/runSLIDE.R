#' running the SLIDE function
#'
#' Run the SLIDE function and formulate the results and the key parameters into a object.

#' @importFrom foreach '%dopar%'
#' @param y_path a string that points to the y vector
#' @param z_path a string that points to the z matrix, user can also input the matrix directly.
#' @param z_matrix the z matrix, user can also input a path to the csv file.
#' @param er_path a string that points to the final er result as an RDS object.
#' @param method which method of SLIDE to run
#' @param do_interacts whether get interaction terms or not
#' @param fdr the fdr threshold
#' @param niter the number of iterations SLIDE run
#' @param spec the threshold to use to choose the latent factors
#' @return A list that contains all key information from SLIDE.
#' @export


runSLIDE <- function(y_path, z_path = NULL, z_matrix, er_path, method = 4, do_interacts=TRUE, fdr = 0.1, niter = 500, spec = 0.1){
  final_res <- NULL
  y <- as.matrix(utils::read.csv(y_path, row.names = 1))
  if (is.null(z_path) == FALSE){z <- as.matrix(utils::read.csv(z_path, row.names = 1))}
  else{z = z_matrix}
  er_res <- readRDS(er_path)
  
  if (er_res$K <= 100){
    f_size = er_res$K
  }else (f_size = 100)
  
  cat("f_size is set as ", f_size, "\n")
  
  SLIDE_res <- SLIDE:::SLIDE(z, y, method = method, do_interacts = do_interacts, betas = NULL, top_prop = NULL, marginals = NULL,
                            spec = spec, fdr = fdr, niter = niter, elbow = FALSE, f_size = f_size, parallel = TRUE, ncore = 10)
  
  SLIDE_param <- c(method, spec, fdr, niter, f_size)
  names(SLIDE_param) <- c("method", "spec", "fdr", "niter", "f_size")

  p1 <- c()
  p2 <- c()
  for (i in SLIDE_res$interaction_vars){
    t <- gsub("Z", "", i)
    t <- strsplit(t, split = "[.]")[[1]]
    p1 <- append(t[[1]], p1)
    p2 <- append(t[[2]], p2)
  }
  p1 <- as.double(p1)
  p2 <- as.double(p2)
  final_res$SLIDE_res <- SLIDE_res
  final_res$SLIDE_param <- SLIDE_param
  final_res$marginal_vals <- SLIDE_res$marginal_vars
  final_res$interaction <- data.frame(p1, p2)
  return(final_res)
}
