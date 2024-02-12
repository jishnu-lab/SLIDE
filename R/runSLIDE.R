#' running the SLIDE function
#'
#' Run the SLIDE function and formulate the results and the key parameters into a object.

#' @importFrom foreach '%dopar%'
#' @param y_path a string that points to the response vector
#' @param y a matrix of the response vector
#' @param z_path a string that points to the z matrix, user can also input the matrix directly.
#' @param z_matrix the z matrix, user can also input a path to the csv file.
#' @param all_latent_factors a object contain all latent factors from getLatentFactors.
#' @param lf_path a string that points to the all_latent_factors as an RDS object.
#' @param method which method of SLIDE to run
#' @param do_interacts whether get interaction terms or not
#' @param fdr the fdr threshold
#' @param niter the number of iterations SLIDE run
#' @param spec the threshold to use to choose the latent factors
#' @param f_size the number of latent factors to split. 
#' @return A list that contains all key information from SLIDE.
#' @export


runSLIDE <- function(y, y_path = NULL, z_path = NULL, z_matrix, all_latent_factors, lf_path = NULL, method = 4, do_interacts=TRUE, fdr = 0.1, niter = 500, spec = 0.1, f_size = NULL){
  final_res <- NULL
  if (!is.null(y_path)){
    y <- as.matrix(utils::read.csv(y_path, row.names = 1))
  }
  
  if (!is.null(z_path)){
    z <- as.matrix(utils::read.csv(z_path, row.names = 1))
    } else {z = z_matrix}
  
  if (!is.null(lf_path)){
    all_latent_factors <- readRDS(lf_path)
    }
  
  if (is.null(f_size)){
    f_size <- calcDefaultFsize(y, all_latent_factors)
  }
  
  cat("f_size is set as ", f_size, "\n")
  if (!is.null(f_size)){
    tryCatch({
      SLIDE_res <- SLIDE(z, y, method = method, do_interacts = do_interacts, betas = NULL, top_prop = NULL, marginals = NULL,
                         spec = spec, fdr = fdr, niter = niter, elbow = FALSE, f_size = f_size, parallel = TRUE, ncore = 10)
    },
    error = function(e){
      warning("An error has occured with SLIDE. Re-running SLIDE with default f_size value.")
    },
    finally = {
      f_size = calcDefaultFsize(y, all_latent_factors)
      cat("Rerunning SLIDE with default f_size as ", f_size, '.\n')
      SLIDE_res <- SLIDE(z, y, method = method, do_interacts = do_interacts, betas = NULL, top_prop = NULL, marginals = NULL,
                         spec = spec, fdr = fdr, niter = niter, elbow = FALSE, f_size = f_size, parallel = TRUE, ncore = 10)
    })
  } else{
    SLIDE_res <- SLIDE(z, y, method = method, do_interacts = do_interacts, betas = NULL, top_prop = NULL, marginals = NULL,
                       spec = spec, fdr = fdr, niter = niter, elbow = FALSE, f_size = f_size, parallel = TRUE, ncore = 10)
  }
  
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
  ##final_res$marginal_vals <- as.numeric(gsub("^[zZ]","",SLIDE_res$marginal_vals))
  final_res$interaction <- data.frame(p1, p2)
  return(final_res)
}
