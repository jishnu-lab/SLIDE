#' calculate Z matrix
#'
#' Calculate Z matrix by importing the PredZ function from EssReg package
#'
#' @param x_path a string that points to the x matrix. 
#' @param x_std a matrix of input data. 
#' @param all_latent_factors a object from getLatentFactors.
#' @param lf_path a string that points to the getLatentFactors result as an RDS object.
#' @param out_path a string that points to the folder where the z matrix will be outputted.
#' @return The z matrix calculated by the PredZ function
#' @export


calcZMatrix <- function(x_std, all_latent_factors, x_path = NULL, lf_path = NULL, out_path){
  if (!is.null(x_path)){
    x <- as.matrix(utils::read.csv(x_path, row.names = 1))
    x <- scale(x, T, T)
  } else if (!is.null(x_std)){
    x <- x_std
  } else{
    stop("No data matrix is inputted.")
  }
  
  if (!is.null(lf_path)){
    er_res <- readRDS(lf_path)
  } else if (!is.null(all_latent_factors)) {
    er_res <- all_latent_factors
  } else{
    stop("No latent factors are being inputted.")
  }
  
  z_matrix <- predZ(x, er_res)
  colnames(z_matrix) <- paste0("Z", c(1:ncol(z_matrix)))
  rownames(z_matrix) <- rownames(x)
  write.csv(z_matrix, paste0(out_path, "z_matrix.csv"), row.names = TRUE)
  return(z_matrix)
}
