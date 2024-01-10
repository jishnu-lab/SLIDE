#' Calculate Z matrix
#'
#' Calculate Z matrix by importing the PredZ function from EssReg package
#'
#' @param x_path a string that points to the x matrix. 
#' @param er_path a string that points to the final er result as an RDS object.
#' @param out_path a string that points to the folder where the z matrix will be outputted.
#' @return The z matrix calculated by the PredZ function
#' @export


CalcZMatrix <- function(x_path, er_path, out_path){
  x_path <- x_path
  x <- as.matrix(utils::read.csv(x_path, row.names = 1))
  x <- scale(x, T, T)
  er_res <- readRDS(er_path)
  
  z_matrix <- EssReg::predZ(x, er_res)
  colnames(z_matrix) <- paste0("Z", c(1:ncol(z_matrix)))
  write.csv(z_matrix, paste0(out_path, "z_matrix.csv"), row.names = TRUE)
  return(z_matrix)
}
