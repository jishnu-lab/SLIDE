#' Reduce sparsity of the data matrix by deleting samples and features that are expressed in few features and samples. 
#'
#' @param data the data matrix
#' @param g_thresh samples with number of zeros greater then g_thresh will be filtered out
#' @param c_thresh features with number of zeros greater then g_thresh will be filtered out
#' @export


ZeroFiltering <- function(data, g_thresh, c_thresh){
  cat("Original dataframe dimension is ", dim(data)[[1]], " by ", dim(data)[[2]], "\n")
  
  #filtered out cells with num zeros greater than a certain threshold
  thresh <- dim(data)[[2]] - g_thresh
  i <- rowSums(data == 0, na.rm=TRUE) <= thresh
  filtered <- data[i, ]
  
  #filtered out genes with num zeros greater than a certain threshold
  thresh <- dim(data)[[1]] - c_thresh
  i <- colSums(data == 0, na.rm=TRUE) <= thresh
  filtered <- filtered[ ,i]
  
  cat("Filtered dataframe dimension is ", nrow(filtered), " by ", ncol(filtered), "\n")
  return(filtered)
}
