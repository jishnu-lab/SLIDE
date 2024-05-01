#' Reduce sparsity of the data matrix by deleting samples and features that are expressed in few features and samples. 
#'
#' @param unfiltered_x the unfiltered x matrix
#' @param unfiltered_y the unfiltered y vector
#' @param col_thresh filtered out columns with number of zeros greater than col_thresh, if entered as nrow of x, no columns will be deleted.
#' @param row_thresh filtered out rows with number of zeros greater than row_thresh, if entered as ncol of x, no rows will be deleted.
#' @return a list including the filtered x matrix and y vector
#' @export


zeroFiltering <- function(x, y, col_thresh, row_thresh){
  cat("Original dataframe dimension is ", dim(x)[[1]], " by ", dim(x)[[2]], "\n")
  
  if (sum(row.names(x) == row.names(y)) != nrow(x)) {stop('Input data matrix and response does not have matching row names.')}
  
  #filtered out cells with num zeros greater than a certain threshold
  #thresh <- dim(x)[[2]] - row_thresh
  i <- rowSums(x == 0, na.rm=TRUE) <= row_thresh
  filtered <- x[i, ]
  
  #filtered out genes with num zeros greater than a certain threshold
  #thresh <- dim(data)[[1]] - col_thresh
  i <- colSums(x == 0, na.rm=TRUE) <= col_thresh
  filtered <- filtered[ ,i]
  
  cat("Filtered dataframe dimension is ", nrow(filtered), " by ", ncol(filtered), "\n")
  
  ### make sure Y has the same names
  filtered_y <- as.matrix(y[row.names(y) %in% row.names(filtered), ])
  colnames(filtered_y) = colnames(y)
  if (nrow(filtered_y) != nrow(filtered)) {stop("Number of samples in data matrix and response do not match post-filtering.")}
  rownames(filtered_y) = rownames(filtered)
  
  filtered_data = list()
  filtered_data$filtered_x <- filtered
  filtered_data$filtered_y <- filtered_y
  return(filtered_data)
}
