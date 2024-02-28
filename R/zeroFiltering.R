#' Reduce sparsity of the data matrix by deleting samples and features that are expressed in few features and samples. 
#'
#' @param unfiltered_x the unfiltered x matrix
#' @param unfiltered_y the unfiltered y matrix
#' @param g_thresh features with number of zeros greater then g_thresh will be filtered out
#' @param c_thresh samples with number of zeros greater than c_thresh will be filtered out
#' @export


zeroFiltering <- function(x, y, g_thresh, c_thresh){
  #input_params <- yaml::yaml.load_file(yaml_path)
  
  #x <- as.matrix(utils::read.csv(input_params$x_path, row.names = 1))
  #y <- as.matrix(utils::read.csv(input_params$y_path, row.names = 1))
  cat("Original dataframe dimension is ", dim(x)[[1]], " by ", dim(x)[[2]], "\n")
  
  if (sum(row.names(x) == row.names(y)) != nrow(x)) {stop('Input data matrix and response does not have matching row names.')}
  
  #filtered out cells with num zeros greater than a certain threshold
  #thresh <- dim(x)[[2]] - c_thresh
  i <- rowSums(x == 0, na.rm=TRUE) <= c_thresh
  filtered <- x[i, ]
  
  #filtered out genes with num zeros greater than a certain threshold
  #thresh <- dim(data)[[1]] - g_thresh
  i <- colSums(x == 0, na.rm=TRUE) <= g_thresh
  filtered <- filtered[ ,i]
  
  cat("Filtered dataframe dimension is ", nrow(filtered), " by ", ncol(filtered), "\n")
  
  ### make sure Y has the same names
  filtered_y <- as.matrix(y[row.names(y) %in% row.names(filtered), ])
  colnames(filtered_y) = colnames(y)
  if (nrow(filtered_y) != nrow(filtered)) {stop("Number of samples in data matrix and response do not match post-filtering.")}
  
  
  # check if output path exists
  if (dir.exists(input_params$out_path)){
    cat("Populating filtered files to ", input_params$out_path, '.\n')
  } else{
    cat ("Output folder not found, creating at ", input_params$out_path, ".\n")
    dir.create(file.path(input_params$out_path), showWarnings = F, recursive = T)
  }
  
  write.csv(filtered, paste0(input_params$out_path, "/filtered_x.csv"))
  write.csv(filtered_y, paste0(input_params$out_path, "/filtered_y.csv"))
  input_params$x_path <-paste0(input_params$out_path, "/filtered_x.csv")
  input_params$y_path <-paste0(input_params$out_path, "/filtered_y.csv")
  
  return(list(input_params, filtered, filtered_y))
}

