#' Check if the input and the yaml files are set correctly. 
#'
#' @param input_params a list of the input parameters from the YAML file.
#' @param pipeline a integer indicates which step of the pipeline is the YAML file for.
#' @export


checkDataParams <- function(input_params, pipeline = NULL){
  #input_params <- yaml::yaml.load_file(yaml_path)
  x <- as.matrix(utils::read.csv(input_params$x_path, row.names = 1))
  y <- as.matrix(utils::read.csv(input_params$y_path, row.names = 1))
  cat("Dimension of input X is:", dim(x)[1], "by", dim(x)[2], ".\n")
  cat("Dimension of input y is:", dim(y)[1], "by", dim(y)[2], ".\n")
  ##################################################################
  ##                          check data                          ##
  ##################################################################
  # check col names
  cat("Checking the format and dimensions of input data and response matrices... \n")
  if (is.null(names(as.data.frame(x))) == TRUE){stop("No feature names found for the input data (X), please add column names...")}
  if (is.null(names(as.data.frame(y))) == TRUE){stop("No column name found for the response data (y), please add column name...")}
  
  # check the dimensions  
  if (dim(x)[1] != length(y)) {stop("The number of samples in input data (X) and response data (y) does not match or data matrix needs to be transposed.")}
  if (dim(y)[1] <= 15) {stop("The number of samples are too low to recieve accurate results.")}
  if (dim(x)[2] >= 4000) {warning("The number of features are over 4000 which can result in long computational time, please reduce the number of features...")}
  
  # check the sample names
  if (sum(rownames(x) != row.names(y)) != 0) {warning("The samples names are not in the same order or do not match for the data matrix and the response vector.\n Note that the row orders of input data (X) and response data (y) needs to be the same.")}
  
  cat("Checking na values in the input data and response matrices... \n")
  # check na values
  if (anyNA(x)) {stop("Input data (X) has NA values...")}
  if (anyNA(y)) {stop("Input response (y) has NA values...")}
  if (anyNA(scale(x, T, T))) {stop("Input data (X) has features with zero standard deviation.")}
  
  # check levels in y
  if (length(unique(y)) < 2) { stop("There is only one unique value in the input response (y).")}
  
  ##################################################################
  ##                 check yaml for main pipeline                 ##
  ##################################################################
  cat("Checking if yaml file is correct for the input data and response matrices... \n")
  if ((length(unique(y)) == 2) & (input_params$y_factor == FALSE)) {
    stop("The response data (y) is binary, set f_factor is set to TRUE. ")
  }
  if ((length(unique(y)) == 2) & (input_params$eval_type == 'corr')) { 
    warning("Only 2 level in response data (y), evaluation type is set as correlation.")
  }
  if ((length(unique(y)) > 2) & (input_params$eval_type == 'auc')) { 
    stop("More than 2 level in response data (y), evaluation type is set as auc. Please change to corr.")
  }
  
  ##################################################################
  ##                          check yaml                          ##
  ##################################################################
  # if (is.null(input_params$k) == FALSE) {
  #   if ((length(y) <= 20) & (input_params$k != length(y))) 
  #   {
  #     warning("The number of data points are low and cross validation is not set to leave-one-out. We recommend setting k as number of samples...")
  #   }
  # }
  

  #################################################################
  ##                pipeline step specific checks                ##
  #################################################################
  if (is.null(pipeline) == FALSE){
    
    if (pipeline == 2){
      if ((length(input_params$lambda) == 0) | (is.null(input_params$lambda) == TRUE)) {
        stop("Lambda should be given multiple values when running pipeline 2...")
      }
      if ((length(input_params$delta) != 1) | (is.null(input_params$lambda) == TRUE)) {
        stop("Too many or no delta value is given...")
      }
    }
    
    if (pipeline == 3){
      if ((length(input_params$lambda) != 1) | (is.null(input_params$lambda) == TRUE)) {
        stop("Too many or no lambda value is given...")
      }
      if ((length(input_params$delta) != 1) | (is.null(input_params$lambda) == TRUE)) {
        stop("Too many or no delta value is given...")
      }
    }
  }
}