#' Check if the input and the yaml files are set correctly. 
#'
#' @param yaml_path a string that points to the yaml file.
#' @param x_path a string that points to the x matrix.
#' @param y_path a string that points to the y vector.
#' @param pipeline a integer indicates which step of the pipeline is the yaml file for.
#' @export


CheckDataParams <- function(yaml_path, x_path, y_path, pipeline){
  x <- as.matrix(utils::read.csv(x_path, row.names = 1))
  y <- as.matrix(utils::read.csv(y_path, row.names = 1))
  er_input <- yaml::yaml.load_file(yaml_path)
  
  ##################################################################
  ##                          check data                          ##
  ##################################################################
  # check col names
  if (is.null(names(as.data.frame(x))) == TRUE){stop("No feature names found for the X, please add column names...")}
  if (is.null(names(as.data.frame(y))) == TRUE){stop("No column name found for the Y, please add column name...")}
  # check the dimensions  
  if (dim(x)[1] != length(y)) {stop("The number of samples of X and Y does not match or data matrix needs to be transposed.")}
  # check the sample names
  if (sum(rownames(x) != row.names(y)) != 0) {warning("The samples names are not in the same order or do not match for the data matrix and the response vector. 
                                                      Note that the row orders of X and Y needs to be the same.")}
  
  ##################################################################
  ##                          check yaml                          ##
  ##################################################################
  if ((length(y) <= 20) & (er_input$k != length(y))) 
  {
    warning("The number of data points are low and cross validation is not set to leave-one-out. We recommend setting k as number of samples...")
  }
  if ((length(unique(y)) == 2) & (er_input$y_factor == FALSE)) {
    warning("The Y is binary, yet f_factor is set to FALSE. ")
  }
  
  #################################################################
  ##                pipeline step specific checks                ##
  #################################################################
  if (pipeline == 2){
    if ((length(er_input$lambda) == 0) | (is.null(er_input$lambda) == TRUE)) {
      stop("Lambda should be given multiple values when running pipeline 2...")
    }
    if ((length(er_input$delta) != 1) | (is.null(er_input$lambda) == TRUE)) {
      stop("Too many or no delta value is given...")
    }
  }
  if (pipeline == 3){
    if ((length(er_input$lambda) != 1) | (is.null(er_input$lambda) == TRUE)) {
      stop("Too many or no lambda value is given...")
    }
    if ((length(er_input$delta) != 1) | (is.null(er_input$lambda) == TRUE)) {
      stop("Too many or no delta value is given...")
    }
  }
}