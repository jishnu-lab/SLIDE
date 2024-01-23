
# usage: result = checkData(yaml_path). Filtered data can be found with result$new_x
#' @export
checkData <- function(x_input = NULL, y_input = NULL) {

  x <- as.matrix(x_input)
  y <- as.matrix(y_input)

  error_messages <- list()
  fixed_data <- list("x" = x, "y" = y)

  # check for incorrectly ordered x and y
  y_correct_order <- match(rownames(x), rownames(y))
  new_y = as.matrix(y[y_correct_order])
  rownames(new_y) <- rownames(x)

  if (all(rownames(new_y) != rownames(y))) {
    error_messages <- append(error_messages, "x and y need to have rownames in the same order \n")
    fixed_data$new_y <- new_y
    fixed_data$new_x <- x
  }

  # check x and y dimensions
  if (nrow(x) != nrow(y)) {
    error_messages = append(error_messages, "x and y need to have same number of rows \n")
  }

  # check x and y for NA values
  if (anyNA(x)) {
    error_messages <- append(error_messages, "x data has NA values \n")
  }
  if (anyNA(y)) {
    error_messages <- append(error_messages, "y has NA values \n")
  }

  # check for NA values in centered/scaled x
  if (anyNA(scale(x, T, T))) {
    zero_cols <- which(apply(x, 2, std) == 0)

    error_messages <- append(error_messages,
                            cat("columns have features with zero standard deviation: ",
                                zero_cols, "\n"))

    filtered_x <- x[, -zero_cols]
    fixed_data$new_x <- filtered_x
  }

  if (length(error_messages) == 0) {
    cat("no NA values or dimension mismatches in input data \n")
  } else {
    for (m in error_messages) {
      cat(m, "\n")
    }
  }

  if (length(fixed_data) > 2) {
    cat(" try running pipeline with returned x and y matrices ")
  }
  return(fixed_data)
}
