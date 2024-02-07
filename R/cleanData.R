#' @export
cleanData <- function(xdata, ydata, mode = 0,
                      quantile_filter = 0.01, remove_zero_median_cols = T,
                      remove_zero_median_rows = T, scale_zeroes_directly = 0.0001,
                     er_input = NULL) {
#   if ( !is.null(xdata)) {

#     res = tryCatch({
#       list("xdata" = as.matrix(xdata), "ydata" = as.matrix(ydata))
#     }, error = function(e) {
#       cat("please pass matrix or dataframe inputs to xdata and ydata arguments\n")
#       return()
#     })
  # make sure in same order
  y_correct_order <- match(rownames(xdata), rownames(ydata))

  # we may not have labeled y values, if not just print a message
  if (any(!is.na(y_correct_order))) {
    ydata <- ydata[y_correct_order]
    # fix ydata rownames
    ydata <- as.matrix(ydata)
    rownames(ydata) <- rownames(xdata)

  } else {
    cat("Y does not have labels. Assuming Y and X have rows in the same order")
  }

  if (mode > 0) {

    if (scale_zeroes_directly > 0) {
      x_zeros = which(x == 0)

      x[x_zeros] = 0.0001
    }
    # remove zero SD too
    zero_sd <- which(apply(xdata, 2, sd) == 0)

    # check if we need to remove any columns
    if (length(zero_sd) > 0) {
      xdata <- xdata[, -zero_sd]
    }

    if (quantile_filter >= 0) {
      # filter quantiles and remove empty rows

      empty_cols <- which(apply(xdata, 2, median) == 0)

      colvars <- apply(xdata, 2, var)

      low_var_cols <- which(colvars < quantile(colvars, quantile_filter))

      if (remove_zero_median_cols) {
        remove_cols <- base::union(empty_cols, low_var_cols)
      } else {
        remove_cols <- low_var_cols
      }

      if (length(remove_cols > 0)) {
        xdata <- xdata[, -remove_cols]
      }

      if (remove_zero_median_rows) {
        # remove empty rows last
        empty_rows <- which(apply(xdata, 1, median) == 0)

        # check if we need to remove any rows
        if (length(empty_rows) > 0) {
          # remove rows
          xdata <- xdata[-empty_rows,]
          ydata <- ydata[-empty_rows]
        }
      }
    }
  }
  
  return(list("x" = xdata, "y" = ydata))
}
