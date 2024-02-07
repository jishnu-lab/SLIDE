#' Categorical To Continuous Data.
#'
#' Transform categorical response vector into continuous values.
#'
#' @param y a response vector of dimension \eqn{n} of categorical values
#' @param order a vector indicating the level ordering of \eqn{y}
#' @return a list containing the old (categorical) \eqn{y} and the new (continuous) \eqn{y}
#' @export

toCont <- function(y, order = NULL) {
  y <- unlist(y)
  if (is.null(order)) {
    uniq_vals <- unique(y)
  } else {
    uniq_vals <- order
  }
  num_uniq <- length(uniq_vals)
  cont_vals <- seq(0, (num_uniq - 1))
  corresp_vals <- cbind(uniq_vals, cont_vals)
  new_y <- plyr::mapvalues(y, from = uniq_vals, to = cont_vals)
  return (list("cat_y" = y,
               "cont_y" = as.numeric(new_y),
               "mapping" = corresp_vals))
}
