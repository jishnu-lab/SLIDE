#' Find value of upsilon.
#'
#' Find upsilon by extracting the fitted values from a linear regression model.
#'
#' @param interaction_var a vector; the names of the significant interaction terms
#' @param interactions a matrix or data frame; the values of the significant interaction terms
#' @param z a matrix or data frame
#' @param marg_var an integer; the index of a marginal variable
#' @return a data frame with 1 column
#' @export

getUpsilon <- function(interaction_var, interactions, z, y, marg_var) {
  if (is.null(interaction_var)) { # if no interaction variables selected, just use marginal
    #cat("              no interaction vars . . . upsilon is marginal variable \n")
    ## create data frame for linear regression
    lin_reg_df <- z[, marg_var] %>%
      as.data.frame()
    colnames(lin_reg_df) <- marg_var
  } else { # if interactions are selected, create a new variable from fitted values of linear regression
    ## create data frame for linear regression
    lin_reg_df <- cbind(z[, marg_var], interactions[, interaction_var]) %>%
      as.data.frame()
    colnames(lin_reg_df) <- c(marg_var, interaction_var)
  }
  ## fit linear regression model and extract fitted values (these are upsilons)
  lin_reg <- stats::lm(y ~ ., data = lin_reg_df)
  lin_reg_fits <- lin_reg$fitted.values
  ## rename column
  upsilon <- data.frame(lin_reg_fits)
  colnames(upsilon) <- paste0("ups", marg_var)
  return (list("upsilon" = upsilon,
               "upsilon_mod" = lin_reg))
}
