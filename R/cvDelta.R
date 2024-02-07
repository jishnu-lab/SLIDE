#' \eqn{\delta} Cross-Validation.
#'
#' Cross validation for choosing \eqn{\delta}. For each delta from the given grids,
#' first split the data into two data sets. Obtain \eqn{I}, \eqn{A_I} and \eqn{C} from data set 1.
#' Then calculate \eqn{A_I \cdot C \cdot A_I^\top} and choose \eqn{\delta} which minimizes the criterion
#' \eqn{A_I \cdot C \cdot A_I^\top - \Sigma(\text{data set 2})}.
#'
#' @importFrom magrittr '%>%'
#' @importFrom foreach '%dopar%'
#' @param raw_x data matrix of dimensions \eqn{n \times p}, not standardized
#' @param fdr_entries a matrix of dimensions \eqn{p \times p} that contains a 1 in positions where the
#' entry is kept in the thresholded version of \eqn{\hat{\Sigma}} and a 0 if not
#' @param deltas_scaled a vector of numerical constants over which to perform the search for the optimal \eqn{\delta}
#' @return the selected optimal \eqn{\delta}
#' @export

cvDelta <- function(raw_x, fdr_entries, deltas_scaled) {
  #### get data matrix dimensions
  n <- nrow(raw_x); p <- ncol(raw_x)

  #### split data into training/validation sets
  samp_ind <- sample(1:n)
  samp_ind <- samp_ind[1:ceiling(2 * n / 3)]
  x_train <- raw_x[samp_ind, ]
  x_val <- raw_x[-samp_ind, ]

  #### calculate the sample correlation matrix for training set
  sigma_train <- cor(x_train);
  sigma_train <- sigma_train * fdr_entries #### control for FDR
  diag(sigma_train) <- 0 #### set diagonal to 0 for findRowMax()
  se_est <- apply(x_train, 2, stats::sd)

  #### calculate the sample correlation matrix for validation set
  sigma_val <- cor(x_val)
  sigma_val <- sigma_val * fdr_entries #### control for FDR

  result_max <- findRowMax(abs(sigma_train))
  max_vals <- result_max$max_vals
  max_inds <- result_max$max_inds

  loss <- c()
  for (i in 1:length(deltas_scaled)) {
    result_fitted <- calFittedSigma(sigma = sigma_train,
                                    delta = deltas_scaled[i],
                                    max_vals = max_vals,
                                    max_inds = max_inds,
                                    se_est = se_est)
    fit_sigma <- result_fitted$fit_sigma
    pure_vec <- result_fitted$pure_vec

    if (is.null(dim(fit_sigma)) && fit_sigma == -1) {
      loss[i] <- Inf
    } else {
      denom <- length(pure_vec) * (length(pure_vec) - 1)
      loss[i] <- 2 * offSum(sigma_val[pure_vec, pure_vec], fit_sigma, se_est[pure_vec]) / denom
    }
  }
  return(deltas_scaled[which.min(loss)])
}
