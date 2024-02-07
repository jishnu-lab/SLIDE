#' \eqn{\lambda} Cross-Validation.
#'
#' Use cross-validation to select \eqn{\lambda} for estimating \eqn{\Omega}.
#' Split the data into two parts and estimate \eqn{C} on both data sets. Then, for each
#' \eqn{lambda}, calculate \eqn{\Omega} on the first dataset and calculate the loss on the second dataset.
#' Find the lambda which minimizes \eqn{C \cdot \Omega - log(|\Omega|)}.
#'
#' @param x a data matrix of dimensions \eqn{n \times p}
#' @param fdr_entries a matrix of dimensions \eqn{p \times p} that contains a 1 in positions where the
#' entry is kept in the thresholded version of \eqn{\hat{\Sigma}} and a 0 if not
#' @param lambdas a vector of numerical constants over which to search for the optimal \eqn{\lambda}
#' @param AI the estimated matrix \eqn{A_I} of dimensions \eqn{p \times K}
#' @param pure_vec a vector of indices of pure variables
#' @param k number of folds for \eqn{k}-fold cross-validation
#' @return the selected optimal \eqn{\lambda}
#' @export

cvLambda <- function(x, fdr_entries, lambdas, AI, pure_vec, k) {
  #### split data matrix into training/validation sets
  samp_ind <- sample(nrow(x), floor(nrow(x) / 2))
  x_train <- x[samp_ind, ]
  x_val <- x[-samp_ind, ]

  #### calculate correlation matrices for training/validation sets
  sigma_train <- crossprod(x_train) / nrow(x_train)
  sigma_train <- sigma_train * fdr_entries #### control for FDR
  sigma_val <- crossprod(x_val) / nrow(x_val)
  sigma_val <- sigma_val * fdr_entries #### control for FDR

  #### calculate C
  C_train <- estC(sigma = sigma_train, AI = AI)
  C_val <- estC(sigma = sigma_val, AI = AI)

  loss <- c()
  for (i in 1:length(lambdas)) {
    Omega <- estOmega(lambdas[i], C_train)
    det_Omega <- det(Omega)
    loss[i] <- ifelse(det_Omega <= 0, Inf, sum(Omega * C_val) - log(det_Omega))
  }
  return(lambdas[which.min(loss)])
}
