#' Estimation \eqn{A_J} (Dantzig).
#'
#' Estimate the matrix \eqn{A_J} with a Dantzig-type estimator.
#'
#' @param C_hat a matrix
#' @param sigma_TJ a matrix
#' @param lambda \eqn{\lambda}, a positive constant
#' @param se_est_J the estimated standard errors of indices in \eqn{J}
#' @return an estimate of matrix \eqn{A_J}
#' @export

estAJDant <- function(C_hat, sigma_TJ, lambda, se_est_J) {
  AJ <- matrix(0, ncol(sigma_TJ), nrow(sigma_TJ))
  for (i in 1:ncol(sigma_TJ)) { ## loop through cols of sigma_TJ (correspond to rows of AJ)
    dantzig_sol <- dantzig(C_hat = C_hat,
                           target_vec = sigma_TJ[, i],
                           lambda = lambda * se_est_J[i])
    if (is.null(dantzig_sol)) {
      return (NULL)
    }
    AJ[i, ] <- dantzig_sol
    if (sum(abs(AJ[i, ])) > 1) {
      AJ[i, ] <- AJ[i, ] / sum(abs(AJ[i, ]))
    }
  }
  return(AJ)
}
