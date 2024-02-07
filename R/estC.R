#' Estimate \eqn{C}.
#'
#' If diagonal = True, estimate only diagonal elements.
#' The estimate for \eqn{C}, \eqn{\hat{C}}, is given by the following:
#' \deqn{\hat{C}_{aa} = \frac{1}{|\hat{I}_a|(|\hat{I}_a| - 1)} \sum_{i,j \in \hat{I}_a, i \neq j} |\hat{\Sigma}_{ij}|}
#' \deqn{\hat{C}_{ab} = \frac{1}{|\hat{I}_a||\hat{I}_b|} \sum_{i \in \hat{I}_a, j \in \hat{I}_b} \hat{A}_{ia}\hat{A}_{ib}\hat{Sigma}_{ij}}
#' for each \eqn{a \in [ \hat{K} ]} and \eqn{b \in [ \hat{K} ], b \neq a}.
#'
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @param AI a matrix of dimensions \eqn{p \times K}
#' @return an estimate of \eqn{C} of dimensions \eqn{K \times K}
#' @export

estC <- function(sigma, AI) {
  K <- ncol(AI)
  C <- diag(0, K, K)
  for (i in 1:K) {
    group_i <- which(AI[, i] != 0)
    sigma_i <- as.matrix(abs(sigma[group_i, group_i]))
    tmp_entry <- sum(sigma_i) - sum(diag(sigma_i))
    C[i, i] <- tmp_entry / (length(group_i) * (length(group_i) - 1))
    if (i < K) {
      for (j in (i + 1):K) {
        group_j <- which(AI[, j] != 0)
        #### adjust the sign for each row
        sigma_ij <- AI[group_i, i] * as.matrix(sigma[group_i, group_j])
        sigma_ij <- t(AI[group_j, j] * t(sigma_ij))
        C[i, j] <- C[j, i] <- sum(sigma_ij) / (length(group_i) * length(group_j))
      }
    }
  }
  return(C)
}
