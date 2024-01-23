#' Find Group Structure.
#'
#' Recover group structure given \eqn{p \times K} matrix, \eqn{A} and perform thresholding.
#'
#' @param A a matrix estimated by Essential Regression
#' @return a list of group indices with sign subpartitions indicating the features included in a given cluster
#' @export

recoverGroup <- function(A) {
  Group <- list()
  for (i in 1:ncol(A)) {
    column <- A[,i]
    posInd <- which(column > 0)
    negInd <- which(column < 0)
    Group[[i]] <- list(pos = posInd, neg = negInd)
  }
  return(Group)
}
