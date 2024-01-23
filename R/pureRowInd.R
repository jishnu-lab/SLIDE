#' Find Indices Of Pure Nodes.
#'
#' For given matrix \eqn{A}, find the row indices which correspond to pure nodes.
#'
#' @param A a matrix of dimensions \eqn{p \times K}
#' @return a vector of indices
#' @export

pureRowInd <- function(A) {
  pureVec <- c()
  for (i in 1:ncol(A)) {
    pureVec <- c(pureVec, which(abs(A[ ,i]) == 1))
  }
  return(pureVec)
}
