#' Linear Program Function.
#'
#' Solve the following using linear programming:
#' \deqn{\beta^{+} - \beta^{-} \leq \lambda + y - \beta^{+} + \beta^{-} \leq \lambda - y}
#' for semipositive \eqn{\beta^{+}, \beta^{-}}
#'
#' @param y response vector of dimension \eqn{K}
#' @param lambda \eqn{\lambda}, a positive constant
#' @return a vector of dimension \eqn{K}
#' @export

LP <- function(y, lambda) {
  K <- length(y)
  cvec <- rep(1, 2 * K)
  bvec <- c(lambda + y, lambda - y, rep(0, 2 * K))
  C <- matrix(0, K, 2 * K)
  for (i in 1:K) {
    indices <- c(i,i + K)
    C[i,indices] = c(1,-1)
  }
  Amat <- rbind(C, -C, diag(-1, nrow = 2 * K))
  LPsol <- linprog::solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  beta <- LPsol[1:K] - LPsol[(K + 1):(2 * K)]
  return (beta)
}
