#' Sign Permutation.
#'
#' For a given integer \eqn{n}, generate its sign permutation.
#'
#' @param n an integer
#' @return matrix of dimensions \eqn{2^n \times n} containing all sign permutations
#' @export

signPerm <- function(n) {
  signPerms <- rep(1, n)
  for (i in 1:n) {
    allCombs <- gtools::combinations(n, i, 1:n)
    for (j in 1:nrow(allCombs)) {
      thisPerm <- rep(1, n)
      thisPerm[allCombs[j, ]] <- -1
      signPerms <- rbind(signPerms, thisPerm)
    }
  }
  return(signPerms)
}
