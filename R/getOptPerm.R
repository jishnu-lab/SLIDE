#' Find Best Column Permutation Matrices.
#'
#' For the given pure node set, find the best column permutation matrix
#' \eqn{A_{perm}} of \eqn{A} such that \eqn{|A_{perm} - B|_F} is minimal for the target \eqn{B}.
#'
#' @param A a matrix
#' @param B a matrix of same dimensions as \eqn{A}
#' @return a list of column and sign permutation matrices
#' @export

getOptPerm <- function(A, B) {
  K <- ncol(A)
  allPerm <- gtools::permutations(K, K, 1:K)
  signPerms <- signPerm(K)
  prevLoss <- norm(A - B, "F")
  optPerm <- vector("list", length <- 2)
  for (i in 1:nrow(allPerm)) {
    # skip the identity permutation
    permi <- methods::as(as.integer(allPerm[i, ]), "pMatrix")
    permutatedA <- A %*% permi
    for (j in 1:nrow(signPerms)) {
      signPermj <- signPerms[j, ]
      # perform sign permutation along columns of newA
      newA <- t(t(permutatedA) * signPermj)
      currLoss <- norm(newA - B, "F")
      if (currLoss <= prevLoss) {
        optPerm[[1]] <- permi
        optPerm[[2]] <- signPermj
        prevLoss <- currLoss
      }
    }
  }
  return(optPerm)
}
