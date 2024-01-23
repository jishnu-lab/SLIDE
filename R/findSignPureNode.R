#' Estimate Signed Subpartitions Of Pure Nodes.
#'
#' Estimate the signed subpartition of pure node sets. If there is an element
#' of a list is empty, then a empty list will be put in that position.
#'
#' @param pure_list a list of pure node indices
#' @param sigma a sample correlation matrix of dimensions \eqn{p \times p}
#' @return a list of signed subpartitions of pure node indices
#' @export

findSignPureNode <- function(pure_list, sigma) {
  signed_pure_list <- list()
  for (i in 1:length(pure_list)) {
    pure_i <- sort(pure_list[[i]])   ### For simulation purpose only.
    if (length(pure_i) != 1) {
      first_pure <- pure_i[1]
      pos <- first_pure
      neg <- c()
      for (j in 2:length(pure_i)) {
        pure_j <- pure_i[j]
        if (sigma[first_pure, pure_j] < 0) {
          neg <- c(neg, pure_j)
        } else {
          pos <- c(pos, pure_j)
        }
      }
      if (length(neg) == 0) {
        neg <- list()
      }
      signed_pure_list[[i]] <- list(pos = pos, neg = neg)
    } else {
      signed_pure_list[[i]] <- list(pos = pure_i, neg = list())
    }
  }
  return(signed_pure_list)
}
