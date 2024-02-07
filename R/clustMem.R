#' Feature Cluster Membership.
#'
#' Function to find the clusters of which a given feature is a member.
#'
#' @importFrom magrittr '%>%'
#' @param feat a feature index
#' @param er_res a list output by the function \link{plainER()}
#' @return a vector of cluster indices
#' @export

clustMem <- function(feat, er_res) {
  er_feats <- readER(er_res)$clusters
  feats_unlist <- lapply(X = er_feats, FUN = unlist)
  feat_mem <- lapply(X = feats_unlist, FUN = function(x){
    if (feat %in% x) {
      return (TRUE)
    }
    return (FALSE)
  })
  membership <- which(feat_mem == TRUE)
  return (membership)
}
