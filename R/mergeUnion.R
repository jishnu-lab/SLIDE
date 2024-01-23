#' Merge Nodes (Union).
#'
#' Merge the provided vector of nodes into the list of vectors of nodes but use
#' union rather than the intersection.
#'
#' @param clusters a list of groups of node indices
#' @param group a vector of node indices
#' @return a list of the merged results
#' @export

mergeUnion <- function(clusters, group) {
  # merge the new group with the previous ones which have common nodes
  if (length(clusters) != 0) {
    common_groups <- sapply(clusters, FUN = function(x, y) {
      length(intersect(x, y))
    }, y = group)
    common_inds <- which(common_groups > 0)
    if (length(common_inds) > 0){
      new_group <- unlist(lapply(common_inds,
                                 FUN = function(x, y){y[[x]]}, y = clusters))
      remain_group <- lapply(which(common_groups == 0),
                             FUN = function(x, y){y[[x]]}, y = clusters)
      clusters <- append(remain_group, list(union(group, new_group)))
      return(clusters)
    }
  }
  clusters <- append(clusters, list(group))
  return(clusters)
}
