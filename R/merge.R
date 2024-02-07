#' Merge Nodes (Intersection). NOT USED ANYMORE
#'
#' Merge a group with other groups containing common nodes.
#'
#' @param clusters a list of groups of node indices
#' @param group a vector of node indices
#' @return a list of the merged results
#' @export

merge <- function(clusters, group) {
  # merge the new group with the previous ones which have common nodes
  if (length(clusters) != 0) {
    for (i in 1:length(clusters)) {
      common_nodes <- intersect(clusters[[i]], group)
      if (length(common_nodes) != 0) {
        clusters[[i]] <- common_nodes
        return(clusters)
      }
    }
  }
  clusters <- append(clusters, list(group))
  return(clusters)
}
