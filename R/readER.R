#' Read Essential Regression Results.
#'
#' Function to read all, pure, and mixed variables from ER results.
#'
#' @importFrom magrittr '%>%'
#' @param er_res a list output by the function \link[EssReg]{plainER} or \link[EssReg]{priorER}
#' @param y a vector of responses
#' @return a list including the clusters, pure variables, mixed variables, the significant
#' @export

readER <- function(x, er_res) {
  var_names <- colnames(x)
  #### get clusters with structure
  clusters <- recoverGroup(er_res$A)
  new_clusters <- NULL
  for (i in 1:length(clusters)) {
    cluster_name <- paste0("Z", i)
    cluster <- clusters[[i]]
    if (length(cluster$pos) > 0) {
      names(cluster$pos) <- var_names[cluster$pos]
    } else {
      cluster$pos <- NULL
    }
    if (length(cluster$neg) > 0) {
      names(cluster$neg) <- var_names[cluster$neg]
    } else {
      cluster$neg <- NULL
    }
    new_clusters[[cluster_name]] <- cluster
  }
  clusters <- new_clusters
  #### get all features included in clusters
  samp_feats <- unlist(clusters) %>% unique()
  #### get pure variables
  pure_vars <- pureRowInd(er_res$A)
  names(pure_vars) <- var_names[pure_vars]
  #### get mixed variables
  mix_vars <- setdiff(samp_feats, pure_vars)
  if (length(mix_vars) > 0) {
    names(mix_vars) <- var_names[mix_vars]
  } else {
    mix_vars <- NULL
  }
  return(list("clusters" = clusters,
              "pure_vars" = pure_vars,
              "mix_vars" = mix_vars))
}
