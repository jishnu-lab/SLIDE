#' Read Essential Regression Results.
#'
#' Function to read all, pure, and mixed variables from ER results.
#'
#' @importFrom magrittr '%>%'
#' @param er_res a list output by the function \link[EssReg]{plainER} or \link[EssReg]{priorER}
#' @param y a vector of responses
#' @return a list including the clusters, pure variables, mixed variables, the significant
#' cluster indices, and the features found in the significant clusters
#' found by ER
#' @export

readER <- function(er_res) {
  A <- er_res$A
  #### get clusters with structure
  clusters <- list()
  for (i in 1:ncol(A)) {
    column <- A[,i]
    posInd <- which(column > 0)
    negInd <- which(column < 0)
    clusters[[i]] <- list(pos = posInd, neg = negInd)
  }
  #### get all features included in clusters
  samp_feats <- as.vector(unlist(clusters)) %>% unique()
  #### get pure variables
  pure_vars <- c()
  for (i in 1:ncol(A)) {
    pure_vars <- c(pure_vars, which(abs(A[ ,i]) == 1))
  }
  #### get mixed variables
  mix_vars <- setdiff(samp_feats, pure_vars)
  return(list("clusters" = clusters,
              "pure_vars" = pure_vars,
              "mix_vars" = mix_vars))
}
