#' Indices To Names.
#'
#' Translate feature indices to feature names and vice versa.
#'
#' @param feats a vector of feature indices or names
#' @param all_names a vector of all feature names
#' @param to_ind a boolean indicating whether to go from indices to names (T) or names to indices (F)
#' @return a vector of feature indices or names depending upon value of \code{toInd}
#' @export

indName <- function(feats, all_names, to_ind = T) {
  if (to_ind) { ## if going from index to name
    where_feat <- match(feats, all_names)
    return (where_feat)
  } else { ## if going from name to index
    which_feat <- all_names[feats]
    return (which_feat)
  }
}
