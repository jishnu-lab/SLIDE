#' Find the best results from \code{niter} iterations of knockoffs.
#'
#' Find the variables that are selected at least \code{spec} proportion of the \code{niter} of the iterations then find the shortest results from all \code{niter} that
#' have the maximum overlap with this list of selected variables.
#' 
#' @param n an integer; sample size
#' @param freq_vars a vector; the names of the most frequent variables
#' @param selected_list a list of the selected variables each iteration of knockoffs
#' @return a vector containing the names of the variables found in the results of the best iteration
#' @export

findOptIter <- function(freq_vars, selected_list) {
  ## look through all iterations and find the maximal number of overlapping variables
  ## between an iteration's selected variables and the most frequent variables
  mm <- max(unlist(lapply(selected_list, function(x){ sum(x %in% freq_vars) })))
  ## find which iterations have this maximal overlap found above
  max_overlap_ind <- which(unlist(lapply(selected_list, function(x){ sum(x %in% freq_vars) })) == mm)
  ## find the length of each iteration in max_overlap_ind
  overlap_list_len <- sapply(max_overlap_ind, function(x){ length(selected_list[[x]]) })
  ## find the shortest iteration
  selected_run <- max_overlap_ind[which.min(overlap_list_len)]
  ## select this iteration's results
  selected_vars <- selected_list[[selected_run]]
  return (selected_vars)
}