#' Check Element In List.
#'
#' Check if an element is in a list. If it does exist in some group, return the group
#' index in that list and its sublist index. Otherwise, return \code{(0,0)}.
#'
#' @param element an element
#' @param groupList a list
#' @return group and sublist index if found, \code{(0, 0)} otherwise
#' @export

checkElement <- function(element, groupList) {
  for (i in 1:length(groupList)) {
    for (j in 1:length(groupList[[i]])) {
      if (element %in% groupList[[i]][[j]])
        return(c(i,j))
    }
  }
  return(c(0,0))
}
