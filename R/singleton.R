#' Check Singletone Node.
#'
#' Check if there exists an element of the given list has length equal to 1.
#'
#' @param pure_list estimated indices of pure nodes
#' @return true or false depending upon existence of singleton element
#' @export

singleton <- function(pure_list) {
  if (length(pure_list) == 0) {
    return(T)
  } else {
    ifelse(sum(sapply(pure_list, FUN = function(x) {length(x)}) == 1) > 0, T, F)
  }
}
