#' Find most frequent variables.
#'
#' Finds the "joint" of the elbow if a histogram were to be made. Loops through all variables and selects those that appear in at least
#' 50% of the iterations or that are before the elbow. In this case, the elbow occurs when the proportion of selections in the iterations
#' changes by more than \code{spec}.
#' 
#' @param spec a numeric constant between 0 and 1; the change in frequency proportion
#' @param niter an integer; the number of times knockoffs was run
#' @param tab_data the frequency table made after running iterated knockoffs
#' @return a vector of names corresponding to the selected columns of \code{z} that are determined to be frequent enough
#' @export

findFrequent <- function(spec, niter, tab_data) {
  if (ncol(tab_data) == 0) {
    return (NULL)
  } else if (ncol(tab_data) == 1) {
    return (names(tab_data))
  }
  ## make frequency table into matrix
  freq_matrix <- data.frame(tab_data[1, ])
  colnames(freq_matrix) <- "freq"
  ## order by frequency
  freq_matrix <- freq_matrix %>% 
    dplyr::arrange(-freq)
  ## add most frequent to selected list
  selected <- rownames(freq_matrix)[1]
  last_freq <- freq_matrix[1, ] / niter
  ## iterate through remaining variables
  for (i in 2:nrow(freq_matrix)) {
    curr_freq <- freq_matrix[i, ] / niter
    freq_diff <- last_freq - curr_freq
    ## if a variable appears more than 50% of the time or has a frequency within
    ## spec of the last selected variable, add it to the selected list
    if (curr_freq > 0.5 || freq_diff < spec) {
      selected <- c(selected, rownames(freq_matrix)[i])
    } else {
      break
    }
    last_freq <- curr_freq
  }
  return (selected)
}
