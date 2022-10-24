#' Find significant upsilons.
#'
#' Use second order knockoffs to find the significant upsilons according to length and frequency.
#'
#' @param upsilon a data frame; the values of upsilon
#' @param y a vector; response values
#' @param spec a numeric constant; the proportion of iterations that a variable must be selected in to be considered significant
#' @param fdr a numeric constant; the target false discovery rate for knockoffs
#' @param niter a numeric constant; the number of times to run knockoffs
#' @param f_size a numeric constant; the target size for the subset (number of columns in each subset)
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a data frame
#' @export

testUpsilon <- function(upsilon, y,  spec = 0.1, fdr = 0.1, niter = 1000, f_size = 100, parallel = TRUE, elbow = FALSE) {
  if (ncol(upsilon) == 1) { # if only one upsilon, just use it
    cat("            only one upsilon . . . using without testing \n")
    print("uplsion error cheking: \n")
    print(colnames(upsilon))
    return (upsilon)
  } else {
    ## find significant upsilons
    selected <- selectShortFreq(z = upsilon,
                                y = y,
                                niter = niter,
                                spec = spec,
                                elbow = elbow,
                                fdr = fdr,
                                f_size = f_size,
                                parallel = parallel)
    if (is.null(selected)) { # if none selected, exit
      cat("      no models selected . . . \n")
      return (NULL)
    }
    selected_upsilon <- upsilon[, selected,drop=F]
  }
  return (selected_upsilon)
}