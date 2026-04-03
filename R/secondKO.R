#' Run aggregated second order knockoffs.
#'
#' Iteratively run knockoffs and report variable selection frequencies.
#'
#' Pre-computes the covariance matrix and SDP solution once, then reuses
#' the cached \code{diag_s} across all \code{niter} iterations. Only the
#' random knockoff sampling step is repeated, avoiding redundant O(p^3)
#' SDP solves.
#'
#' @param z a matrix or data frame of measured data of dimensions \eqn{n} by \eqn{p}
#' @param y a matrix or data frame of responses of dimensions \eqn{n} by 1
#' @param niter a numeric constant; the number of times to run knockoffs
#' @param fdr a numeric constant; the target false discovery rate for knockoffs
#' @param parallel a boolean flag; run iterations in parallel or sequentially
#' @return a list including the names of the selected columns of \code{z} that are determined to be significant and a table of their selection frequencies
#' @export

secondKO <- function(z, y, niter = 100, fdr = 0.2, parallel = TRUE) {
  ## Z and Y sanity check
  if (is.null(z) || is.null(y)) {
    stop("Z and Y must not be empty.")
  }

  ## Pre-compute the SDP solution ONCE (the expensive O(p^3) part).
  ## knockoff::create.second_order internally:
  ##   1. mu    <- colMeans(X)
  ##   2. Sigma <- cov(X)
  ##   3. diag_s <- solve SDP on Sigma   (expensive, deterministic)
  ##   4. sample knockoffs from diag_s   (cheap, random)
  ## Steps 1-3 are identical across all niter iterations, so we cache them
  ## and only repeat step 4.

  z_mat <- as.matrix(z)
  mu <- colMeans(z_mat)
  Sigma <- cov(z_mat)

  ## Verify Sigma is positive-definite; apply shrinkage if not
  if (!knockoff:::is_posdef(Sigma)) {
    if (!requireNamespace("corpcor", quietly = TRUE)) {
      stop("corpcor is not installed", call. = FALSE)
    }
    Sigma <- tryCatch(
      suppressWarnings(
        matrix(as.numeric(corpcor::cov.shrink(z_mat, verbose = FALSE)),
               nrow = ncol(z_mat))
      ),
      error = function(e) {
        stop("Shrinkage estimation of covariance matrix failed: ", e$message,
             call. = FALSE)
      }
    )
  }

  ## Solve SDP once with fallback to equi method
  diag_s <- tryCatch(
    knockoff::create.solve_asdp(Sigma),
    error = function(e) {
      warning("ASDP solver failed, falling back to equi method: ", e$message,
              immediate. = TRUE)
      knockoff::create.solve_equi(Sigma)
    }
  )

  ## Build a closure that generates knockoffs from the cached diag_s.
  ## knockoff::create.gaussian accepts a diag_s argument that bypasses
  ## the internal SDP solve, so we leverage the package's own implementation.
  create_knockoffs_cached <- function(X) {
    knockoff::create.gaussian(X, mu, Sigma, diag_s = diag_s)
  }

  if (parallel == TRUE) { ## parallel computing
    selected_list <- foreach::foreach(i = 1:niter) %dopar% {
      result <- knockoff::knockoff.filter(X = z,
                                          y = as.matrix(y),
                                          knockoffs = create_knockoffs_cached,
                                          statistic = knockoff::stat.glmnet_lambdasmax,
                                          offset = 0,
                                          fdr = fdr)
      names(result$selected)
    }
  } else { ## sequential computing
    selected_list <- foreach::foreach(i = 1:niter) %do% {
      result <- knockoff::knockoff.filter(X = z,
                                          y = y,
                                          knockoffs = create_knockoffs_cached,
                                          statistic = knockoff::stat.glmnet_lambdasmax,
                                          offset = 0,
                                          fdr = fdr)
      names(result$selected)
    }
  }

  ## create frequency table
  tab_data <- table(unlist(selected_list))
  return(list(selected_list = selected_list,
              tab_data = as.matrix(t(tab_data))))
}
