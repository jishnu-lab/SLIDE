#' Significant \eqn{\beta}s.
#'
#' Find top \eqn{\beta}s according to cutoff percentage.
#'
#' @importFrom magrittr '%>%'
#' @param betas a vector of estimates for \eqn{\beta}s
#' @param cutoff a numeric constant (proportion) for thresholding (default = 0.10)
#' @return a list of beta indices organized by sign
#' @export

sigBetas <- function(betas, cutoff) {
  num_betas <- length(betas)
  betas <- data.frame(seq(1, num_betas),
                      betas)
  colnames(betas) <- c("ind", "beta")
  num_sel <- ceiling(num_betas * cutoff / 2)
  neg <- betas %>%
    dplyr::filter(beta < 0)
  neg <- neg[order(neg$beta, decreasing = FALSE), ]
  pos <- betas %>%
    dplyr::filter(beta > 0)
  pos <- pos[order(pos$beta, decreasing = TRUE), ]

  #### get rid of rownames
  rownames(pos) <- NULL
  rownames(neg) <- NULL

  #### select positive betas
  num_sel <- min(num_sel, nrow(neg))
  pos_sel <- pos[1:num_sel, ]

  #### select negative betas
  num_sel <- min(num_sel, nrow(pos))
  neg_sel <- neg[1:num_sel, ]

  sig_betas <- list("pos_betas" = pos,
                    "neg_betas" = neg,
                    "pos_sig" = pos_sel[, 1],
                    "neg_sig" = neg_sel[, 1])
  return (sig_betas)
}
