#' Method A: Top Betas.
#'
#' Using the provided proportion, retrieve the top \eqn{\beta}s by absolute size and identify the corresponding latent \eqn{Z}s.
#'
#' @param betas a vector of values for the \eqn{\beta}s
#' @param top_prop a numeric constant: the proportion of \eqn{\beta}s to consider significant
#' @return a vector of \eqn{Z} indices
#' @export

getTopBetas <- function(betas, top_prop) {
  ## find number of betas in top_prop proportion
  num_k <- length(betas) * 0.01 * top_prop
  ## round to integer
  num_k <- round(num_k)
  ## make an even integer
  if (num_k %%2 != 0) {
    num_k <- num_k + 1
  }
  ## divide in half for positive and negative betas
  num_k <- num_k / 2
  ## get top positive betas
  pos_beta <- order(betas, decreasing = TRUE)[1:num_k]
  ## get top negative betas
  neg_beta <- order(betas, decreasing = FALSE)[1:num_k]
  ## concatenate
  top_betas <- c(pos_beta, neg_beta)
  return (top_betas)
}