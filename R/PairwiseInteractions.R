#' PairwiseInteractions
#'
#' @param index_list a string that points to the x matrix. 
#' @param mat a string that points to the final er result as an RDS object.
#' @return a matrix of pairwise interactions
#' @export


PairwiseInteractions <- function(index_list, mat) {
  num_cols <- ncol(mat)
  index_combinations <- expand.grid(index_list, seq_len(num_cols))
  col_names <- paste0(colnames(mat)[index_combinations[, 1]], ".", colnames(mat)[index_combinations[, 2]])
  interaction_mat <- mat[, index_combinations[, 1]] * mat[, index_combinations[, 2]]
  colnames(interaction_mat) <- col_names
  return(list(interaction=as.data.frame(interaction_mat)))
}
