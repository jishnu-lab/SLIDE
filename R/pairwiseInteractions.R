#' PairwiseInteractions
#'
#' @param index_list a string that points to the x matrix.
#' @param mat a string that points to the final er result as an RDS object.
#' @return a matrix of pairwise interactions
#' @export


pairwiseInteractions <- function(index_list, mat) {
  num_cols <- ncol(mat)
  index_combinations <- expand.grid(seq_len(num_cols),index_list,stringsAsFactors=F)
  temp <- index_combinations$Var1
  index_combinations$Var1 <- index_combinations$Var2
  index_combinations$Var2 <- temp

  col_names <- paste0(colnames(mat)[as.numeric(index_combinations[, 1])], ".", colnames(mat)[as.numeric(index_combinations[, 2])])
  interaction_mat <- mat[, as.numeric(index_combinations[, 1])] * mat[, as.numeric(index_combinations[, 2])]
  if(is.null(dim(interaction_mat))){
    interaction_mat <- matrix(interaction_mat,nrow=1)}

  colnames(interaction_mat) <- col_names
  return(list(interaction=as.data.frame(interaction_mat)))
}
