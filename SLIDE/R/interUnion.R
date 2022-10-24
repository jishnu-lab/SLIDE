#' Construct interaction terms according to union logic. 
#'
#' Given a list of significant marginal variables, construct their interaction terms with all other
#' variables without repetitions (i.e. if \eqn{Z1.Z2} is already made, do not create \eqn{Z2.Z1}).
#' 
#' @param marginal_vars a vector of significant marginal variable indices
#' @param z a matrix or data frame; the values for all variables
#' @return a list including the significant marginal variables that were used as a vector and a matrix
#' containing the values of the constructed interaction terms
#' @export

interUnion <- function(marginal_vars, z) {
  ## cast data to matrix
  z <- as.matrix(z)

  ## give the data matrix column names if they do not already exist
  if (is.null(colnames(z)) || !grepl("Z",colnames(z)[1])) {
    colnames(z) <- paste0('Z', seq(1, ncol(z)))
  }
  if(length(marginal_vars)==0){stop("marginal vars cannot be NA!")}
  ## extract the values of the supplied significant marginals
  # print("where the error happnes")
  # print(colnames(z))
  # print(marginal_vars)
  # print("where the error happnes")
  # print(colnames(z))
  # print(marginal_vars)
  
  
  
  if(grepl("Z",marginal_vars[1]) & grepl("Z",colnames(z)[1])){margs <- z[,  marginal_vars]}
  if(!grepl("Z",marginal_vars[1]) & !grepl("Z",colnames(z)[1])){margs <- z[,  marginal_vars]} 
  if(!grepl("Z",marginal_vars[1]) & grepl("Z",colnames(z)[1])){margs <- z[, paste0('Z', marginal_vars)]} 
  if(grepl("Z",marginal_vars[1]) & !grepl("Z",colnames(z)[1])){margs <- z[, gsub('Z',"",marginal_vars)]}
  
  
  
  ## form them into a matrix
  margs <- matrix(margs, nrow = nrow(z))
  ## extract their indices from their names
  marg_var_inds <- gsub("Z", "", marginal_vars)
  ## rename the columns of this matrix
  colnames(margs) <- paste0("Z", marg_var_inds)

  ## initialize results vectors
  used <- NULL # vector of already used marginals
  inter_terms <- NULL
  ## loop through supplied significant marginals
  for (i in 1:ncol(margs)) {
    ## get the values of one marginal
    marg_col <- margs[, i]
    ## get this marginal's name
    colname <- colnames(margs)[i]
    ## add current marginal to vector of used
    used <- c(used, colname)
    cross_term_inds <- colnames(z)
    ## remove marginals that have been used
    cross_term_inds <- cross_term_inds[!cross_term_inds %in% used]
    if (nrow(margs) == 1) { ## if only one significant marginal
      ## create interaction term
      add_terms <- marg_col * as.matrix(z[, cross_term_inds])
      ## form into matrix
      add_terms <- matrix(add_terms, nrow = 1)
    } else {
      ## create interaction terms (automatically matrix)
      add_terms <- diag(marg_col) %*% as.matrix(z[, cross_term_inds])
    }
    ## rename matrix columns
    colnames(add_terms) <- paste0(colname, ".", cross_term_inds)
    ## add interactions with this marginal to matrix of all interactions
    inter_terms <- cbind(inter_terms, add_terms)
  }

  if(length(grep("NA",inter_terms))!=0){stop("marginal vars cannot be NA!")}
  return(list("marginal" = margs,
              "interaction" = inter_terms))
}
