#' Output the top features in each latent factors
#'
#' Output the top features in each latent factors by taking a union of features with highest A loadings and spearman correlations or univariate AUCs
#'
#' @param x_path a string that points to the x matrix.
#' @param y_path a string that points to the y vector.
#' @param er_path a string that points to the final er result as an RDS object.
#' @param out_path a string that points to a folder where txt files of each latent factors will be outputted.
#' @param SLIDE_res the object outputted by the runSLIDE function. 
#' @param num_top_feats the number of top features to get from each criteria such as A loadings or AUCs.
#' @param condition either "corr" or "auc", indicates which one to use to get top features. 
#' @return The SLIDE_res object with feature results added to one of the slots.
#' @export


GetTopFeatures <- function(x_path, y_path, er_path, out_path, SLIDE_res, num_top_feats = 10, condition){
  
  x <- as.matrix(utils::read.csv(x_path, row.names = 1))
  y <- as.matrix(utils::read.csv(y_path, row.names = 1))
  er_res <- readRDS(er_path)
  ks <- union(union(unique(SLIDE_res$interaction$p1), unique(SLIDE_res$interaction$p2)), unique(SLIDE_res$marginal_vals))
  
  if (is.null(ks) == TRUE){stop('The SLIDE_res input is not formatted correctly. Please re-run the runSLIDE function...')}
  if ("auc" == condition & length(unique(y[, 1])) != 2){stop('Only 2 levels allowed for y when condition = "auc".')}
  
  A <- er_res$A[, ks]
  
  gene_names <- colnames(x)
  
  temp <- NULL
  
  for (i in 1:ncol(A)){
    AUCs <- c()
    signs <- c()
    corrs <- c()
    idx <- which(A[, i] != 0)
    A_loading <- abs(A[, i][-which(A[, i] == 0)])
    names <- gene_names[idx]
    
    for (j in 1:length(idx)) { ## loop through variables with big enough loadings
      corr <- cor(x[,idx[j]],y,method = "spearman")
      sign <- sign(corr)
      
      if (condition == "auc"){
        AUC <- pROC::auc(y[,1], x[, idx[j]], levels = sort(unique(y[, 1])), direction = "<")
        AUCs <- c(AUCs, AUC)
      }
      
      corrs <- c(corrs, corr)
      signs <- c(signs, sign)
    }
    color <- recode(signs, "-1" = "Blue", "0" = "White", "1"= "Red")
    if (condition == "auc"){
      
      df <- data.frame(names, A_loading, AUCs, corrs, color)
      df <- df[order(-df$A_loading), ]
      top <- df[1:num_top_feats, ]
      
      df <- df[order(-df$AUCs), ]
      bot <- df[1:num_top_feats, ]
    }else if(condition == "corr"){
      
      df <- data.frame(names, A_loading, corrs, color)
      df <- df[order(-df$A_loading), ]
      top <- df[1:num_top_feats, ]
      
      df <- df[order(-abs(df$corrs)), ]
      bot <- df[1:num_top_feats, ]
    }
    final <- unique(rbind(top, bot))
    temp[[i]] <- final
    write.table(final, file = paste0(out_path, "/gene_list_",
                                     colnames(A)[i], ".txt"), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
  }
  names(temp) <- colnames(A)
  SLIDE_res$feature_res <- temp
  return(SLIDE_res)
}
