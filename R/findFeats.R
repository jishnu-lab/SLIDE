#' Find the features of interest in a given latent factor.
#'
#' For a given latent factor of interest, extract the features with the highest loadings in the
#' \eqn{A} matrix calculated by Essential Regression, find the Spearman correlation between the response (\eqn{y}) 
#' and these features, and find the p-values on the coefficients of these features from
#' a simple linear regression of \eqn{y} on these features (univariate). 
#' 
#' @param x a matrix or data frame of numeric values; the features
#' @param y a vector or data frame with one column of numeric values; the true response
#' @param z an integer; the column index of the latent factor of interest
#' @param er_res a list; the results from running Essential Regression
#' @param thresh a numeric value; the threshold to use to determine significance for the absolute loadings of each feature in \eqn{A}
#' @param p_thresh a numeric value; the threshold to use to determine significance for the univariate p-values
#' @return a list including the pure variables, mixed variables, features with significant p-values, features with 
#' significantly large loadings, correlations between response and features, and the univariate p-values for the 
#' latent factor of interest
#' @export

findFeats <- function(x, y, z, er_res, thresh = 0.01, p_thresh = 0.1) {
  ## scale the data
  x_std <- scale(x, T, T)
  y_std <- scale(y, T, T)
  ## extract the feature names
  feat_names <- colnames(x)
  ## "parse" the ER results to get cluster memberships
  clusters <- EssReg::readER(x = x, er_res = er_res)
  ## extract the A matrix
  a <- er_res$A

  ## find which features have absolute loadings in A that are greater than or
  ## equal to the specified cutoff
  top_feats <- which(abs(a[, z]) >= thresh)
  ## extract the loadings
  loadings <- a[top_feats, z]
  ## translate the indices of these significant features to their names
  top_feats <- indName(top_feats, feat_names, F)
  ## rename the loadings
  names(loadings) <- top_feats
  
  ## initialize results vectors
  corrs <- NULL
  p_vals <- NULL
  ## loop through features with significantly large loadings
  for (i in 1:length(top_feats)) { 
    ## get correlation between y and feature
    corr <- cor(y_std, x_std[, top_feats[i]], method = "spearman") 
    ## rename correlations
    names(corr) <- top_feats[i]
    ## append to results vector
    corrs <- c(corrs, corr)
    
    ## get p-value in linear model between one feature and y
    p_val <- summary(lm(y_std ~ x_std[, top_feats[i]]))$coefficients[, 4][2] 
    ## rename p-value
    names(p_val) <- top_feats[i]
    ## append to results vectof
    p_vals <- c(p_vals, p_val)
  }

  ## extract pure features
  pure_vars_inds <- which(abs(loadings) == 1)
  ## filter to just pure features in Z of interest
  pure_vars_in_clust <- loadings[pure_vars_inds] 
  ## filter mixed variables to just those in Z of interest
  mix_vars_in_clust <- loadings[-pure_vars_inds]
  ## find features with significant p-values
  sig_inds <- which(p_vals < p_thresh)
  ## filter p-values to just those that are significant
  sig_vars_in_clust <- top_feats[sig_inds]
  ## get just mixed features with significant p-values
  sig_mix_vars_in_clust <- setdiff(sig_vars_in_clust, pure_vars_in_clust)

  return (list("pure_vars" = pure_vars_in_clust,
               "mixed_vars" = mix_vars_in_clust,
               "pvals<p_thresh" = sig_mix_vars_in_clust,
               "loadings>=thresh" = loadings,
               "correlations_y.feat" = corrs,
               "pvals" = p_vals))
}
