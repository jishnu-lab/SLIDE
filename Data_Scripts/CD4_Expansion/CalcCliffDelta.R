library(effsize)

# function to calculate Cliff's Delta value for the Z-factor matrix from multiple methods.
CalcCliffDelta <- function(sig_z, y, comb, sig_idx = NULL){
  #sig_z, the z matrix with significant LFs subsetted
  #y, the response vector
  #comb, a list of combinations for pairwise comparison (need to match the levels in y)
  all_cds <- data.frame()

  for (i in 1:dim(sig_z)[[2]]){
    col <- sig_z[ , i]
    df1 <- as.data.frame(col)
    df1['stages'] <- y[, 1]
    if (is.null(sig_idx) == TRUE){
      p <- ggplot(df1, aes(x=factor(stages), y=col)) + geom_boxplot(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(i)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    } else{
      p <- ggplot(df1, aes(x=factor(stages), y=col)) + geom_boxplot(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(sig_idx[[i]]))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    }
    print(p)
    #comb <- list(c(0, 1), c(1, 2), c(0, 2))
    if (is.null(sig_idx) == TRUE){
      idx <- rep(i, length(comb))
    } else {
      idx <- rep(sig_idx[[i]], length(comb))
    }
    b1 <- c()
    b2 <- c()
    deltas <- c()
    for (j in 1:length(comb)){
      s <- comb[[j]]
      b1 <- append(b1, s[1])
      b2 <- append(b2, s[2])
      res = cliff.delta(df1[df1['stages']==s[1], ]$col, df1[df1['stages']==s[2], ]$col, return.dm=TRUE)
      deltas <- append(deltas, res$estimate)
    }
    tmp_df <- do.call(rbind.data.frame, Map("c", idx, b1, b2, deltas))
    colnames(tmp_df) <- c("LF", "b1", "b2", "deltas")
    all_cds <- rbind(all_cds, tmp_df)
  }
  return(all_cds)
}


CalcNullCliffDelta <- function(z, comb, sig_idx){
  rand_z = z[, -c(sig_idx)]
  rand_idx = sample(ncol(rand_z), size = length(sig_idx), replace = FALSE)
  rand_z = rand_z[, rand_idx]
  cd <- CalcCliffDelta(rand_z, y, comb, sig_idx = rand_idx)
  return(cd)
}


CalcNullPVal <- function(real_cd, null_cd){
  real_cd <- real_cd[(real_cd["b1"] != 1 ), ]
  null_cd <- null_cd[(null_cd["b1"] != 1 ), ]
  res <- wilcox.test(abs(real_cd$deltas), abs(null_cd$deltas))
  p_val <- res$p.value
  return(p_val)
}

CliffDeltaPermute <- function(z, nrep, SLIDE_cd, mofa_cd, comb, sig_idx){
  SLIDE_pvals = c()
  mofa_pvals = c()
  null_cds = data.frame(matrix(ncol = 4))
  for (rep in 1:nrep){
    print(rep)
    null_cd <- CalcNullCliffDelta(z, comb = comb, sig_idx = sig_idx)
    SLIDE_pval = CalcNullPVal(SLIDE_cd, null_cd)
    mofa_pval = CalcNullPVal(mofa_cd, null_cd)
    SLIDE_pvals = append(SLIDE_pvals, SLIDE_pval)
    mofa_pvals = append(mofa_pvals, mofa_pval)
    colnames(null_cds) = colnames(null_cd)
    null_cds = rbind(null_cds, null_cd)
  }
  null_cds = null_cds[-1, ]
  return(list(SLIDE_pvals, mofa_pvals, null_cds))
}





