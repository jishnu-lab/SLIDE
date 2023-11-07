z <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/z_matrix.csv",
                        row.names = 1))

y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Data/y.csv",
                        row.names = 1))

sig_idx <- c(21, 3, 1, 17, 14, 8, 16)
sig_z <- z[ ,sig_idx]
print(colnames(sig_z))


# box plot, no abs(), no gridline/ no background. rename x and ys.
for (i in 1:dim(sig_z)[[2]]){
  col <- sig_z[ , i]
  df1 <- as.data.frame(col)
  df1['stages'] <- y[, 1]
  #ggplot(df1, aes(x=factor(stages), y=col)) + geom_violin(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(sig_idx[[i]]))
  print(ggplot(df1, aes(x=factor(stages), y=col)) + geom_boxplot(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(sig_idx[[i]])) + ylim(0, 2.5))
}


# perform mann whitney u test
idx_single <- which(y == 0)
idx_small <- which(y == 1)
idx_med <- which(y == 2)
for (i in 1:dim(sig_z)[[2]]){
  print(sig_idx[[i]])
  df2 = NULL
  cat("This is for cluster " , sig_idx[[i]], "...")
  col <- sig_z[ ,i][idx_single]
  df2 <- as.data.frame(col)
  df2["group"] <- c(rep("single", length(col)))
  col <-sig_z[ ,i][idx_small]
  tmp <- as.data.frame(col)
  tmp["group"] <- c(rep("small", length(col)))
  df2 <- rbind(df2, tmp)
  col <-sig_z[ ,i][idx_med]
  tmp <- as.data.frame(col)
  tmp["group"] <- c(rep("medium", length(col)))
  df2 <- rbind(df2, tmp)

  kruskal.test(col ~ factor(group), data = df2)
}

