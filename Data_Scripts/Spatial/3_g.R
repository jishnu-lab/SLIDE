library(ggplot2)
z <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/z_matrix.csv",
                        row.names = 1))

y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Data/y.csv",
                        row.names = 1))

sig_idx <- c(21, 3, 1, 17, 14, 8, 16)
sig_z <- z[ ,sig_idx]
print(colnames(sig_z))


# box plot, no abs(), no gridline/ no background. rename x and ys.
for (i in 1:dim(sig_z)[[2]]){
  #flip the sign for Z1 and Z17 because the positively correlated gene (with Y) has negative loading in A matrix
  if ((colnames(sig_z)[i] == "Z1") || (colnames(sig_z)[i] == "Z17")){
    col <- sig_z[ , i] * -1
  } else {
    col <- sig_z[ , i]
  }
  df1 <- as.data.frame(col)
  df1['stages'] <- y[, 1]
  #print(ggplot(df1, aes(x=factor(stages), y=col)) + geom_violin(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(sig_idx[[i]])))
  print(ggplot(df1, aes(x=factor(stages), y=col)) + geom_boxplot(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(sig_idx[[i]]))+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")))
  # write.csv(df1, paste0("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/3g/",
  #                       colnames(sig_z)[i], ".csv" ))
}



# perform mann whitney u test
idx_single <- which(y == 0)
idx_small <- which(y == 1)
idx_med <- which(y == 2)
for (i in 1:dim(sig_z)[[2]]){
  df2 = NULL
  cat("This is for cluster " , sig_idx[[i]], "...\n")
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

  print(kruskal.test(col ~ factor(group), data = df2))
}

