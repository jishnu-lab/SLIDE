library(ggplot2)

z <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/Results/results_w/z_matrix.csv",
                        row.names = 1))

y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/Data/expansion_y.csv",
                        row.names = 1))

sig_idx <- c(19, 28, 9, 29)
sig_z <- z[ ,sig_idx]
print(colnames(sig_z))


# box plot, no abs(), no gridline/ no background. rename x and ys.
for (i in 1:dim(sig_z)[[2]]){
  if (colnames(sig_z)[[i]] %in% c("Z19", "Z9")){
    col <- sig_z[ , i] * -1
  }else{
    col <- sig_z[ , i]
  }
  df1 <- as.data.frame(col)
  #df1['col'] <- -1 * df1["col"]
  df1['col'] <- df1["col"]
  df1['stages'] <- y[, 1]
  #print(ggplot(df1, aes(x=factor(stages), y=col)) + geom_boxplot(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(sig_idx[[i]])))
  print(ggplot(df1, aes(x=factor(stages), y=col)) + geom_boxplot(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(sig_idx[[i]]))
        +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black")))
}



# toy = df1[df1['stages'] == 0, ]['col'][, 1]
# min(df1[df1['stages'] == 0, ]['col'][, 1])


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
  #wilcox.test(single_col, small_col, paired =  FALSE)
  #wilcox.test(small_col, med_col, paired = FALSE)
  #wilcox.test(single_col, med_col, paired = FALSE)
}
