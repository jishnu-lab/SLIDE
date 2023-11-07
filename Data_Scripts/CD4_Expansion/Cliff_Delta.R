library(effsize)


##################################################################
##                 Import All Loading Matricies                 ##
##################################################################
z <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/Results/results_w/z_matrix.csv",
                        row.names = 1))

y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/Data/expansion_y.csv",
                        row.names = 1))

mofa <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/Benchmark/MOFA-VAE/MOFA_facros_v2.rds")

scVI <- read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/Benchmark/SCVI/scvi_factors.csv",
                           sep = " ", row.names = 1)
row.names(scVI) <- scVI$cells
scVI$cells <- NULL

##################################################################
##                      Check LF Direction                      ##
##################################################################
# Check_LF_Direction <- function(z, y){
#   directions = c()
#   for (n in 1:ncol(z)){
#     col = z[, n]
#     sign = sign(cor(col, y, method = "spearman"))
#     directions = append(directions, sign)
#   }
#   return(directions)
# }
# 
# sig_mofa <- c(1, 2, 3, 7)
# sig_mofa <- mofa[ ,sig_mofa]
# mofa_dir <- Check_LF_Direction(sig_mofa, y)
# 
# sig_idx <- c(19, 28, 9) # chosen by SLIDE
# sig_z <- z[, sig_idx]
# SLIDE_dir <- Check_LF_Direction(sig_z, y)

#################################################################
##                    Calculate Cliff Delta                    ##
#################################################################


Calc_Cliff_Delta <- function(sig_z){
  all_cds <- data.frame()
  # box plot, no abs(), no gridline/ no background. rename x and ys.
  for (i in 1:dim(sig_z)[[2]]){
    col <- sig_z[ , i]
    df1 <- as.data.frame(col)
    df1['stages'] <- y[, 1]
    p <- ggplot(df1, aes(x=factor(stages), y=col)) + geom_boxplot(aes(fill = factor(stages))) + scale_fill_brewer(palette="Dark2") + ggtitle(as.character(i)) + ylim(-2, 2) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(p)
    comb <- list(c(0, 1), c(1, 2), c(0, 2))
    idx <- rep(i, length(comb))
    b1 <- c()
    b2 <- c()
    deltas <- c()
    for (j in 1:3){
      s <- comb[[j]]
      b1 <- append(b1, s[1])
      b2 <- append(b2, s[2])
      res = cliff.delta(df1[df1['stages']==s[1], ]$col, df1[df1['stages']==s[2], ]$col, return.dm=TRUE)
      deltas <- append(deltas, res$estimate)
    }
    tmp_df <- do.call(rbind.data.frame, Map("c", idx, b1, b2, deltas))
    colnames(tmp_df) <- c("n_col", "b1", "b2", "deltas")
    all_cds <- rbind(all_cds, tmp_df)
  }
  return(all_cds)
}


sig_mofa <- c(1, 2, 3, 7)
sig_mofa <- mofa[ ,sig_mofa]

sig_idx <- c(19, 28, 9) # chosen by SLIDE
sig_z <- z[, sig_idx]
#sig_z <- scVI
mofa_cd <- Calc_Cliff_Delta(sig_mofa)
scVI_cd <- Calc_Cliff_Delta(scVI)
SLIDE_cd <- Calc_Cliff_Delta(sig_z)

#mofa_cds <- all_cds

# perform mann whitney u test
idx_single <- which(y == 0)
idx_small <- which(y == 1)
idx_med <- which(y == 2)
for (i in 1:dim(sig_z)[[2]]){
  df2 = NULL
  cat("This is for column " , i, "...\n")
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


##################################################################
##                         Scatter Plot                         ##
##################################################################
# randomly choose 3 out of 4 LF for MOFA, picking 1 2 3
# choose 3 most significant LF for scVI, picking 1 3 4

plot_df <- SLIDE_cd
plot_df['method'] <- rep("SLIDE", nrow(SLIDE_cd))
mofa_cd['method'] <- rep("mofa", nrow(mofa_cd))
plot_df <- rbind(plot_df, mofa_cd[1:9, ])
scVI_cd['method'] <- rep("scVI", nrow(mofa_cd))
plot_df <- rbind(plot_df, scVI_cd[scVI_cd["n_col"] != 2, ])

#plot_df <- plot_df[(plot_df["b1"] != 0 ) &(plot_df["b1"] != 2 ), ]


p <- ggplot(plot_df, aes(x=factor(method), y=abs(deltas))) + geom_boxplot(aes(fill = factor(method))) + scale_fill_brewer(palette="Set2") + geom_jitter(color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



