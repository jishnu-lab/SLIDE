setwd("/Users/xiaoh/Desktop/Research/Hierarchical_ER/Data_Analyses/")
source("GetTopBenchmarkLF.R")
source("CalcCliffDelta.R")

library(effsize)
library(ggplot2)


##################################################################
##                 Import All Loading Matricies                 ##
##################################################################
z <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/z_matrix.csv",
                        row.names = 1))

y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Data/y.csv",
                        row.names = 1))

mofa <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/benchmark/MOFA-VAE/MOFA_facros.rds")

#################################################################
##                    Calculate Cliff Delta                    ##
#################################################################
comb <- list(c(1, 2), c(2, 3), c(1, 3))

# get the most significant LFs from MOFA+ to size match with SLIDE LF
sig_mofa_idx <- GetTopBenchmarkLF(mofa, y)
sig_mofa_idx <- sig_mofa_idx[1:7] # 7 is the size of LF from SLIDE 2, 1, 4, 3, 7, 6, 5
sig_mofa <- mofa[ ,sig_mofa_idx]


mofa_cd <- CalcCliffDelta(sig_mofa, y, comb, sig_idx = sig_mofa_idx) #n_col: the LF number, b1 and b2: pairwise comparison bars, #deltas: CLiff Delta Value

sig_z_idx <- c(14, 21, 1, 8, 16, 17, 3)
sig_z <- z[, sig_z_idx]
SLIDE_cd <- CalcCliffDelta(sig_z, y, comb, sig_idx = sig_z_idx)

#################################################################
##                     Permutation Testing                     ##
#################################################################

permute_pvals <- CliffDeltaPermute(z, 30, SLIDE_cd, mofa_cd, comb, sig_idx = sig_z_idx) # need to modify the CalcNullPDelta function
median(permute_pvals[[1]])
median(permute_pvals[[2]])

all_ps = as.data.frame(permute_pvals)
colnames(all_ps) = c("SLIDE", "MOFA+")
#write.csv(all_ps, "p_vals.csv")

##################################################################
##                         Scatter Plot                         ##
##################################################################
# randomly choose 3 out of 4 LF for MOFA, picking 1 2 3
# choose 3 most significant LF for scVI, picking 1 3 4
#Null_cd <- permute_pvals[[3]]
plot_df <- SLIDE_cd
plot_df['method'] <- rep("SLIDE", nrow(SLIDE_cd))
mofa_cd['method'] <- rep("mofa", nrow(mofa_cd))
#Null_cd['method'] <- rep("null", nrow(Null_cd))
#plot_df <- rbind(plot_df, mofa_cd, Null_cd)
plot_df <- rbind(plot_df, mofa_cd)
#plot_df['method_idx'] <- c(rep(2, 21), rep(4, 21), rep(6, nrow(Null_cd)))
plot_df['method_idx'] <- c(rep(2, 21), rep(4, 21))
toy_df <- plot_df[(plot_df["b1"] != 2 ), ]
write.csv(toy_df, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/no_scale/3_h.csv")
#result <- wilcox.test(abs(toy_df[toy_df["method"] == "SLIDE", ]$deltas), abs(toy_df[toy_df["method"] == "mofa", ]$deltas))

p <- ggplot(toy_df, aes(x = method_idx, y = abs(deltas))) + geom_point(aes(color=method)) + xlim(1, 7) + scale_color_manual(values=c('#cd9701', '#c77cff')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#median for SLIDE
SLIDE_med <- median(c(abs(toy_df[toy_df["method_idx"] == 2, ]['deltas']))$deltas)
mofa_med <- median(c(abs(toy_df[toy_df["method_idx"] == 4, ]['deltas']))$deltas)
