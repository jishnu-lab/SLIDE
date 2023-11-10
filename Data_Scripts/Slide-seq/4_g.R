setwd("/Users/xiaoh/Desktop/Research/Hierarchical_ER/Data_Analyses/")
source("CalcCliffDelta.R")
source("GetPairwiseComb.R")
source("GetTopBenchmarkLF.R")

library(ggplot2)


##################################################################
##                 Import All Loading Matricies                 ##
##################################################################
z <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/SlideSeq/ER/050423/SLIDE_Results/z_matrix.csv",
                        row.names = 1))

y <- as.matrix(read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/SlideSeq/ER/050423/Data/y.csv",
                        row.names = 1))

mofa <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/SlideSeq/MOFA_VAE/MOFA_facros_v2.rds")

#################################################################
##                    Calculate Cliff Delta                    ##
#################################################################
comb <- GetPairwiseComb(y)

mofa_idx <- GetTopBenchmarkLF(mofa, y)
mofa_idx <- mofa_idx[1:6]# the 6 comes from the # of LF from SLIDE
sig_mofa <- mofa[, mofa_idx]
mofa_cd <- CalcCliffDelta(sig_mofa, y, comb, sig_idx = mofa_idx) #n_col: the LF number, b1 and b2: pairwise comparison bars, #deltas: CLiff Delta Value

sig_z_idx <- c(21, 24, 25, 28, 29, 7)
sig_z <- z[, sig_z_idx]
SLIDE_cd <- CalcCliffDelta(sig_z, y, comb, sig_idx = sig_z_idx)

#################################################################
##                     Permutation Testing                     ##
#################################################################

permute_pvals <- CliffDeltaPermute(z, 30, SLIDE_cd, mofa_cd, comb, sig_idx = sig_z_idx)
median(permute_pvals[[1]])
median(permute_pvals[[2]])

all_ps = as.data.frame(permute_pvals)
colnames(all_ps) = c("SLIDE", "MOFA+")
#write.csv(all_ps, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/SlideSeq/ER/050423/SLIDE_Results/p_vals.csv")

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
#plot_df['method_idx'] <- c(rep(2, nrow(SLIDE_cd)), rep(4, nrow(mofa_cd)), rep(6, nrow(Null_cd)))
plot_df['method_idx'] <- c(rep(2, nrow(SLIDE_cd)), rep(4, nrow(mofa_cd)))

toy_df <- plot_df[(plot_df["b1"] != 1 ), ]
write.csv(toy_df, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/SlideSeq/ER/050423/4_g.csv")

#result <- wilcox.test(abs(toy_df[toy_df["method"] == "SLIDE", ]$deltas), abs(toy_df[toy_df["method"] == "mofa", ]$deltas))


p <- ggplot(toy_df, aes(x = method_idx, y = abs(deltas))) + geom_point(aes(color=method)) + xlim(1, 7) + scale_color_manual(values=c('#cd9701', '#c77cff')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#median for SLIDE
SLIDE_med <- median(c(abs(toy_df[toy_df["method_idx"] == 2, ]['deltas']))$deltas)
mofa_med <- median(c(abs(toy_df[toy_df["method_idx"] == 4, ]['deltas']))$deltas)

###
# subset Z matrix to exclude the SLIDE LFs
# for n reps:
#   randomly picking size-matched LFs from the SLIDE Z matrix
#   null = calculate Cliff's Delta for the randomly chosen LFs
#   p1 = Mannwhitney(SLIDE_cd, null) 01 and 02
#   p2 = Mannwhitney(mofa_cd, null)
#   append p1 to p1_list
#   append p2 to p2_list
#
# final_p1 = median(p1_list)
# final_p2 = median(p2_list)


###
# for each comparison in SLIDE:
#  for nreps:
#   p_hat = randomly chose one LF to perform either 01 or 02 Cliff delta calculation
#   appent p_hat to p_hat_list
#
# mannwhitney(int, p_hat_lit) int is one scalar from the original CD test

#################################################################
##                     Permutation Testing                     ##
#################################################################

permute_pvals <- CliffDeltaPermute(z, 30, SLIDE_cd, mofa_cd, comb, sig_idx = sig_z_idx)
median(permute_pvals[[1]])
median(permute_pvals[[2]])

all_ps = as.data.frame(permute_pvals)
colnames(all_ps) = c("SLIDE", "MOFA+")
write.csv(all_ps, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/SlideSeq/ER/050423/SLIDE_Results/p_vals.csv")

