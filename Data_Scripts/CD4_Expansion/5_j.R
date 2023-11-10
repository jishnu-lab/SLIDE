library(effsize)
library(ggplot2)
source("CalcCliffDelta.R")

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


#################################################################
##                    Calculate Cliff Delta                    ##
#################################################################

comb <- list(c(0, 1), c(1, 2), c(0, 2))
sig_mofa_idx <- c(1, 2, 3, 7)
sig_mofa <- mofa[ ,sig_mofa_idx]
MOFA_cd <- CalcCliffDelta(sig_mofa, y, comb, sig_idx = sig_mofa_idx)


sig_idx <- c(19, 28, 9) # chosen by SLIDE
sig_z <- z[, sig_idx]
#sig_z <- scVI
SLIDE_cd <- CalcCliffDelta(sig_z, y, comb, sig_idx = sig_idx)


scVI_cd <- CalcCliffDelta(scVI, y, comb, sig_idx = NULL)



#################################################################
##                     Permutation Testing                     ##
#################################################################

permute_pvals <- CliffDeltaPermute(z, 30, SLIDE_cd, MOFA_cd, comb, sig_idx = sig_idx) # need to modify the CalcNullPDelta function
median(permute_pvals[[1]])
median(permute_pvals[[2]])

all_ps = as.data.frame(permute_pvals)
colnames(all_ps) = c("SLIDE", "MOFA+")
write.csv(all_ps, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/p_vals.csv")

permute_pvals <- CliffDeltaPermute(z, 30, SLIDE_cd, scVI_cd, comb, sig_idx = sig_z_idx) # need to modify the CalcNullPDelta function
median(permute_pvals[[1]])
median(permute_pvals[[2]])

all_ps = as.data.frame(permute_pvals)
colnames(all_ps) = c("SLIDE", "scVI")
write.csv(all_ps, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/p_vals_scVI.csv")


##################################################################
##                         Scatter Plot                         ##
##################################################################
# randomly choose 3 out of 4 LF for MOFA, picking 1 2 3
# choose 3 most significant LF for scVI, picking 1 3 4
Null_cd <- permute_pvals[[3]]
plot_df <- SLIDE_cd
plot_df['method'] <- rep("SLIDE", nrow(SLIDE_cd))
MOFA_cd['method'] <- rep("mofa", nrow(MOFA_cd))
Null_cd['method'] <- rep("null", nrow(Null_cd))
plot_df <- rbind(plot_df, MOFA_cd[1:9, ])
scVI_cd['method'] <- rep("scVI", nrow(MOFA_cd))
#plot_df <- rbind(plot_df, scVI_cd[scVI_cd["LF"] != 2, ], Null_cd)
plot_df <- rbind(plot_df, scVI_cd[scVI_cd["LF"] != 2, ])
#plot_df['method_idx'] <- c(rep(2, 9), rep(4, 9), rep(6, 9), rep(8, nrow(Null_cd)))
plot_df['method_idx'] <- c(rep(2, 9), rep(4, 9), rep(6, 9))

# removing 1 vs 2
toy_df <- plot_df[(plot_df["b1"] != 1 ), ]
write.csv(toy_df, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/5j.csv")

# perform mann whitney u test
#wilcox.test(abs(toy_df[toy_df["method"] == "SLIDE", ]$deltas), abs(toy_df[toy_df["method"] == "mofa", ]$deltas))


p <- ggplot(toy_df, aes(x = method_idx, y = abs(deltas))) + geom_point(aes(color=method)) + xlim(1, 7) + scale_color_manual(values=c('#cd9701', '#fa9c1b', '#c77cff', '#356920')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#median for SLIDE
SLIDE_med <- median(c(abs(plot_df[plot_df["method_idx"] == 2, ]['deltas']))$deltas)
mofa_med <- median(c(abs(plot_df[plot_df["method_idx"] == 4, ]['deltas']))$deltas)
vae_med <- median(c(abs(plot_df[plot_df["method_idx"] == 6, ]['deltas']))$deltas)


