library(doParallel)
library(dplyr)

files <- list.files("/Users/xiaoh/Desktop/Research/Hierarchical_ER/R")
source_files <- paste0("/Users/xiaoh/Desktop/Research/Hierarchical_ER/R/", files)
sapply(source_files, source)

final_res <- NULL
x <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Data/ER_X.csv", row.names = 1))
z <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/z_matrix.csv", row.names = 1)) ## standardized
y <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Data/y.csv", row.names = 1))
er_res <- readRDS('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/final_er_output.rds')

SLIDE_res <- SLIDE(z, y, method = 4, do_interacts = TRUE, betas = NULL, top_prop = NULL, marginals = NULL,
                   spec = 0.3, fdr = 0.1, niter = 5000, elbow = FALSE, f_size = 21, parallel = TRUE, ncore = 10)

print(SLIDE_res)
#feature_res <- findFeats(x, y, z, er_res, thresh = 0.01, p_thresh = 0.1)
SLIDE_param <- c(4, 0.3, 0.1, 5000, 21)
names(SLIDE_param) <- c("method", "spec", "fdr", "niter", "f_size")


ks <- c(21, 3, 14, 1, 8, 16, 17)
temp <- NULL
#names <-
for (i in 1:length(ks)){
  feature_res <-SigGenes_A(er_res, ks[[i]], as.matrix(x), matrix(y,ncol = 1), thre=0.1, negcol="blue", posCol="red",orderbycol=T)
  feature_res <- feature_res[order(feature_res[, 1], decreasing = TRUE), ]
  temp[[i]] <- feature_res
}

names(temp) <- ks
final_res$SLIDE_res <- SLIDE_res
final_res$feature_res <- temp
final_res$SLIDE_param <- SLIDE_param
final_res$feature_thresh <- 0.1

saveRDS(final_res, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/SLIDE_Res.rds")
