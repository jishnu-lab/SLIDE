library(doParallel)
library(dplyr)

files <- list.files("/Users/xiaoh/Desktop/Research/Hierarchical_ER/R")
source_files <- paste0("/Users/xiaoh/Desktop/Research/Hierarchical_ER/R/", files)
sapply(source_files, source)

final_res <- NULL
x <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Data/sig_match.csv", row.names = 1))
z <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Results/results_w/z_matrix.csv", row.names = 1)) ## standardized
y <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Data/Y.csv", row.names = 1))
er_res <- readRDS('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Results/results_w/final_er_output.rds')

print(er_res$K)
SLIDE_res <- SLIDE(z, y, method = 4, do_interacts = TRUE, betas = NULL, top_prop = NULL, marginals = NULL,
                   spec = 0.3, fdr = 0.1, niter = 1000, elbow = FALSE, f_size = 100, parallel = TRUE, ncore = 10)

print(SLIDE_res)
SLIDE_param <- c(4, 0.3, 0.1, 1000, 100)
names(SLIDE_param) <- c("method", "spec", "fdr", "niter", "f_size")


ks <- c(137, 151, 166)
temp <- NULL
for (i in 1:length(ks)){
  feature_res <-SigGenes_A(er_res, ks[[i]], as.matrix(x), matrix(y,ncol = 1), thre=0.05, negcol="blue", posCol="red",orderbycol=T)
  feature_res <- feature_res[order(feature_res[, 1], decreasing = TRUE), ]
  temp[[i]] <- feature_res
}

names(temp) <- ks
final_res$SLIDE_res <- SLIDE_res
final_res$feature_res <- temp
final_res$SLIDE_param <- SLIDE_param
final_res$feature_thresh <- 0.05

#saveRDS(final_res, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Results/results_w/SLIDE_Res.rds")

#Res <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_091322/Result/w_scale/SLIDE_Res.rds")
Res <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Joglekar_CD4/Expansion/HER_091822/Results/results_w/SLIDE_res.rds")

for ( i in length(Res$feature_res)){
  print("new list")
  df = Res$feature_res[[i]][order(abs(Res$feature_res[[i]]$correlation), decreasing = TRUE), ]
  print(df)
}

