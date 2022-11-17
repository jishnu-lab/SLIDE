#library(EssReg)

x_path <- "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Data/sig_match.csv"
x <- as.matrix(utils::read.csv(x_path, row.names = 1))
x <- scale(x, T, T)
er_res <- readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Results/results_w/final_er_output.rds")

z_matrix <- predZ(x, er_res)
colnames(z_matrix) <- paste0("Z", c(1:ncol(z_matrix)))
write.csv(z_matrix, "/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/Pre_Post_Treatment/HER_091022/Results/results_w/z_matrix.csv", row.names = TRUE, col.names = TRUE)
