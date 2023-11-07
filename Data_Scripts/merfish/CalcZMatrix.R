library(EssReg)

x_path <- "x.csv"
x <- as.matrix(utils::read.csv(x_path, row.names = 1))
x <- scale(x, T, T)
er_res <- readRDS("./results/final_er_output.rds")

z_matrix <- predZ(x, er_res)
colnames(z_matrix) <- paste0("Z", c(1:ncol(z_matrix)))
write.csv(z_matrix, "./results/z_matrix.csv", row.names = TRUE, col.names = TRUE)
