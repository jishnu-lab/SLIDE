
#mofa <- as.data.frame(readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/benchmark/MOFA-VAE/MOFA_facros.rds"))
#vae <- as.data.frame(readRDS("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/benchmark/MOFA-VAE/VAE.rds"))


#l1 <- mofa[order(mofa$Factor5, decreasing=TRUE),, drop = FALSE]
#l1 <- vae[order(vae$Z7, decreasing=TRUE),, drop = FALSE]
y <- as.matrix(utils::read.csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Poholek_Spatial/Spatial_ER_083022/Data/y.csv", row.names = 1))

p_vals <- c()
for (i in 1:ncol(mofa)){
  reg <- lm(y ~ mofa[, i])
  y_hat <- reg$fitted.values
  p_val <- cor.test(y, y_hat, method = "spearman")$p.value
  p_vals <- append(p_vals, p_val)
}


