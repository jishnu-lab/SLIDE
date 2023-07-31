## VAE_Merfish
library(keras)
library(ggplot2)
library(SLIDE)
library(dplyr)
library(keras)

input_dim  <-241
latent_dim <- 13








x <- read.table("/ix/djishnu/Javad/SLIDEbench_merfish/x.csv",
                row.names = 1,
                sep=",",
                header = T)

y <- read.table("/ix/djishnu/Javad/SLIDEbench_merfish/y.csv",
                row.names = 1,
                sep = ",",
                header=T)

colnames(y) <-"y"




Z_ER <- read.table("/ix/djishnu/Javad/SLIDEbench_merfish/z_matrix.csv",header = T,
                   row.names = 1,sep = ",")


data_norm<- scale(x,T,T)



encoder <- keras_model_sequential() %>%
  layer_dense(units = latent_dim*2, activation = "sigmoid", input_shape = c(input_dim), kernel_initializer = "glorot_uniform")
#layer_dense(units = latent_dim*2)



decoder <- keras_model_sequential() %>%
  layer_dense(units = c(latent_dim), activation = "relu", input_shape = c(latent_dim)) %>%
  layer_dense(units = 128, activation = "relu") %>%
  layer_dense(units = input_dim, activation = "sigmoid")
## no sampling
vae <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output[, 1:latent_dim]))
# sampling




vae %>% compile(optimizer = "adam", loss = "CosineSimilarity")

# Train the VAE

vae_fit <- vae %>% fit(data_norm, data_norm, epochs = 100, batch_size = 24, validation_split = 0.2)
encoder_output <- stats::predict(encoder, data_norm)

sumVAE <- summary(lm(as.matrix(y)~encoder_output[,1:latent_dim]))
vAERsqaured <- sumVAE$r.squared

loading <- encoder %>% keras::get_weights()
vaeoadings <- encoder %>% keras::get_weights()
vaeoadings <- vaeoadings[[1]][,1:latent_dim]

row.names(vaeoadings) <- colnames(x)
colnames(vaeoadings) <- paste0("Z",1:latent_dim)








saveRDS(vaeoadings,"/ix/djishnu/Javad/SLIDEbench_merfish/MOFA_VAE/VAE.rds")





################################################################################
#Mofa
################################################################################

library(MOFA2)
data_norm <- list()
data_norm$data$view_1 <-  t(scale(x,T,T))
MOFAmodel <- create_mofa(data_norm$data)
MOFAmodel <- prepare_mofa(MOFAmodel)
MOFAmodel <- run_mofa(MOFAmodel, use_basilisk = TRUE)
#CD4_MOFA <- load_model(file = "CD4-Mofa.hdf5",remove_inactive_factors = FALSE)
factors <- get_factors(MOFAmodel, factors = "all")
factors <- factors$group1
MofaLoading <- MOFA2::get_weights(MOFAmodel)$view_1

saveRDS(factors,"/ix/djishnu/Javad/SLIDEbench_merfish/MOFA_VAE/MOFA_facros_v2.rds")
saveRDS(MofaLoading,"/ix/djishnu/Javad/SLIDEbench_merfish/MOFA_VAE/MOFA_loading_v2.rds")
