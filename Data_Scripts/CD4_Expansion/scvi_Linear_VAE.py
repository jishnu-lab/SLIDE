############################ Loading packages ##################################
import scvi
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
############################Reading Annot Data##################################

adata = sc.read("/ix/djishnu/Javad/CDexpansion/raw/CD4.h5ad")
adata.layers["counts"] = adata.X.copy()
scvi.model.LinearSCVI.setup_anndata(adata, layer="counts")
model = scvi.model.LinearSCVI(adata, n_latent=4)
model.train(max_epochs=250, plan_kwargs={"lr": 5e-3}, check_val_every_n_epoch=10)

############################Saving the factors and loadings#####################

## Factors
Z_hat = model.get_latent_representation()
np.savetxt("/ix/djishnu/Javad/CDexpansion/scvi_factors.csv",Z_hat)

## Loadings
loadings = model.get_loadings()
loadings.to_csv("/ix/djishnu/Javad/CDexpansion/scvi_loadings.csv")
################################################################################
