############################ Loading packages ##################################
import scvi
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pickle
############################Reading Annot Data##################################

adata = sc.read("/ix/djishnu/Javad/CDexpansion/raw/CD4.h5ad")
scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata, n_latent=4)

# train the model
model.train()
latent = model.get_latent_representation()


with open("/ix/djishnu/Javad/CDexpansion/SCVI/scvi_factors.pickle", "wb") as file:
    pickle.dump(latent, file)
