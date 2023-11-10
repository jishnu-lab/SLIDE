#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:17:35 2022

@author: xiaoh
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu 

#%%
# def GetPolyFit(gene, x, y, space):
#     poly_est = np.polyfit(x[gene].to_numpy(), y["MRSS"].to_numpy(), 3)
#     y_poly_pred = np.polyval(poly_est, np.linspace(space[0], space[1], space[2]))
#     return y_poly_pred

# def Plot(gene, x, y, y_poly_pred, space):    
#     plt.plot(x[gene].to_numpy(), y["MRSS"].to_numpy(), "bo", color = "black")
#     plt.plot(np.linspace(space[0], space[1], space[2]), y_poly_pred, "--")
#     plt.xlabel("Gene Expression")
#     plt.ylabel("MRSS")
#     plt.title(gene)
#     plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/" + gene+".pdf", 
#                 bbox_inches = "tight")
#     plt.show()

def main(genes, x, y):
    #plot the gene expressions with respect to MRSS
    exp = dict()
    for gene in genes:
        print(f"gene name: {gene}")
        assert gene in list(x.columns), print("Gene name error...")
        exp[gene] = x[gene].to_numpy()
        coef, p = spearmanr(x[gene].to_numpy(), y["MRSS"].to_numpy())
        print(f"spearman correlation coef:{coef}")
        print(f"spearman correlation p_val:{p}")

        plt.plot(x[gene].to_numpy(), y["MRSS"].to_numpy(), "bo", color = "black")
        plt.xlabel("Gene Expression")
        plt.ylabel("MRSS")
        # plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/" + gene+".pdf", 
        #             bbox_inches = "tight")
        plt.show()
    exp['MRSS'] = y["MRSS"].to_numpy()
    exp = pd.DataFrame.from_dict(exp)
    return exp
    
x = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/Var50_mtrp.csv",
                index_col=0)

y = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/SkinScore_MRSS.csv",
                index_col=0)

genes = ["c.3.WIF1", "c.6.CCL19", "c.3.MRPS6", "c.16.SAT1", "c.1.APOE", "c.7.S100A9",\
        "c.3.SAA1", "c.5.IGFBP5"]
ssc_g = main(genes, x, y)
ssc_g.to_csv("2_g_h.csv")


#%%
def GetSigZ(sig_idx, z):
    #subset z matrix to significant factors
    str_idx = []
    for i in range(len(sig_idx)):
        t = "Z" + str(sig_idx[i])
        str_idx.append(t)
    sig_z = z.loc[:, str_idx]
    return sig_z

z = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/z_matrix.csv",
                index_col=0)

# z = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/benchmark/MOFA_VAE/Mofa_factors.csv",
#                 index_col=0)


y = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/SkinScore_MRSS.csv",
                index_col=0)


sig_idx =[6, 10, 12, 47, 56, 85, 99, 77, 90]
sig_z = GetSigZ(sig_idx, z)

exp = dict()
for i in range(sig_z.shape[1]):
    print(f"plotting for LF {sig_z.columns[i]}")
    if (sig_z.columns[i] == "Z10") or (sig_z.columns[i] == "Z85"):
        exp[sig_z.columns[i]] = (sig_z.iloc[:, i].to_numpy()) * -1
    else:
        exp[sig_z.columns[i]] = (sig_z.iloc[:, i].to_numpy())
    plt.plot(exp[sig_z.columns[i]], y["MRSS"].to_numpy(), "bo", color = "black")
    plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/" +sig_z.columns[i]+".pdf", 
                bbox_inches = "tight")
    plt.show()
    coef, p = spearmanr(exp[sig_z.columns[i]], y["MRSS"].to_numpy())
    print(f"spearman correlation coef:{coef}")
    print(f"spearman correlation p_val:{p}")
exp['MRSS'] = y["MRSS"].to_numpy()
exp = pd.DataFrame.from_dict(exp)

exp.to_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/2_i/2_i.csv")


# i = 0
# # plt.plot(z.iloc[:, i].to_numpy(), y["MRSS"].to_numpy(), "bo", color = "black")
# # plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/" + str_idx[i]+".pdf", 
# #             bbox_inches = "tight")
# plt.show()

# coef, p = spearmanr(z.iloc[:, i].to_numpy(), y["MRSS"].to_numpy())
