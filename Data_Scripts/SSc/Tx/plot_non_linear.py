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
#%%
def GetPolyFit(gene, x, y, space):
    poly_est = np.polyfit(x[gene].to_numpy(), y["MRSS"].to_numpy(), 3)
    y_poly_pred = np.polyval(poly_est, np.linspace(space[0], space[1], space[2]))
    return y_poly_pred

def Plot(gene, x, y, y_poly_pred, space):    
    plt.plot(x[gene].to_numpy(), y["MRSS"].to_numpy(), "bo", color = "black")
    plt.plot(np.linspace(space[0], space[1], space[2]), y_poly_pred, "--")
    plt.xlabel("Gene Expression")
    plt.ylabel("MRSS")
    plt.title(gene)
    plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/" + gene+".pdf", 
                bbox_inches = "tight")
    plt.show()

    
    
x = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/Var50_mtrp.csv",
                index_col=0)

y = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/SkinScore_MRSS.csv",
                index_col=0)
    
#%%
gene = "c.3.WIF1"
assert gene in list(x.columns), print("Gene name error...")

#space = [0, 2.5, 50]
#y_poly_pred = GetPolyFit(gene, x, y, space)
#Plot(gene, x, y, y_poly_pred, space)

plt.plot(x[gene].to_numpy(), y["MRSS"].to_numpy(), "bo", color = "black")
plt.xlabel("Gene Expression")
plt.ylabel("MRSS")
plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/" + gene+".pdf", 
            bbox_inches = "tight")
plt.show()

#%%
# calculate R2
gene = "c.16.SAT1"
assert gene in list(x.columns), print("Gene name error...")

coef, p = spearmanr(x[gene].to_numpy(), y["MRSS"].to_numpy())
print(coef)
print(p)


#%%

z = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/z_matrix.csv",
                index_col=0)

y = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/SkinScore_MRSS.csv",
                index_col=0)

#sig_idx =[6, 10, 12, 47, 56, 85, 99]
sig_idx =[77, 90]
str_idx = []
for i in range(len(sig_idx)):
    t = "Z" + str(sig_idx[i])
    str_idx.append(t)
    
sig_z = z.loc[:, str_idx]


i = 0
plt.plot(sig_z.iloc[:, i].to_numpy(), y["MRSS"].to_numpy(), "bo", color = "black")
plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/" + str_idx[i]+".pdf", 
            bbox_inches = "tight")
plt.show()

#%% 
# Calc R2


z = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/z_matrix.csv",
                index_col=0)

y = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/SkinScore_MRSS.csv",
                index_col=0)

sig_idx =[6, 10, 12, 47, 56, 85, 99, 77, 90]
#sig_idx =[77, 90]
str_idx = []
for i in range(len(sig_idx)):
    t = "Z" + str(sig_idx[i])
    str_idx.append(t)
    
sig_z = z.loc[:, str_idx]

i = 8
print(sig_idx[i])
coef, p = spearmanr(sig_z.iloc[:, i].to_numpy(), y["MRSS"].to_numpy())
print(coef)
print(p)

