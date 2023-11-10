import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu 

#%% functions
def GetSigZ(sig_idx, z):
    #subset z matrix to significant factors
    str_idx = []
    for i in range(len(sig_idx)):
        t = "Z" + str(sig_idx[i])
        str_idx.append(t)
    sig_z = z.loc[:, str_idx]
    return sig_z

def CalcRho(sig_z, y):
    stats = dict()
    for i in range(sig_z.shape[1]):
        coef, p = spearmanr(sig_z.iloc[:, i].to_numpy(), y["MRSS"].to_numpy())
        #print(i)
        #print(coef)
        #print(p)
        stats[i] = [coef, p]
    return stats

def MannWhitney(SLIDE, MOFA):
    coef, p = mannwhitneyu(SLIDE, MOFA)
    #print(coef)
    #print(p)
    return p


#%% calculate effect size
z = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/z_matrix.csv",
                index_col=0)
mofa_z = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/benchmark/MOFA_VAE/Mofa_factors.csv",
                index_col=0)
vae_z = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/benchmark/MOFA_VAE/VAE_factors.csv")

y = pd.read_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Data/SkinScore_MRSS.csv",
                index_col=0)

sig_idx =[6, 10, 12, 47, 56, 85, 99, 77, 90]
sig_z = GetSigZ(sig_idx, z)
SLIDE_stats = CalcRho(sig_z, y)
MOFA_stats  = CalcRho(mofa_z, y)
VAE_stats = CalcRho(vae_z, y)

#plot a dot plot from the stats
SLIDE_stats = pd.DataFrame.from_dict(SLIDE_stats, orient='index', columns = ["coef", "p"])
SLIDE_stats['method_idx'] = np.repeat(2, len(SLIDE_stats))

MOFA_stats = pd.DataFrame.from_dict(MOFA_stats, orient='index', columns = ["coef", "p"])
MOFA_stats['method_idx'] = np.repeat(4, len(MOFA_stats))

VAE_stats = pd.DataFrame.from_dict(VAE_stats, orient='index', columns = ["coef", "p"])
VAE_stats['method_idx'] = np.repeat(6, len(VAE_stats))

plot_df = pd.concat([SLIDE_stats, MOFA_stats, VAE_stats])
plot_df.to_csv("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/2_k.csv")

# plot a dot plot with x as method_idx and y as coef with x from 1 to 5
plt.figure(figsize=(5, 7))
plt.plot(SLIDE_stats['method_idx'], SLIDE_stats['coef'].abs(), "bo", color = "#c77cff")
plt.plot([1.5, 2.5], [SLIDE_stats['coef'].abs().median(), SLIDE_stats['coef'].abs().median()], color='black', linestyle='-', linewidth=2)
plt.plot(MOFA_stats['method_idx'], MOFA_stats['coef'].abs(), "bo", color = "#cd9701")
plt.plot([3.5, 4.5], [MOFA_stats['coef'].abs().median(), MOFA_stats['coef'].abs().median()], color='black', linestyle='-', linewidth=2)
plt.plot(VAE_stats['method_idx'], VAE_stats['coef'].abs(), "bo", color = "#f761cc")
plt.plot([5.5, 6.5], [VAE_stats['coef'].abs().median(), VAE_stats['coef'].abs().median()], color='black', linestyle='-', linewidth=2)
plt.xlim(1, 7)
# plt.savefig("/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/SLIDE/Lafyatis_SSc/All_Cell_Type/HER_081222/Results/rho_SLIDE_MOFA_VAE_V2.pdf", 
#             bbox_inches = "tight")
plt.show()


#%%

nreps = 50
slide_ps = []
mofa_ps = []
vae_ps = []

for rep in range(nreps):
    print(rep)
    # get random Z and the rho for that Z
    rand_z = z.drop(z.columns[sig_idx], axis = 1)
    rand_z = rand_z.sample(n=len(sig_idx), axis=1).reset_index(drop=True)
    print(rand_z.columns)
    rand_stats = CalcRho(rand_z, y)
    rand_stats = pd.DataFrame.from_dict(rand_stats, orient='index', columns = ["coef", "p"])
    # get the mannwhitneyu p value from SLIDE, MOFA and random
    #rand_stats.fillna(0, inplace=True)
    slide_p = MannWhitney(abs(pd.DataFrame(SLIDE_stats).iloc[:, 0]), abs(pd.DataFrame(rand_stats).iloc[:, 0]))
    mofa_p = MannWhitney(abs(pd.DataFrame(MOFA_stats).iloc[:, 0]), abs(pd.DataFrame(rand_stats).iloc[:, 0]))
    vae_p = MannWhitney(abs(pd.DataFrame(VAE_stats).iloc[:, 0]), abs(pd.DataFrame(rand_stats).iloc[:, 0]))
    slide_ps.append(slide_p)
    mofa_ps.append(mofa_p)
    vae_ps.append(vae_p)

np.median(slide_ps)
np.median(mofa_ps)


