# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Lafyatis/All_Cell_Type/HER_041422/data/Var50df.csv', index_col = 0)
exp = data.sum(axis = 0)


hist = plt.hist(exp, bins = 20, edgecolor = 'black')
#plt.axvline(x = 5, color = 'r', linestyle = '-')
plt.xlabel("Count")
plt.ylabel("# Genes")
plt.show()  


filtered = exp[exp > 5].index
len(exp[exp <= 5])

data = data[data.columns.intersection(list(filtered))]

data.to_csv('/Users/xiaoh/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Research_Files/Lafyatis/All_Cell_Type/HER_041422/data/Var50_thresh5.csv')
