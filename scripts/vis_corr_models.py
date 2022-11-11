'''
This script creates plots looking at the correlations between various tumor and response-specific signatures at the TCR level
'''
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, spearmanr
def GKDE(x,y,z=None):
    xy = np.vstack([x, y])
    kernel = gaussian_kde(xy,weights=z)
    z = kernel(xy)
    r = np.argsort(z)
    x ,y, z = x[r], y[r], z[r]
    return x,y,z,kernel,r

df = pd.read_csv('tcrs_scored.csv')

x_col,y_col = 'AB1','AB1_RS'
# x_col,y_col = 'RENCA','RENCA_RS'
idx = np.random.choice(len(df),50000,replace=False)
x,y,c,_,_ = GKDE(np.array(df[x_col].iloc[idx]),
                 np.array(df[y_col].iloc[idx]))
                 #np.array(df['freq'].iloc[idx]))
plt.scatter(x,y,c=c,cmap='jet')
plt.xlabel(x_col,fontsize=16)
plt.ylabel(y_col,fontsize=16)
plt.tight_layout()
spearmanr(np.array(df[x_col]),np.array(df[y_col]))
