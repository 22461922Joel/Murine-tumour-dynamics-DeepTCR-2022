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

color_dict = {'RS':'g','NR':'r'}
#continuous
model_sel = 'AB1' #choose which model to plot signatures
sig_sel = 'AB1' #specify the signature one wants to quantify

df_sel = df[df['Model']==model_sel]
df_sel['w_'+sig_sel] = df_sel[sig_sel]*df_sel['freq']
df_agg = df_sel.groupby(['Sample_ID','Mouse','Timepoint','Response']).agg({'w_'+sig_sel:'sum'}).reset_index()
df_agg.sort_values(by=['Timepoint','Response'],inplace=True)

#plot signature over time points
plt.figure()
sns.boxplot(data=df_agg,x='Timepoint',y='w_'+sig_sel)
plt.ylabel(sig_sel,fontsize=16)
plt.xlabel('Timepoint',fontsize=16)
plt.tight_layout()

#plot signature over time points stratified by response
plt.figure()
sns.boxplot(data=df_agg,x='Timepoint',y='w_'+sig_sel,hue='Response',palette=color_dict)
plt.ylabel(sig_sel,fontsize=16)
plt.xlabel('Timepoint',fontsize=16)
plt.tight_layout()