'''
This script takes models fit by model_classifier.py and response_classifier.py and applies them at a per sequence level to assign tumor and resposne specific signatures to all TCRs in all samples.
'''
from DeepTCR.DeepTCR import DeepTCR_WF
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
import os

if os.path.exists('tcrs_scored.csv'):
    df = pd.read_csv('tcrs_scored.csv')
else:
    df = pd.read_csv('../data/data.tsv',sep='\t')
    df['Sample_ID'] = 'Model'+'_'+ df['Model'] +\
                      '_'+ 'Mouse' +'_'+df['Mouse'].astype(str) + \
                      '_' + 'Timepoint'+'_'+df['Timepoint'].astype(str)

    df_count = df.groupby(['Sample_ID']).agg({'count':'sum'}).reset_index()
    count_dict = dict(zip(df_count['Sample_ID'],df_count['count']))
    df['count_sum'] = df['Sample_ID'].map(count_dict)
    df['count_sum'] = df['count_sum'].astype(int)
    df['freq'] = df['count']/df['count_sum']

beta_sequences = np.array(df['CDR3'])
v_beta = np.array(df['v_gene'])
j_beta = np.array(df['j_gene'])

DTCR = DeepTCR_WF('model')
out = DTCR.Sample_Inference(sample_labels=None,
                      beta_sequences=beta_sequences,
                      v_beta=v_beta,
                      j_beta=j_beta,
                      batch_size=10000,
                            models=['model_'+str(x) for x in np.random.choice(25,10,replace=False)])

df['AB1'] = out[:,0]
df['RENCA'] = out[:,1]

DTCR = DeepTCR_WF('response_AB1_0')
out = DTCR.Sample_Inference(sample_labels=None,
                      beta_sequences=beta_sequences,
                      v_beta=v_beta,
                      j_beta=j_beta,
                      batch_size=10000,
                            # models=None)
                            models=['model_'+str(x) for x in np.random.choice(25,10,replace=False)])

df['AB1_NR'] = out[:,0]
df['AB1_RS'] = out[:,1]

DTCR = DeepTCR_WF('response_RENCA_4')
out = DTCR.Sample_Inference(sample_labels=None,
                      beta_sequences=beta_sequences,
                      v_beta=v_beta,
                      j_beta=j_beta,
                      batch_size=10000,
                            # models=None)
                            models=['model_'+str(x) for x in np.random.choice(25,10,replace=False)])

df['RENCA_NR'] = out[:,0]
df['RENCA_RS'] = out[:,1]

df.to_csv('tcrs_scored.csv',index=False)


df_ab1 = df[df['Model']=='AB1']

idx = np.random.choice(len(df_ab1),50000,replace=False)
x,y,c,_,_ = GKDE(np.array(df_ab1['AB1'].iloc[idx]),
                 np.array(df_ab1['AB1_RS'].iloc[idx]),
                 np.array(df_ab1['freq'].iloc[idx]))
plt.scatter(x,y,c=c,cmap='jet')
plt.xlabel('AB1')
plt.ylabel('AB1_RS')

spearmanr(np.array(df['AB1']),np.array(df['AB1_RS']))

x_col,y_col = 'AB1','AB1_RS'
x_col,y_col = 'RENCA','RENCA_RS'
idx = np.random.choice(len(df),10000,replace=False)
x,y,c,_,_ = GKDE(np.array(df[x_col].iloc[idx]),
                 np.array(df[y_col].iloc[idx]))
                 #np.array(df['freq'].iloc[idx]))
plt.scatter(x,y,c=c,cmap='jet')
plt.xlabel(x_col)
plt.ylabel(y_col)
spearmanr(np.array(df[x_col]),np.array(df[y_col]))
