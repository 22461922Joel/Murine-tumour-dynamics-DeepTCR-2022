'''
This script takes a model fit to distinguish AB1 from RENCA TCRs on timepoint 0 and assesses the predictive signature at timepoints 2,4,6.
'''
from DeepTCR.DeepTCR import DeepTCR_WF
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

df = pd.read_csv('../data/data.tsv',sep='\t')
df['Sample_ID'] = 'Model'+'_'+ df['Model'] +\
                  '_'+ 'Mouse' +'_'+df['Mouse'].astype(str) + \
                  '_' + 'Timepoint'+'_'+df['Timepoint'].astype(str)

df_count = df.groupby(['Sample_ID']).agg({'count':'sum'}).reset_index()
count_dict = dict(zip(df_count['Sample_ID'],df_count['count']))
df['count_sum'] = df['Sample_ID'].map(count_dict)
df['count_sum'] = df['count_sum'].astype(int)
df['freq'] = df['count']/df['count_sum']

timepoint_train = 0

DFs = []
for timepoint_inf in [2,4,6]:
    df_sel = df[df['Timepoint'] == timepoint_inf]
    beta_sequences = np.array(df_sel['CDR3'])
    v_beta = np.array(df_sel['v_gene'])
    j_beta = np.array(df_sel['j_gene'])
    counts = np.array(df_sel['count'])
    sample_labels = np.array(df_sel['Sample_ID'])
    class_labels = np.array(df_sel['Model'])

    DTCR = DeepTCR_WF('model_'+str(timepoint_train))
    DTCR.Sample_Inference(sample_labels=sample_labels,
                          beta_sequences=beta_sequences,
                          v_beta=v_beta,
                          j_beta=j_beta,
                          counts=counts,
                          models=['model_' + str(x) for x in np.random.choice(25, 10, replace=False)])

    df_label = df_sel[['Sample_ID','Model']].drop_duplicates()
    label_dict = dict(zip(df_label['Sample_ID'],df_label['Model']))
    df_pred = DTCR.Inference_Pred_Dict['AB1'][:]
    df_pred['Model'] = df_pred['Samples'].map(label_dict)
    df_pred['Model_bin'] = 0
    df_pred['Model_bin'][df_pred['Model']=='AB1']=1
    DFs.append(df_pred)

fig,ax = plt.subplots()
for df_pred,timepoint_inf in zip(DFs,[2,4,6]):
    roc_score = roc_auc_score(df_pred['Model_bin'],df_pred['Pred'])
    fpr,tpr,_ = roc_curve(df_pred['Model_bin'],df_pred['Pred'])
    plt.plot(fpr,tpr,lw=2,label='%s (area = %0.4f)' % (str(timepoint_inf), roc_score))
plt.xlabel('False Positive Rate',fontsize=16)
plt.ylabel('True Positive Rate',fontsize=16)
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.legend(loc='lower right')
plt.tight_layout()

DFs = pd.concat(DFs)
DFs.to_csv('Supp4B.csv',index=False)
