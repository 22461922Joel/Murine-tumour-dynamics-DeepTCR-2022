'''
This script fits a model to distinguish between TCR-seq samples from the AB1 vs RENCA tumors.
'''
from DeepTCR.DeepTCR import DeepTCR_WF
import pandas as pd
import numpy as np
import pickle

df = pd.read_csv('../data/data.tsv',sep='\t')
df['Sample_ID'] = 'Model'+'_'+ df['Model'] +\
                  '_'+ 'Mouse' +'_'+df['Mouse'].astype(str) + \
                  '_' + 'Timepoint'+'_'+df['Timepoint'].astype(str)

timepoint = 0
df_sel  = df[df['Timepoint']==timepoint]

beta_sequences = np.array(df_sel['CDR3'])
v_beta = np.array(df_sel['v_gene'])
j_beta = np.array(df_sel['j_gene'])
counts = np.array(df_sel['count'])
sample_labels = np.array(df_sel['Sample_ID'])
class_labels = np.array(df_sel['Model'])
time_point = np.array(df_sel['Timepoint'])

DTCR = DeepTCR_WF('model_'+str(timepoint))
DTCR.Load_Data(beta_sequences=beta_sequences,v_beta=v_beta,j_beta=j_beta,
               counts=counts,sample_labels=sample_labels,class_labels=class_labels)

folds = 25
LOO = 6
epochs_min = 10
size_of_net = 'small'
num_concepts=64
hinge_loss_t = 0.3
train_loss_min=0.1
seeds = np.array(range(folds))
graph_seed = 0
subsample = None

DTCR.Monte_Carlo_CrossVal(folds=folds,LOO=LOO,
                          subsample=subsample,subsample_by_freq=True,subsample_valid_test=True,
                          epochs_min=epochs_min,
                          combine_train_valid=True,train_loss_min=train_loss_min,
                          hinge_loss_t=hinge_loss_t,graph_seed=graph_seed,size_of_net=size_of_net,
                          seeds=seeds,num_concepts=num_concepts)

with open('model_'+str(timepoint)+'_preds.pkl','wb') as f:
    pickle.dump([
        DTCR.beta_sequences,
        DTCR.v_beta,
        DTCR.j_beta,
        DTCR.counts,
        DTCR.freq,
        DTCR.sample_id,
        DTCR.class_id,
        time_point,
        DTCR.predicted
    ],f,protocol=4)

DTCR.DFs_pred['AB1'].to_csv('AB1_sample_preds_'+timepoint+'.csv',index=False)
DTCR.DFs_pred['RENCA'].to_csv('RENCA_sample_preds_'+timepoint+'.csv',index=False)

DTCR.AUC_Curve()
DTCR.Representative_Sequences(top_seq=50,make_seq_logos=False)

class_sel = 'RENCA'
DTCR.Residue_Sensitivity_Logo(beta_sequences=np.array(DTCR.Rep_Seq[class_sel]['beta'].iloc[0:25]),
                              v_beta = np.array(DTCR.Rep_Seq[class_sel]['v_beta'].iloc[0:25]),
                              j_beta = np.array(DTCR.Rep_Seq[class_sel]['j_beta'].iloc[0:25]),
                              class_sel=class_sel,
                              models=['model_'+str(x) for x in np.random.choice(25,10,replace=False)],
                              Load_Prev_Data=False,
                              min_size=0.5)