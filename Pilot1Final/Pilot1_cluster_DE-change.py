#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import os
from anndata import read_h5ad
from sklearn.neighbors import NearestNeighbors


import warnings
warnings.filterwarnings('ignore')
import sys

sys.path.append('/data/yosef2/users/chenling/scVI/')
sys.path.append('/data/yosef2/users/chenling/tabula-sapiens/')
from annotation.utils import *
import scvi
print(scvi.__version__)

import os
os.getcwd()

import logging
import os
import pickle
from MulticoreTSNE import MulticoreTSNE as TSNE

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from hyperopt import hp

from scvi.inference import UnsupervisedTrainer, AlternateSemiSupervisedTrainer, SemiSupervisedTrainer
from scvi.models import VAE, SCANVI
from copy import deepcopy

from scvi.dataset.anndataset import AnnDatasetFromAnnData
import scanpy as sc
from anndata import read_h5ad
from anndata import AnnData

import matplotlib

n_epochs = 150
show_plot = True
test_mode = False
use_batches = True
use_cuda = True
lr = 1e-3
retrain=True


# In[2]:


rawdata_path = '../../TabulaSapiens/AnnotationsRound1/data/'
tenx = read_h5ad(rawdata_path + 'tabula-sapiens-10X-pilot-filtered.h5ad')
ss2 = read_h5ad(rawdata_path + 'tabula-sapiens-facs-pilot-filtered.h5ad')
raw = tenx.concatenate(ss2)
raw.obs.index = [x[:-2] for x in raw.obs.index]


# In[ ]:


ann = []
data_path = '/data/yosef2/users/chenling/TabulaSapiens/AnnotationFinal/'
tissue = 'blood'
all_files = os.listdir(data_path)
files = []
for c in ['immune','epithelial','stromal','endothelial','mixed']:
    files += [x for x in all_files if (('h5ad' in x and 'TSP' in x) and c in x and 'all' not in x)]

   
for f in files:
    data = read_h5ad(data_path + f)
    ann.append(data.obs)


# In[ ]:


ann = pd.concat(ann, axis=0, sort=False)


# In[ ]:


ann


# In[ ]:


temp = pd.concat([raw.obs[['batch']], ann], axis=1, sort=False)
assert np.sum(raw.obs.index == temp.index)==len(temp.index)
raw.obs = temp
temp.loc[temp['tissue_method'].astype(str)=='nan','tissue_method'] = 'blood_10x'
raw.obs.replace(np.nan, 'nan', regex=True, inplace=True)


# In[ ]:

_, raw.obs['batch_id'] = np.unique(temp['tissue_method'].values.astype(str),return_inverse=True)

# In[ ]:


all_dataset = AnnDatasetFromAnnData(raw, 'batch_id')


# In[ ]:


all_dataset.n_batches


# In[ ]:


os.getcwd()


# In[ ]:


model_file = '../scVImodels/AnnotationFinal.scVI.pkl'
posterior = get_scvi_posterior(all_dataset, model_file)
latent, _, _ = posterior.get_latent()


# In[ ]:


# Regenerating Compartment Prediction

compartment_pred = pd.read_csv(data_path + 'compartment_seeds_prediction.csv',index_col=0)
raw.obs = pd.concat([compartment_pred['predicted'],raw.obs], axis=1)
raw.obs.rename_axis({'predicted':'Compartment Prediction'}, inplace=True, axis=1)

pancreas_pred = pd.read_csv(data_path + 'pancreas_seeds_prediction.csv',index_col=0)
raw.obs = pd.concat([pancreas_pred['predicted'], raw.obs], axis=1)
raw.obs.rename_axis({'predicted':'Pancreas Prediction'}, inplace=True, axis=1)


# In[ ]:


# using scVI latent space in scanpy

raw.obsm["X_scvi"] = latent
sc.pp.neighbors(raw, n_neighbors=20, n_pcs=30, use_rep="X_scvi")
sc.tl.umap(raw)

X = raw.obsm['X_scvi']
nbrs = NearestNeighbors(n_neighbors=30).fit(X)
distances, indices = nbrs.kneighbors(X)

smooth_comp = []
compartments = deepcopy(raw.obs['Compartment Prediction'].values.astype(str))
for i,x in enumerate(compartments):
    res = clustercomp(compartments[indices[i][1:]])
    smooth_comp.append(res)
    

raw.obs['Smoothed Compartment Prediction'] = smooth_comp


# # Generating DE results

# In[ ]:


temp = np.unique(ann.loc[ann['tissue']=='blood','leiden']).astype(int)
len(temp) == np.max(temp)+1


# In[ ]:


np.unique(raw.obs['tissue'])


# In[ ]:


np.unique(raw.obs['Smoothed Compartment Prediction'])


# In[ ]:


temp = [x+'_'+y for x,y in zip(raw.obs['tissue'], raw.obs['Smoothed Compartment Prediction'])]


# In[ ]:


raw.obs['organ_compartment'] = temp


# In[ ]:


def DEbyOrganCompartment(raw, full, label, organs, compartments,  split_subset=None):
    pred_label = np.asarray(raw.obs[label].values)
    pred_celltype, pred_label = np.unique(pred_label, return_inverse=True)
    for organ in np.unique(raw.obs[organs]):
        if split_subset is None: 
            split_subset = np.unique(raw.obs.loc[raw.obs[organs].values == organ, compartments])
        writer = pd.ExcelWriter(organ+'.DE.change.xlsx', engine='xlsxwriter')
        res = []
        clust = []
        for compartment in split_subset:
            print(compartment)
            subset = np.logical_and(raw.obs[compartments].values == compartment,
                                   raw.obs[organs].values == organ)
            de_res, de_cluster = full.one_vs_all_degenes(cell_labels=pred_label,
                                                         subset=subset, mode='change')
            filtered = []
            for i, ct in enumerate(pred_celltype[de_cluster]):
                x = de_res[i]
#                 filt = np.logical_and(np.abs(x['bayes_factor'].values) > 1.3, x['raw_mean1'].values > 0.1)
                filt = x['proba_de']>0.90
                filtered.append(x.loc[filt])
#                 res = x.loc[filt]
#                 if len(ct)>31:
#                     ct = ct[:31]
#                 res.to_excel(writer, sheet_name=ct)
#                 res['celltype'] = ct
#             writer.save()
            temp = []
            for i,x in zip(de_cluster, filtered):
#                 logfc =np.log(x['scale1']/x['scale2']+1)
                logfc = x['median']
                clust_id = np.repeat(i, len(logfc))
                y = pd.DataFrame([clust_id, logfc.values, logfc.index], index=['id','logFC','Genenames']).T
                temp.append(y)
            if len(temp)>0:
                res = pd.concat(temp, axis=1)
                res.to_excel(writer, sheet_name=compartment)
        writer.save()


# In[ ]:


DEbyOrganCompartment(raw=raw, full=posterior, label='leiden', organs='tissue',
               compartments ='Smoothed Compartment Prediction')

