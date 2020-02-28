# peripheral neuron 'CL:0000111'
# basal cell 'CL:0000646'
# also only treating cells that all algorithms agree on as seed labels

import sys
sys.path.append('/data/yosef2/users/chenling/scVI/')
from scvi.dataset import AnnDatasetFromAnnData
sys.path.append('/data/yosef2/users/chenling/tabula-sapiens/')
from annotation.utils import *
from copy import deepcopy
import scanpy as sc
from anndata import read_h5ad

data_path = '/data/scratch/users/chenling/TabulaSapiens/AnnotationFinal/'
# change after moving things to scratch
model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/'
csv = sys.argv[1]
csv = model_path + 'chenling_refinement.csv'
csv = pd.read_table(csv, sep=',', skiprows=2)
adata = data_path + 'CombineCompartments.scANVI.all.count.h5ad'
adata = read_h5ad(adata)

# map tissue_tech
tissue_tech_list = ['Bladder_ss2', 'Blood_ss2', 'Lung_ss2', 'Muscle_ss2',
       'PancreasEndocrine_ss2', 'PancreasExocrine_ss2', 'bladder_10x',
       'blood_10x', 'endopancreas_10x', 'exopancreas1_10x',
       'exopancreas2_10x', 'lung_10x', 'muscle_10x']

batch_id = [tissue_tech_list.index(x) for x in adata.obs['Tissue Tech']]
adata.obs['batch'] = batch_id

data = AnnDatasetFromAnnData(adata, batch_label='batch')
model_file = model_path + 'scVImodels/AnnotationFinal.scVI.pkl'
posterior = get_scvi_posterior(data, model_file)
latent, _, _ = posterior.get_latent()

label_names = ['scANVI Prediction By Organ', 'scANVI Prediction By Organ-Compartment',
               'scANVI Prediction By Organ with Alternate Training',
               'scANVI Prediction By Organ-Compartment with Alternate Training']
pairwise = []
for i in range(len(label_names)):
    for j in np.arange(i+1, len(label_names)):
        pairwise.append(adata.obs[label_names[i]].values.astype(str)==adata.obs[label_names[j]].values.astype(str))

pairwise = np.asarray(pairwise)
count = np.sum(pairwise, axis=0)

labels = adata.obs['scANVI Prediction By Organ'].values.astype(str)
labels[count<6] = 'unassigned'
len(np.unique(adata.obs['Manual Annotation Round 2']))
# note: requiring all methods to agree did not miss additional cell types other than the ones that was missed already by 'scANVI Prediction by Organ'

# include manually annotated cell types that were not predicted in the last round
missed = []
for x in np.unique(adata.obs['Manual Annotation Round 2']):
    if x not in np.unique(labels):
        missed.append(x)
        labels[adata.obs['Manual Annotation Round 2']==x] = x


# include newly annotated cell types
assert adata.shape[0] == np.sum(adata.obs.index.values == csv['index'].values)
new_celltype = np.unique(csv[csv.columns[1]])
new_celltype = new_celltype[new_celltype!='unassigned']
for x in new_celltype:
    labels[csv[csv.columns[1]]==x] = x
#
#
# use input_ann for labels
data.cell_types, data.labels = np.unique(labels, return_inverse=True)
data.labels = data.labels.reshape(len(data.labels), 1)
data.n_labels = len(data.cell_types)
#
# generate prediction for all cells at once
full, pred_celltype = scanvi_pred(data, model_file, model_path+'scanvi/all.manual.semi.20.pkl',nlabels=20)
x, y = np.unique(pred_celltype, return_counts=True)
_, pred_celltype = scanvi_pred(data, model_file, nlabels=20)

missed = []
for c in np.unique(data.cell_types):
    if c not in np.unique(pred_celltype):
        print(c)
        missed.append(c)

adata.obs['new_pred'] = pred_celltype
adata.write('temp.h5ad')
# generate prediction for each organ separately
per_organ_pred = np.zeros(adata.shape[0])
per_organ_pred = per_organ_pred.astype(str)

for organ in np.unique(adata.obs['tissue']):
    adata_sub = adata[adata.obs['tissue'].values == organ]
    data_subset = AnnDatasetFromAnnData(adata_sub, batch_label='batch')
    model_file = '%sscVImodels/vae.%s.ann2.pkl' % (model_path, organ.lower())
    data_subset.cell_types, data_subset.labels = np.unique(labels[adata.obs['tissue'].values == organ], return_inverse=True)
    data_subset.labels = data_subset.labels.reshape(len(data_subset.labels), 1)
    data_subset.n_labels = len(data_subset.cell_types)
    posterior = get_scvi_posterior(data_subset, model_file)
    _, pred_celltype = scanvi_pred(data_subset, model_file, nlabels=20)
    per_organ_pred[adata.obs['tissue'] == organ] = pred_celltype


adata = read_h5ad('temp.h5ad')
adata.obs['new_pred_per_organ'] = per_organ_pred
adata.write('temp.h5ad')
# when using all cells from all tissues there are a lot more cell types that are missing from the prediction
# alternate misses more cell types
# 50 annotated per cell type v.s. 30 annotated per cell type
# what are they manually annotated as?

