import sys
tissue = sys.argv[1]
sys.path.append('/data/yosef2/users/chenling/scVI/')
from scvi.dataset import AnnDatasetFromAnnData
sys.path.append('/data/yosef2/users/chenling/tabula-sapiens/')
from annotation.utils import *
from copy import deepcopy
import scanpy as sc
from anndata import read_h5ad
from sklearn.neighbors import NearestNeighbors

data_path = '/data/yosef2/users/chenling/TabulaSapiens'

rawdata_path = '%s/AnnotationsRound1/data/'%data_path
capital_tissue = tissue[0].upper() + tissue[1:]


if tissue =='pancreas':
    rawdata_path = '../../AnnotationsRound1/data/'
    tenxendo = read_h5ad(rawdata_path + 'tabula-sapiens-10X-pilot-filtered-%s.h5ad' % 'endopancreas')
    ss2endo = read_h5ad(rawdata_path + 'tabula-sapiens-facs-pilot-filtered-%s.h5ad' % 'Endopancreas')
    tenxexo = read_h5ad(rawdata_path + 'tabula-sapiens-10X-pilot-filtered-%s.h5ad' % 'exopancreas')
    tenxexo1 = read_h5ad(rawdata_path + 'tabula-sapiens-10X-pilot-filtered-%s.h5ad' % 'exopancreas1')
    ss2exo = read_h5ad(rawdata_path + 'tabula-sapiens-facs-pilot-filtered-%s.h5ad' % 'Exopancreas')
    raw = tenxendo.concatenate(ss2endo, tenxexo, tenxexo1, ss2exo)
    raw.obs.index = [x[:-2] for x in raw.obs.index]
else:
    tenx = read_h5ad(rawdata_path + 'tabula-sapiens-10X-pilot-filtered-%s.h5ad' % tissue)
    ss2 = read_h5ad(rawdata_path + 'tabula-sapiens-facs-pilot-filtered-%s.h5ad' % capital_tissue)
    raw = tenx.concatenate(ss2)
    raw.obs.index = [x[:-2] for x in raw.obs.index]


markers = {'Epithelial': {'CDH1',  'CLDN4',  'EPCAM',  'GCG'},
            'Endothelial': {'CA4', 'CDH5', 'CLDN5', 'LYVE1', 'PECAM1', 'VWF'},
            'Stromal': {'BGN', 'DCN', 'COL1A2'},
            'Immune': {'LCP1', 'PTPRC', 'RAC2'},
            'Pancreas': {'CELA3A', 'CPA1'}}

adata = deepcopy(raw)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)



for x in markers.keys():
    # sc.tl.score_genes(adata, ctrl_size=100, gene_list=markers[x],score_name = x)
    adata.obs[x] = np.sum(adata[:,list(markers[x])].X, axis=1)
    indicator = x+'_label'
    adata.obs[indicator] = (adata.obs[x].values> (2*len(markers[x])))
    print(x, np.sum(adata.obs[indicator]))


data = AnnDatasetFromAnnData(raw, 'batch')
model_file = '%s/scVImodels/vae.%s.%s.ann2.pkl' % (data_path, tissue, 'combined')
posterior = get_scvi_posterior(data, model_file)
latent, _, _ = posterior.get_latent()


all_indicators = pd.concat([adata.obs[x+'_label'] for x in markers.keys()], axis = 1)
unique_celltypes = ((all_indicators == True).sum(axis=1) <=1 )

labels = np.repeat(0, data.X.shape[0])
for i, c in enumerate(markers.keys()):
    idx = np.where(
        np.logical_and(
        unique_celltypes,
        adata.obs[c+'_label'].values==True))[0]
    labels[idx] = i+1

data.cell_types = ['unassigned']+list(markers.keys())
data.labels = labels.reshape(len(labels),1)
data.n_labels = len(data.cell_types)

# generate prediction for all cells at once
full, comp_pred = scanvi_pred(data, model_file, alternate=True)
comp_pred = np.asarray(comp_pred)

nbrs = NearestNeighbors(n_neighbors=10).fit(latent)
distances, indices = nbrs.kneighbors(latent)

smooth_comp = []
for i, x in enumerate(comp_pred):
    res = clustercomp(comp_pred[indices[i][1:]])
    smooth_comp.append(res)

meta = pd.DataFrame(index = adata.obs.index)
meta['Compartment_Prediction'] = comp_pred
meta['Compartment_Smoothed'] = smooth_comp
print( np.unique(smooth_comp, return_counts=True))
meta.to_csv(data_path+'/meta/%s.compartment.csv'%tissue)