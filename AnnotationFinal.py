import sys
tissue = sys.argv[1]

sys.path.append('/data/yosef2/users/chenling/scVI/')
from scvi.dataset import AnnDatasetFromAnnData
sys.path.append('/data/yosef2/users/chenling/tabula-sapiens/')
from annotation.utils import *
from copy import deepcopy
import scanpy as sc
from anndata import read_h5ad

data_path = '/data/yosef2/users/chenling/TabulaSapiens'

files = os.listdir('%s/CompartmentSplit/'%data_path)
h5ad = ['%s/CompartmentSplit/'%data_path + x for x in files if (tissue in x)]
compartments = [x.split('.')[1] for x in h5ad]
adata = [read_h5ad(x) for x in h5ad]
adata = adata[0].concatenate(adata[1:])
adata.obs.index = [x[:-2] for x in adata.obs.index]

rawdata_path = '%s/AnnotationsRound1/data/'%data_path
capital_tissue = tissue[0].upper() + tissue[1:]

tenx = read_h5ad(rawdata_path + 'tabula-sapiens-10X-pilot-filtered-%s.h5ad' % tissue)
ss2 = read_h5ad(rawdata_path + 'tabula-sapiens-facs-pilot-filtered-%s.h5ad' % capital_tissue)
raw = tenx.concatenate(ss2)

temp = pd.concat([raw.obs[['10X_plate', 'tissue', 'batch', 'n_genes']],
                  adata.obs[['G2M_score', 'S_score', 'input_ann2',
                             'leiden_scvi_split', 'scanvi_ann', 'smooth_comp']]],
                 axis=1, sort=False)
temp['input_ann2'].fillna('unassigned', inplace=True)
temp['smooth_comp'].fillna('nan', inplace=True)
assert np.sum(raw.obs.index == temp.index) == len(temp.index)
raw.obs = temp
data = AnnDatasetFromAnnData(raw, batch_label='batch')
assert len(np.unique(data.batch_indices.ravel())) == 2

model_file = '%s/scVImodels/vae.%s.%s.ann2.pkl' % (data_path, tissue, 'combined')
posterior = get_scvi_posterior(data, model_file)
latent, _, _ = posterior.get_latent()

# using scVI latent space in scanpy
raw.obsm["X_scvi"] = latent
sc.pp.neighbors(raw, n_neighbors=10, n_pcs=30, use_rep="X_scvi")
sc.tl.umap(raw)
sc.tl.leiden(raw, key_added="leiden_scvi", resolution=1)

# use input_ann for labels
data.cell_types, data.labels = np.unique(raw.obs['input_ann2'], return_inverse=True)
data.labels = data.labels.reshape(len(data.labels), 1)
data.n_labels = len(data.cell_types)

# generate prediction for all cells at once
scanvi_model_file = '%s/scVImodels/scanvi.%s.%s.ann2.pkl' % (data_path, tissue, 'combined')
full, pred_celltype = scanvi_pred(data, model_file, scanvi_model_file)
raw.obs['pred2'] = pred_celltype

scanvi_model_file = '%s/scVImodels/scanvi.%s.%s.alternate.ann2.pkl' % (data_path, tissue, 'combined')
full, pred_celltype = scanvi_pred(data, model_file, scanvi_model_file, alternate=True)
raw.obs['pred2_alternate'] = pred_celltype

if tissue is 'muscle':
    compartments = ['Stromal', 'Muscle', 'Endothelial', 'Immune']
else:
    compartments = ['Stromal', 'Epithelial', 'Endothelial', 'Immune']

# training scanVI for each compartment individually
raw.obs['pred_by_compartment'] = 'nan'
for compartment in compartments:
    data_subset = deepcopy(data)
    data_subset.update_cells(raw.obs['smooth_comp'].values == compartment)
    labels = [data.cell_types[i] for i in data_subset.labels.ravel()]
    data_subset.cell_types, labels = np.unique(labels, return_inverse=True)
    data_subset.labels = labels.reshape(len(labels), 1)
    scanvi_model_file = '%s/scVImodels/scanvi.%s.%s.ann2.pkl' % (data_path, tissue, compartment)
    full_subset, pred_subset = scanvi_pred(data_subset, model_file, scanvi_model_file)
    raw.obs.loc[raw.obs['smooth_comp'] == compartment, 'pred_by_compartment'] = np.asarray(pred_subset)


# training scanVI for each compartment individually
raw.obs['pred_by_compartment_alternate'] = 'nan'


for compartment in compartments:
    data_subset = deepcopy(data)
    data_subset.update_cells(raw.obs['smooth_comp'].values == compartment)
    labels = [data.cell_types[i] for i in data_subset.labels.ravel()]
    data_subset.cell_types, labels = np.unique(labels, return_inverse=True)
    data_subset.labels = labels.reshape(len(labels), 1)
    scanvi_model_file = '%s/scVImodels/scanvi.%s.%s.alternate.ann2.pkl' % (data_path, tissue, compartment)
    full_subset, pred_subset = scanvi_pred(data_subset, model_file, scanvi_model_file, alternate=True)
    raw.obs.loc[raw.obs['smooth_comp'] == compartment, 'pred_by_compartment_alternate'] = np.asarray(pred_subset)
    raw.obs['pred_by_compartment_alternate'].fillna('nan', inplace=True)

miss1 = check_missing_celltypes(raw, 'pred2')
miss2 = check_missing_celltypes(raw, 'pred_by_compartment')
miss3 = check_missing_celltypes(raw, 'pred2_alternate')
raw.obs['pred_by_compartment_alternate'].fillna('nan', inplace=True)
miss4 = check_missing_celltypes(raw, 'pred_by_compartment_alternate')
print(miss1, "\n", miss2, "\n", miss3, "\n", miss4)


# saving objects for raw counts, imputed counts and scanpy normalized counts
raw.write_h5ad('CombineCompartments.scANVI.%s.h5ad' % tissue)
imputed = posterior.sequential().imputation(transform_batch=[0, 1])
new_adata = deepcopy(raw)
new_adata.X = imputed
new_adata.write_h5ad('CombineCompartments.scANVI.%s.imputed.h5ad' % tissue)
new_adata = deepcopy(raw)
new_adata.X = imputed
sc.pp.normalize_total(new_adata, target_sum=1e4)
sc.pp.log1p(new_adata)
sc.pp.scale(adata, max_value=10)
new_adata.write_h5ad('CombineCompartments.scANVI.%s.scanpy_norm.h5ad' % tissue)

# Differential Expression
DEbyCompartment(raw, posterior, 'pred_by_compartment', 'smooth_comp', split_subset=compartments,
                filename=data_path + '/de_res/%s.%s.1.xlsx' % (tissue, 'combined'))

DEbyCompartment(raw, posterior, 'pred2', 'smooth_comp', split_subset=compartments,
                filename=data_path + '/de_res/%s.%s.2.xlsx' % (tissue, 'combined'))

DEbyCompartment(raw, posterior, 'pred_by_compartment_alternate', 'smooth_comp', split_subset=compartments,
                filename=data_path + '/de_res/%s.%s.3.xlsx' % (tissue, 'combined'))

DEbyCompartment(raw, posterior, 'pred2_alternate', 'smooth_comp', split_subset=compartments,
                filename=data_path + '/de_res/%s.%s.4.xlsx' % (tissue, 'combined'))
