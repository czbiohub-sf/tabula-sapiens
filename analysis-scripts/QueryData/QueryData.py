import sys
# query = sys.argv[1]
# tissue = sys.argv[2]
retrain = False
sys.path.append('/data/yosef2/users/chenling/scVI/')
from scvi.dataset import AnnDatasetFromAnnData
sys.path.append('/data/yosef2/users/chenling/tabula-sapiens/')
from utils.annotations_utils import *
from anndata import read_h5ad, read_mtx
from scanpy import read_10x_mtx

import obonet
obo = '../TabulaSapiens/data/cl.obo.txt'
f = open(obo, "r")
co = obonet.read_obo(f)
f.close()

ontology_dict = {}
for x in co.nodes.keys():
    ontology_dict[co.nodes[x]['name']] = x


# query = '../TabulaSapiensData/PublicDatasets/lung/lung.droplet.h5ad'
# dir = 'lung_atlas'

# query = '../TabulaSapiensData/PublicDatasets/lung/droplet_reyfman_normal_lung.20200319.RC4.h5ad' # 25246
# query_barcode = '../TabulaSapiensData/PublicDatasets/lung/reyfman.cellid.csv'
# dir = 'droplet_reyfman'

query = '../TabulaSapiensData/PublicDatasets/lung/dropseq_barga_normal_lung.20200319.RC4.h5ad' #5916 #6007
query_barcode = '../TabulaSapiensData/PublicDatasets/lung/barga.cellid.csv'
dir = 'dropseq_barga'

query = read_h5ad(query)

# check that the genes in the new dataset is the same as the old one for reyfman and barga dataset
if dir !='lung_atlas':
    if dir == 'dropseq_barga':
        temp = query.obs.index
    if dir == 'droplet_reyfman':
        temp = [x.split('_')[2] for x in query.obs.index]
    a,b = np.unique(temp, return_counts=True)
    np.sum(b>1) / len(b)
    query_barcode = pd.read_csv(query_barcode)
    filter = np.asarray([x in query.obs.index.values for x in query_barcode['x']])
    count = query.raw.X[filter]
    query.X = count.astype(int)
    query.obs = query.obs.loc[query_barcode[filter]['x']]
    query.obs.index == query_barcode['x'][filter]

csv = query.X
csv = pd.DataFrame(csv)
csv.columns = query.var.index.values
csv.index = query.obs.index.values
csv.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/count.csv'%dir)
labels = query.obs['free_annotation']
labels.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/labels.csv'%dir)

# load atlas data
data_path = '/data/scratch/users/chenling/TabulaSapiens/AnnotationFinal/'
# model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/droplet_reyfman/'
# res_path = 'droplet_reyfman/'
# model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/dropseq_barga/'
# res_path = 'dropseq_barga/'
model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/lung_atlas/'
res_path = 'lung_atlas/'
# model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/Tal/'
# res_path = 'Tal/'

if not os.path.exists(model_path):
    os.mkdir(model_path)
if not os.path.exists(model_path+'scanvi/'):
    os.mkdir(model_path+'scanvi/')
if not os.path.exists(res_path):
    os.mkdir(res_path)

adata = data_path + 'CombineCompartments.scANVI.all.count.h5ad'
adata = read_h5ad(adata)
model_file = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/AnnotationFinal.scVI.pkl'
lung_model_file = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/AnnotationFinal.lung.scVI.pkl'
lung_combined_model_file = model_path + 'AnnotationFinal.lung_combined.scVI.pkl'

# filter genes in query to be the same as the atlas data
def match_genes(query, adata):
    combined = adata.concatenate(query, join='outer')
    genes_combined = combined.var
    genes_combined['indices'] = np.arange(combined.shape[1])
    shared = pd.concat([adata.var, genes_combined], join='inner', axis=1)
    assert np.sum(shared.index==adata.var.index) == adata.shape[1]
    combined = combined[:, shared['indices'].values]
    return combined[combined.obs['batch']=='1',:]

matched_query = match_genes(query, adata)
query_data = AnnDatasetFromAnnData(matched_query)

###################################################################
# train model with all datasets
###################################################################
tissue_tech_list = ['Bladder_ss2', 'Blood_ss2', 'Lung_ss2', 'Muscle_ss2',
   'PancreasEndocrine_ss2', 'PancreasExocrine_ss2', 'bladder_10x',
   'blood_10x', 'endopancreas_10x', 'exopancreas1_10x',
   'exopancreas2_10x', 'lung_10x', 'muscle_10x']

batch_id = [tissue_tech_list.index(x) for x in adata.obs['Tissue Tech']]
adata.obs['batch'] = batch_id

train_data = AnnDatasetFromAnnData(adata, batch_label='batch')
labels = adata.obs['Manual Annotation Round 2']
posterior = get_scvi_posterior(train_data, model_file)

train_data.cell_types, train_data.labels = np.unique(labels, return_inverse=True)
train_data.labels = train_data.labels.reshape(len(train_data.labels), 1)
train_data.n_labels = len(train_data.cell_types)
full, pred_celltype = scanvi_pred(train_data, model_file,
                              model_path+'scanvi/all.manual.semi.20.pkl',
                              nlabels=30, retrain=True)
###################################################################
# train model with just the lung dataset
###################################################################
lung_data = adata[adata.obs['tissue']=='Lung']
csv = lung_data.X
csv = pd.DataFrame(csv.todense())
csv.columns = lung_data.var.index.values
csv.index = lung_data.obs.index.values
dir = 'pilot1'
csv.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/count.csv'%dir)
labels = query.obs['free_annotation']
labels.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/labels.csv'%dir)

tissue_tech_list = ['Lung_ss2','lung_10x']
batch_id = [tissue_tech_list.index(x) for x in lung_data.obs['Tissue Tech']]
lung_data.obs['batch'] = batch_id
train_data = AnnDatasetFromAnnData(lung_data, batch_label='batch')
labels = lung_data.obs['Manual Annotation Round 2']

posterior = get_scvi_posterior(train_data, lung_model_file)

train_data.cell_types, train_data.labels = np.unique(labels, return_inverse=True)
train_data.labels = train_data.labels.reshape(len(train_data.labels), 1)
train_data.n_labels = len(train_data.cell_types)
full, pred_celltype = scanvi_pred(train_data, lung_model_file,
                              model_path+'scanvi/lung.manual.semi.20.pkl',
                              nlabels=30, retrain=True)
###################################################################
# retrain entire model with just the lung and the new dataset
###################################################################
from anndata import AnnData
query = AnnData(X=query_data.X, obs=query.obs, var=pd.DataFrame(index=query_data.gene_names))
lung_combined = lung_data.concatenate(query)
lung_combined.obs['Tissue Tech'].fillna('Query_10X', inplace=True)
lung_combined.obs['Manual Annotation Round 2'].fillna('unassigned', inplace=True)
tissue_tech_list = list(np.unique(lung_combined.obs['Tissue Tech']))
batch_id = [tissue_tech_list.index(x) for x in lung_combined.obs['Tissue Tech']]
lung_combined.obs['batch'] = batch_id
train_data = AnnDatasetFromAnnData(lung_combined, batch_label='batch')
labels = lung_combined.obs['Manual Annotation Round 2']

posterior = get_scvi_posterior(train_data, lung_combined_model_file, retrain=False)

train_data.cell_types, train_data.labels = np.unique(labels, return_inverse=True)
train_data.labels = train_data.labels.reshape(len(train_data.labels), 1)
train_data.n_labels = len(train_data.cell_types)
full, pred_celltype = scanvi_pred(train_data, lung_combined_model_file,
                              model_path+'scanvi/lung_combined.manual.semi.20.pkl',
                              nlabels=30, retrain=True)


# manually adjust so the dimensions in the model matches
celltype1 = np.unique(adata.obs['Manual Annotation Round 2'])
celltype2 = np.unique(lung_data.obs['Manual Annotation Round 2'])
lung_combined.obs['Manual Annotation Round 2'].fillna('unassigned', inplace=True)
celltype3 = np.unique(lung_combined.obs['Manual Annotation Round 2'])

query_celltype_tal = [celltype3[i] for i in pred_celltype]

query_data.cell_types = np.unique(labels)

query_data.n_labels = 70
query_data.n_batches = 13
_, query_celltype1 = scanvi_pred(query_data, model_file,
                            model_path+'scanvi/all.manual.semi.20.pkl',
                            forward_only=True)

query_celltype1 = [celltype1[i] for i in query_celltype1]

query_data.n_labels = 44
query_data.n_batches = 2
_, query_celltype2 = scanvi_pred(query_data, lung_model_file,
                            model_path+'scanvi/lung.manual.semi.20.pkl',
                            forward_only=True)

query_celltype2 = [celltype2[i] for i in query_celltype2]

query_data.n_labels = 44
query_data.n_batches = 3
_, query_celltype3 = scanvi_pred(query_data, lung_combined_model_file,
                            model_path+'scanvi/lung_combined.manual.semi.20.pkl',
                            forward_only=True)

query_celltype3 = [celltype3[i] for i in query_celltype3]




np.save(file=res_path+'query_celltype1.npy', arr=query_celltype1)
np.save(file=res_path+'query_celltype2.npy', arr=query_celltype2)
np.save(file=res_path+'query_celltype3.npy', arr=query_celltype3)
query.obs.to_csv(res_path + 'obs.csv')

# # temp = pd.crosstab(labels, query_celltype)
#
# for idx, query_celltype in enumerate([query_celltype1, query_celltype2, query_celltype3]):
#     temp = confusion_matrix(labels, query_celltype)
#     names = np.sort(np.unique(np.concatenate([labels, query_celltype])))
#
#     res = []
#     for i,x in enumerate(names):
#         for j,y in enumerate(names):
#             if temp[i,j]!=0:
#                 res.append([x, y, temp[i,j]])
#
#     res = pd.DataFrame(res)
#     res = res[res[2]>10]
#     res = res.sort_values(2)
#     res.to_csv(res_path+'pretrained_model_performance.%i.csv'%idx)
#
# # find CL terms

import scanpy as sc

query_celltype1 = np.load(res_path + 'query_celltype1.npy')
query_celltype2 = np.load(res_path + 'query_celltype2.npy')
query_celltype3 = np.load(res_path + 'query_celltype3.npy')

query.obs['pred1'] = query_celltype1
query.obs['pred2'] = query_celltype2
query.obs['pred3'] = query_celltype3

train_data = AnnDatasetFromAnnData(query)
posterior = get_scvi_posterior(train_data, res_path+'scVI.pkl')
latent, _, _ = posterior.get_latent()

query.obsm["X_scvi"] = latent
sc.pp.neighbors(query, n_neighbors=20, n_pcs=30, use_rep="X_scvi")
sc.tl.umap(query)
query.write(res_path + 'predicted.h5ad')

lung_combined.obsm["X_scvi"] = latent
sc.pp.neighbors(lung_combined, n_neighbors=20, n_pcs=30, use_rep="X_scvi")
sc.tl.umap(lung_combined)
lung_combined.obs['pred'] = query_celltype_tal
lung_combined.obs = lung_combined.obs[['Manual Annotation Round 2','Tissue Tech', 'batch','n_genes', 'pred']]
qc = sc.pp.calculate_qc_metrics(lung_combined)[0]
temp = pd.concat([lung_combined.obs, qc], axis=1)
lung_combined.obs = temp
lung_combined.write(res_path + 'predicted.h5ad')
