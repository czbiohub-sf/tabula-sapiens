import sys
retrain = False
sys.path.append('/data/yosef2/users/chenling/scVI/')
from scvi.dataset import AnnDatasetFromAnnData
sys.path.append('/data/yosef2/users/chenling/tabula-sapiens/')
from utils.annotations_utils import *
from anndata import read_h5ad

# query = '../TabulaSapiensData/PublicDatasets/lung/lung.droplet.h5ad'
# dir = 'lung_atlas'

query = '../TabulaSapiensData/PublicDatasets/lung/droplet_reyfman_normal_lung.20200319.RC4.h5ad' # 25246
query_barcode = '../TabulaSapiensData/PublicDatasets/lung/reyfman.cellid.csv'
dir = 'droplet_reyfman'

# query = '../TabulaSapiensData/PublicDatasets/lung/dropseq_barga_normal_lung.20200319.RC4.h5ad' #5916 #6007
# query_barcode = '../TabulaSapiensData/PublicDatasets/lung/barga.cellid.csv'
# dir = 'dropseq_barga'
#
query = read_h5ad(query)

# get the raw count data rather than the normalized data saved in query.X
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

if scipy.sparse.issparse(csv):
    csv = pd.DataFrame(csv.todense())
else:
    csv = pd.DataFrame(csv)

csv.columns = query.var.index.values
csv.index = query.obs.index.values
csv.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/count.csv'%dir)
labels = query.obs['free_annotation']
labels.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/labels.csv'%dir)

# load atlas data
data_path = '/data/scratch/users/chenling/TabulaSapiens/AnnotationFinal/'
if dir == 'lung_atlas':
    model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/lung_atlas/'
    res_path = 'lung_atlas/'
elif dir == 'droplet_reyfman':
    model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/droplet_reyfman/'
    res_path = 'droplet_reyfman/'
elif dir == 'dropseq_barga':
    model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/dropseq_barga/'
    res_path = 'dropseq_barga/'

if not os.path.exists(model_path):
    os.mkdir(model_path)
if not os.path.exists(model_path+'scanvi/'):
    os.mkdir(model_path+'scanvi/')
if not os.path.exists(res_path):
    os.mkdir(res_path)

adata = data_path + 'CombineCompartments.scANVI.all.count.h5ad'
adata = read_h5ad(adata)
model_file = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/AnnotationFinal.scVI.pkl'
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
lung_data = adata[adata.obs['tissue']=='Lung']

lung_combined = lung_data.concatenate(matched_query)
lung_combined.obs['Tissue Tech'].fillna('Query_10X', inplace=True)
lung_combined.obs['Manual Annotation Round 2'].fillna('unassigned', inplace=True)
tissue_tech_list = list(np.unique(lung_combined.obs['Tissue Tech']))
batch_id = [tissue_tech_list.index(x) for x in lung_combined.obs['Tissue Tech']]
lung_combined.obs['batch'] = batch_id

train_data = AnnDatasetFromAnnData(lung_combined, batch_label='batch')
posterior = get_scvi_posterior(train_data, lung_combined_model_file)
imputed = posterior.sequential().imputation()

csv = lung_combined.X.todense()

if scipy.sparse.issparse(csv):
    csv = pd.DataFrame(csv.todense())
else:
    csv = pd.DataFrame(csv)

csv.columns = lung_combined.var.index.values
csv.index = lung_combined.obs.index.values
csv.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/combined.count.csv'%dir)

csv=imputed
csv = pd.DataFrame(csv)
csv.columns = lung_combined.var.index.values
csv.index = lung_combined.obs.index.values
csv.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/combined.aligned.count.csv'%dir)

labels = lung_combined.obs['scANVI Prediction By Organ']
labels.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/combined.labels.csv'%dir)

batch = lung_combined.obs['Tissue Tech']
batch.to_csv('/data/yosef2/users/chenling/scRNAseq_Benchmark/Snakemake/lung_dataset/%s/combined.batch.csv'%dir)
