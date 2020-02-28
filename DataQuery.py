import sys
query = sys.argv[1]
tissue = sys.argv[2]
retrain = False
sys.path.append('/data/yosef2/users/chenling/scVI/')
from scvi.dataset import AnnDatasetFromAnnData
sys.path.append('/data/yosef2/users/chenling/tabula-sapiens/')
from annotation.utils import *
from anndata import read_h5ad
from sklearn.metrics import confusion_matrix


import obonet
from networkx import ancestors
obo = '../TabulaSapiens/data/cl.obo.txt'
f = open(obo, "r")
co = obonet.read_obo(f)
f.close()

ontology_dict = {}
for x in co.nodes.keys():
    ontology_dict[co.nodes[x]['name']] = x


query = '../TabulaSapiens/lungTravaglini/lung.droplet.h5ad'
# check that the genes in the new dataset is the same as the old one
query = read_h5ad(query)

# load atlas data
data_path = '/data/scratch/users/chenling/TabulaSapiens/AnnotationFinal/'
model_path = '/data/yosef2/users/chenling/TabulaSapiens/scVImodels/'

adata = data_path + 'CombineCompartments.scANVI.all.count.h5ad'
adata = read_h5ad(adata)
labels = adata.obs['Manual Annotation Round 2']
model_file = model_path + 'AnnotationFinal.scVI.pkl'
if retrain:
    # map tissue_tech
    tissue_tech_list = ['Bladder_ss2', 'Blood_ss2', 'Lung_ss2', 'Muscle_ss2',
           'PancreasEndocrine_ss2', 'PancreasExocrine_ss2', 'bladder_10x',
           'blood_10x', 'endopancreas_10x', 'exopancreas1_10x',
           'exopancreas2_10x', 'lung_10x', 'muscle_10x']

    batch_id = [tissue_tech_list.index(x) for x in adata.obs['Tissue Tech']]
    adata.obs['batch'] = batch_id

    data = AnnDatasetFromAnnData(adata, batch_label='batch')

    posterior = get_scvi_posterior(data, model_file)
    latent, _, _ = posterior.get_latent()

    # use input_ann for labels
    data.cell_types, data.labels = np.unique(labels, return_inverse=True)
    data.labels = data.labels.reshape(len(data.labels), 1)
    data.n_labels = len(data.cell_types)
    #
    # generate prediction for all cells at once
    full, pred_celltype = scanvi_pred(data, model_file, model_path+'scanvi/all.manual.semi.20.pkl',nlabels=20, retrain=True)
    x, y = np.unique(pred_celltype, return_counts=True)
    missed = []
    for c in np.unique(data.cell_types):
        if c not in np.unique(pred_celltype):
            print(c)
            missed.append(c)

# filter genes in query to be the same as the atlas data

combined = adata.concatenate(query, join='outer')
genes_combined = combined.var
genes_combined['indices'] = np.arange(combined.shape[1])
shared = pd.concat([adata.var, genes_combined], join='inner', axis=1)
assert np.sum(shared.index==adata.var.index) == adata.shape[1]
combined = combined[:, shared['indices'].values]

query_data = AnnDatasetFromAnnData(combined[combined.obs['batch']=='1',:])
# manually adjust so the dimensions in the model matches
query_data.n_labels = 70
query_data.n_batches = 13
query_data.cell_types = np.unique(labels)
_, query_celltype = scanvi_pred(query_data, model_file,
                                model_path+'scanvi/all.manual.semi.20.pkl',
                                forward_only=True)

np.save(file='query_celltype.npy', arr=query_celltype)

query_celltype = np.load('query_celltype.npy')
labels = pd.read_csv('../TabulaSapiens/lungTravaglini/droplet_normal_lung_blood_P1-3_metadata.csv')
labels = labels['free_annotation'].values
# temp = pd.crosstab(labels, query_celltype)
temp = confusion_matrix(labels, query_celltype)
names = np.sort(np.unique(np.concatenate([labels, query_celltype])))

res = []
for i,x in enumerate(names):
    for j,y in enumerate(names):
        if temp[i,j]!=0:
            res.append([x, y, temp[i,j]])

res = pd.DataFrame(res)
res = res[res[2]>10]
res = res.sort_values(2)
res.to_csv('pretrained_model_performance.csv')
# find CL terms
ontology_map = pd.read_csv('../TabulaSapiens/lungTravaglini/lungatlas_freeannotation_ontology.csv' )
ontology_map = ontology_map.T[1:]
ontology_map = ontology_map.rename(columns={0:'id'})
from annotation.celltype_dict import celltype_dict
query_ontology = []
for x in query_celltype:
    if x in ontology_dict.keys():
        query_ontology.append(ontology_dict[x])
    else:
        query_ontology.append(celltype_dict[x])

lung_ontology_dict = [ ]
ontology_map.index = ontology_map.index.str.strip().str.upper()
labels_ontology = []
missed = []
for x in labels:
    if x.upper() in ontology_map.index.values:
        labels_ontology.append(ontology_map.loc[x.upper(),'id'])
    else:
        missed.append(x)

np.unique(missed)
len(np.unique(labels_ontology))

agreement = []
agreement_by_type = {}
for x in np.unique(labels_ontology):
    agreement_by_type[x] = []
for x,y in zip(query_ontology, labels_ontology):
    if y is not 'None':
        if x == y:
            agreement.append(2)
            agreement_by_type[y].append(2)
        elif x in ancestors(co,y) or y in ancestors(co,x):
                agreement.append(1)
                agreement_by_type[y].append(1)
        else:
            agreement.append(0)

agreement = np.asarray(agreement)
np.mean(agreement>=2)
np.mean(agreement>=1)

# weighted accuracy
for x in agreement_by_type.keys():
    if x != 'None':
        print(co.nodes[x]['name'], ',', np.mean(np.asarray(agreement_by_type[x])==2),',', np.mean(np.asarray(agreement_by_type[x])>=1))
