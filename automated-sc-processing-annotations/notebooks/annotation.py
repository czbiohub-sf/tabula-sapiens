import anndata
import numpy as np
import scanpy as sc
import scvi 
import scanorama
import os
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier


import OnClass as oc
from OnClass.utils import read_data, write_anndata_data
from OnClass.OnClassModel import OnClassModel


def save_results(adata,
              results_adata_path, 
              obs_keys=[],
              obsm_keys=[], 
              compression='gzip'):
  """
  If results_adata_path exists, will load and save results into it
  Else, will save adata to results_adata_path
  
  Parameters
  ----------
  adata
  	adata with results in it
  results_adata_path
  	path to save results. If it already exists, will load and save data to it
  obs_keys
   	obs keys to save
  obsm_keys
  	obsm keys to save
  compression
  	If enabled, will save with compression. Smaller file sizes, but longer save times
  """
  if os.path.exists(results_adata_path):
    results = anndata.read(results_adata_path)
    for key in obs_keys:
      results.obs[key]=adata.obs[key]
    for key in obsm_keys:
      results.obsm[key]=adata.obsm[key]
    results.write(results_adata_path, compression)
  else:
    adata.write(results_adata_path, compression)
  
def process_query(ref_adata,
                  query_adata,
                  ref_labels_key, 
                  ref_batch_key,
                  query_labels_key, 
                  query_batch_key,
                  unknown_celltype_label,
                  join = 'inner',
                  train_from_scratch = True,
                  inplace = True,
                  training_method = 'online' # online or offline
                  ):
  ref_genes = ref_adata.var_names.to_numpy().astype('str')
  query_genes = query_adata.var_names.to_numpy().astype('str')

  if inplace is not True:
    query_adata = query_adata.copy()

  if len(set(query_genes)) != len(query_genes):
    print('Warning: Your genes are not unique.')
  
  if set(ref_genes).issubset(set(query_genes)):
    print('All ref genes are in query dataset. Can use pretrained models')
  else:
    print('Not all reference genese are in query dataset. Retraining models')
    training_method = 'offline'
  
  ref_adata.obs['_labels'] = ref_adata.obs[ref_labels_key]
  ref_adata.obs['_batch'] = ref_adata.obs[ref_batch_key]

  query_adata.obs['_batch'] = query_adata.obs[query_batch_key].astype('str') + '_query'
  
  query_adata.obs['_labels'] = 'unknown'
  known_cell_idx = np.where(query_adata.obs[query_labels_key] != unknown_celltype_label)[0]
  if len(known_cell_idx) != 0:
    query_adata.obs['_labels'][known_cell_idx] = query_adata.obs[query_labels_key][known_cell_idx]

  if training_method == 'online':
    adata = query_adata[:,ref_adata.var_names].copy()
  elif training_method == 'offline':
    adata = ref_adata.concatenate(query_adata, join = join)
    adata.obs_names = [name[:-2] for name in adata.obs_names]
  
  adata.layers["counts"] = adata.X.copy()
  sc.pp.normalize_total(adata, target_sum=1e4)
  sc.pp.log1p(adata)
  sc.pp.scale(adata, max_value=10, zero_center=False)

  sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
  )

  scvi.data.setup_anndata(adata, 
                          batch_key = '_batch', 
                          labels_key = '_labels',
                          layer='counts')
  return adata  

def subsample_dataset(train_data, 
                      labels_key, 
                      n_samples_per_label=100,
                      n_total_samples=None):

    sample_idx = []
    labels, counts = np.unique(train_data.obs[labels_key], return_counts=True)
    print(labels)
    if n_total_samples is not None:
      if n_total_samples > train_data.n_obs:
        print("n_total_samples exceeds size of dataset. Setting to input adata size of {} cells".format(train_data.n_obs))
        return train_data.obs_names
      
      n_samples_per_label = int(n_total_samples/len(labels))
    print("Sampling {} per label".format(n_samples_per_label))
    
    for i, label in enumerate(labels):
        label_locs = np.where(train_data.obs[labels_key] == label)[0]
        if counts[i] < n_samples_per_label:
            sample_idx.append(label_locs)
        else:
            label_subset = np.random.choice(label_locs, n_samples_per_label, replace=False)
            sample_idx.append(label_subset)
    sample_idx = np.concatenate(sample_idx)

    if n_total_samples is not None:
      if len(sample_idx)<n_total_samples:
        all_idx=np.arange(train_data.n_obs)
        remaining_idx = list(set(all_idx)-set(sample_idx))
        remaining_samples=n_total_samples -len(sample_idx)
        label_subset = np.random.choice(remaining_idx, remaining_samples, replace=False)
        sample_idx = np.concatenate((sample_idx,label_subset))
    return train_data.obs_names[sample_idx]

def run_rf_on_hvg(adata,
                  train_idx,
                  test_idx, 
                  save_key = 'rf_on_hvg_pred'):
  train_X = adata[train_idx].layers['counts']
  train_Y = adata[train_idx].obs['_labels'].to_numpy()
  print("Training random forest classifier with {} cells".format(len(train_Y)))
  test_X = adata[test_idx].layers['counts']
  rf = RandomForestClassifier()
  rf.fit(train_X, train_Y)
  rf_pred = rf.predict(test_X)
  adata.obs[save_key] = 'na'
  adata.obs[save_key][test_idx] = rf_pred
 

def run_onclass(adata, 
                ref_adata_path, 
                query_adata_path,
                cl_obo_file, 
                cl_ontology_file,
                cell_ontology_obs_key='cell_ontology_id',
                save_key='onclass_pred',
                n_hidden=500,
                max_iter=20,
                save_model='onclass_model',
                shard_size=50000):
  oc = OnClassModel()
  tp2emb, tp2i, i2tp = oc.EmbedCellTypes(dim=500,
                                         cell_type_network_file=cl_ontology_file,
                                         use_pretrain=None)
  train_X, train_genes, train_Y = read_data(feature_file=ref_adata_path, 
                                            tp2i = tp2i,
                                            AnnData_label=cell_ontology_obs_key)
  oc.train(train_X, 
           train_Y,
           tp2emb,
           train_genes,
           nhidden=[n_hidden],
           max_iter=max_iter,
           use_pretrain = None,
           save_model = save_model)
  test_X, test_genes, test_adata = read_data(feature_file=query_adata_path,
                                               tp2i = tp2i, 
                                               return_AnnData=True)

  test_adata.obs['onclass_pred'] = 'na'
  print(test_adata.n_obs)
  if test_adata.n_obs>shard_size:
    for i in range(0,test_adata.n_obs, shard_size):
      tmp_Anndata = test_adata[i:i+shard_size]
      tmp_X = test_X[i:i+shard_size]
      test_label = oc.predict(tmp_X, test_genes,log_transform=True)
      onclass_pred = write_anndata_data(test_label,
                         tmp_Anndata, 
                         i2tp,
                         name_mapping_file=cl_obo_file)
      onclass_pred = onclass_pred.obs['OnClass_annotation_ontology_name']
      test_adata.obs[save_key][i:i+shard_size] = onclass_pred
  else:
    test_label = oc.predict(test_X, test_genes,log_transform=True)
    onclass_pred = write_anndata_data(test_label,
                                      test_adata, 
                                      i2tp,
                                      name_mapping_file=cl_obo_file)
    onclass_pred = onclass_pred.obs['OnClass_annotation_ontology_name']
    test_adata.obs[save_key] = onclass_pred

  adata.obs['onclass_pred'] = 'na'
  adata.obs['onclass_pred'][test_adata.obs_names]=test_adata.obs['onclass_pred']
  return adata

def run_svm_on_hvg(adata, train_idx, test_idx):
  train_X = adata[train_idx].layers['counts']
  test_X = adata[test_idx].layers['counts']
  print(train_X.shape)
  train_Y = adata[train_idx].obs['_labels'].to_numpy()

  clf = svm.LinearSVC()
  clf.fit(train_X, train_Y)
  svm_pred = clf.predict(test_X)
  
  #save_results
  adata.obs['svm_pred'] = 'na'
  adata.obs['svm_pred'][test_idx] = svm_pred
  return adata

def run_scvi(adata, 
             n_latent = 50,
             n_layers = 3,
             dropout_rate = 0.1,
             dispersion='gene',
             max_epochs = None,
             batch_size = 1024,
             num_workers = 4,
             var_subset_type = 'inner_join',
             save_folder = None,
             overwrite=True
             ):
  model = scvi.model.SCVI(adata, 
                          n_latent = n_latent,
                          n_layers = n_layers,
                          dropout_rate=dropout_rate,
                          dispersion=dispersion,
                         )
  
  model.train(max_epochs = max_epochs, 
              batch_size = batch_size,
              num_workers = num_workers)
  
  adata.obsm['X_scvi'] = model.get_latent_representation()
  if save_folder is not None:
    save_folder = os.path.join(save_folder, 'scvi_model_4batches')
    model.save(save_folder, overwrite=overwrite)
  
def run_knn_on_scvi(adata,
                    train_idx,
                    test_idx, 
                    obsm_key='X_scvi',
                    result_key='knn_on_scvi_pred'):
  
  if obsm_key not in adata.obsm.keys(): 
    print('Please train scVI first or pass in a valid obsm_key.')
  
  train_X = adata[train_idx].obsm[obsm_key]
  test_X = adata[test_idx].obsm[obsm_key]

  train_Y = adata[train_idx].obs['_labels'].to_numpy()
  knn = KNeighborsClassifier(n_neighbors = 15, weights='uniform')
  knn.fit(train_X, train_Y)
  knn_pred= knn.predict(test_X)

  #save_results
  adata.obs[result_key] = 'na'
  adata.obs[result_key][test_idx] = knn_pred

def run_knn_on_scanorama(adata, 
                         train_idx, 
                         test_idx ,
                         obsm_key='X_scanorama',
                         result_key='knn_on_scanorama_pred'):
  
  if obsm_key not in adata.obsm.keys(): 
    print('Please run scanorama first or pass in a valid obsm_key.')

  train_X = adata[train_idx].obsm[obsm_key]
  test_X = adata[test_idx].obsm[obsm_key]

  train_Y = adata[train_idx].obs['_labels'].to_numpy()
  knn = KNeighborsClassifier(n_neighbors = 15, weights='uniform')
  knn.fit(train_X, train_Y)
  knn_pred= knn.predict(test_X)

  #save_results
  adata.obs[result_key] = 'na'
  adata.obs[result_key][test_idx] = knn_pred
  

def run_scanorama(adata, batch_key='_batch'):
  adatas  = [adata[adata.obs[batch_key] == i] for i in np.unique(adata.obs[batch_key])]
  integrated = scanorama.integrate_scanpy(adatas, dimred=50)
  integrated = np.concatenate(integrated)
  tmp_adata = anndata.concat(adatas)
  tmp_adata.obsm['X_scanorama'] = integrated
  
  cell_order = adata.obs_names
  tmp_adata = tmp_adata[cell_order].copy()
  adata.obsm['X_scanorama'] = tmp_adata.obsm['X_scanorama']
  return adata

def test():
  print('tres')

def run_bbknn(adata, batch_key='_batch'):
  print('Running bbknn')
  sc.tl.pca(adata, svd_solver='arpack')
  sc.external.pp.bbknn(adata,
                       batch_key=batch_key,
                       approx=True, 
                       metric='angular',
                       n_pcs=20,
                       trim=None, 
                       n_trees=10,
                       use_faiss=True,
                       set_op_mix_ratio=1.0,
                       local_connectivity=1)
  return adata


def run_knn_on_bbknn(adata,
                     train_idx, 
                     test_idx,
                     labels_key='_labels',
                     result_key='knn_on_bbknn_pred'):
  print('Classifying with knn on bbknn distances')
  distances = adata.obsp['distances']
  
  #change this later
  #will fail if ref_n_obs is not sequential
  ref_n_obs = len(train_idx)

  train_Y = adata[train_idx].obs[labels_key].to_numpy()
  train_distances = distances[:ref_n_obs, :ref_n_obs]

  knn = KNeighborsClassifier(n_neighbors = 2, metric='precomputed')
  knn.fit(train_distances, y = train_Y)

  test_distances = distances[ref_n_obs:, :ref_n_obs]
  knn_pred = knn.predict(test_distances)

  #save_results
  adata.obs[result_key] = 'na'
  adata.obs[result_key][test_idx] = knn_pred
  
def run_scanvi(adata, 
               n_layers=3,
               dropout_rate=0.2,
               n_classifier_layers=1,
               classifier_dropout=0.4, 
               n_epochs_unsupervised=0,
               n_epochs_semisupervised=None,
               n_latent=100,
               batch_size=1024,
               num_workers=4,
               n_epochs_kl_warmup=20,
               n_samples_per_label=100,
               save_folder=None,
               overwrite=True):
  if n_epochs_semisupervised is None:
    n_epochs_semisupervised = np.min(
      [round((20000 / adata.shape[0]) * 400), 400]
    )
  arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",  
  )
  scanvae_model_kwargs = arches_params.copy()
  scanvae_model_kwargs.update({'classifier_parameters':{'n_layers':n_classifier_layers, 
                                                        "dropout_rate":classifier_dropout}})

  model = scvi.model.SCANVI(adata,
                          unlabeled_category = 'unknown', 
                          n_layers=n_layers, 
                          encode_covariates=True,
                          dropout_rate=dropout_rate,
                          n_latent=n_latent, 
                          vae_model_kwargs=arches_params,
                          scanvae_model_kwargs=scanvae_model_kwargs
                          )
  model.train(n_epochs_unsupervised=n_epochs_unsupervised, 
            n_epochs_semisupervised=n_epochs_semisupervised,
            n_epochs_kl_warmup=n_epochs_kl_warmup, 
            batch_size=batch_size, 
            num_workers=num_workers, 
            n_samples_per_label=n_samples_per_label
            )
  adata.obsm['X_scanvi'] = model.get_latent_representation()
  adata.obs['scanvi_pred'] = model.predict()

  if save_folder is not None:
    save_folder = os.path.join(save_folder, 'scanvi_model')
    model.save(save_folder, overwrite=overwrite)
  
